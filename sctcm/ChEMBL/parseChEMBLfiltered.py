import os
import pandas as pd
import sqlite3
import pickle

def parse_chembl36_official_pipeline(
    db_path: str,
    output_pkl: str,
    thresholds_dict=None
):
    # ===================== 默认阈值 =====================
    DEFAULT_THRESHOLDS = {
        'none': 6,          # 1 μM
        'NHR': 7,           # 100 nM
        'GPCR': 7,          # 100 nM
        'Ion Channel': 5,   # 10 μM
        'Kinase': 7.53,     # 30 nM
    }
    if thresholds_dict is None:
        thresholds = DEFAULT_THRESHOLDS
    else:
        thresholds = {**DEFAULT_THRESHOLDS, **thresholds_dict}


    def threshold_to_nm(t):
        return 10 ** (9 - t)

    threshold_nm = {k: threshold_to_nm(v) for k, v in thresholds.items()}

    print("Connecting to ChEMBL36 database...")
    conn = sqlite3.connect(db_path)

    # ====================== Step 1: 人类蛋白组件 ======================
    print("Step 1/7: Get human protein components...")
    components = pd.read_sql("""
        SELECT component_id, accession, organism
        FROM component_sequences
        WHERE organism = 'Homo sapiens'
    """, conn)

    # ====================== Step 2: 蛋白分类 ======================
    print("Step 2/7: Filter allowed target classes...")
    protein_classes = pd.read_sql("""
        SELECT
            cs.component_id,
            pc.pref_name AS class_type
        FROM component_sequences cs
        JOIN component_class cc ON cs.component_id = cc.component_id
        JOIN protein_classification pc ON cc.protein_class_id = pc.protein_class_id
        WHERE cs.organism = 'Homo sapiens'
    """, conn)

    allowed_classes = {
        'kinase', 'transferase', 'phosphatase', 'protease', 'enzyme',
        'transcription factor', 'epigenetic regulator', 'transporter',
        'channel', 'membrane receptor', 'nuclear receptor', 'secreted protein'
    }

    def is_allowed(class_type):
        return any(c in allowed_classes for c in str(class_type).lower().split())

    protein_classes = protein_classes[protein_classes['class_type'].apply(is_allowed)]
    components = components[components['component_id'].isin(protein_classes['component_id'])]

    # 标注蛋白家族
    def assign_family(class_type):
        s = str(class_type).lower()
        if 'nuclear receptor' in s or 'nhr' in s:
            return 'NHR'
        elif 'gpcr' in s or 'membrane receptor' in s:
            return 'GPCR'
        elif 'ion channel' in s or 'channel' in s:
            return 'Ion Channel'
        elif 'kinase' in s:
            return 'Kinase'
        else:
            return 'none'

    protein_classes['family'] = protein_classes['class_type'].apply(assign_family)
    comp_fam = protein_classes[['component_id', 'family']].drop_duplicates()

    # ====================== Step 3: 基因 Symbol + ATC 分类 ======================
    print("Step 3/7: Map gene symbols & ATC codes...")
    
    # 基因symbol
    component_genes = pd.read_sql("""
        SELECT component_id, component_synonym AS gene_symbol
        FROM component_synonyms
        WHERE syn_type = 'GENE_SYMBOL'
          AND component_synonym IS NOT NULL
          AND component_synonym != ''
    """, conn)
    

    molecule_atc = pd.read_sql("""
        SELECT 
            molregno, 
            level5 as atc_code
        FROM molecule_atc_classification
    """, conn)

    components = components.merge(component_genes, on="component_id", how="inner")
    components = components.merge(comp_fam, on="component_id", how="left")
    components['family'] = components['family'].fillna('none')

    # ====================== Step 4: 活性数据 ======================
    print("Step 4/7: Filter activities (EC50/IC50/Ki/Kd)...")
    activities = pd.read_sql("""
        SELECT
            act.activity_id,
            act.assay_id,
            act.molregno AS molecule_id,
            act.standard_value AS value,
            td.tid
        FROM activities act
        JOIN assays ass ON act.assay_id = ass.assay_id
        JOIN target_dictionary td ON ass.tid = td.tid
        WHERE act.standard_type IN ('EC50', 'IC50', 'Ki', 'Kd')
          AND act.standard_units = 'nM'
          AND act.standard_relation IN ('=', '<', '<=')
          AND act.standard_value > 0
          AND td.organism = 'Homo sapiens'
    """, conn)
    activities['value'] = activities['value'].astype(float)

    # ====================== Step 5: Target → Component ======================
    print("Step 5/7: Link targets to components...")
    target_components = pd.read_sql("SELECT tid, component_id FROM target_components", conn)
    activities = activities.merge(target_components, on="tid", how="inner")
    activities = activities[activities['component_id'].isin(components['component_id'])]

    # 按家族阈值过滤
    activities = activities.merge(components[['component_id', 'family']], on="component_id")
    activities['threshold_nm'] = activities['family'].map(threshold_nm)
    activities = activities[activities['value'] <= activities['threshold_nm']]

    # ====================== Step 6: 药物信息 + ATC ======================
    print("Step 6/7: Merge drug information...")
    molecules = pd.read_sql("""
        SELECT 
            molregno, 
            chembl_id, 
            pref_name 
        FROM molecule_dictionary
    """, conn)
    molecules.rename(columns={"molregno": "molecule_id"}, inplace=True)
    activities = activities.merge(molecules, on="molecule_id", how="inner")
    activities = activities.merge(molecule_atc, left_on="molecule_id", right_on="molregno", how="left")
    

    activities['atc_code'] = activities['atc_code'].fillna('UNKNOWN')

    # ====================== Step 7: 生成 嵌套字典 ======================
    print("Step 7/7: Build NESTED drug dictionary (drug2cell compatible)...")
    
    df = activities.merge(
        components[["component_id", "gene_symbol"]],
        on="component_id",
        how="inner"
    )
    df = df[["chembl_id", "pref_name", "atc_code", "gene_symbol"]].drop_duplicates()

    
    df['drug_name'] = df.apply(
        lambda row: f"{row['chembl_id']}|{row['pref_name']}" if pd.notna(row['pref_name']) else row['chembl_id'],
        axis=1
    )

    # ==============================================
    # 核心：生成嵌套字典 ATC → Drug → Targets
    # ==============================================
    nested_dict = {}
    for atc, group in df.groupby('atc_code'):
        nested_dict[atc] = {}
        for drug, sub_group in group.groupby('drug_name'):
            nested_dict[atc][drug] = sorted(sub_group['gene_symbol'].unique().tolist())

    # 保存
    with open(output_pkl, "wb") as f:
        pickle.dump(nested_dict, f)

    # ====================== 统计输出 ======================
    print("\n" + "="*60)
    print("📊 最终统计（嵌套字典版本 · drug2cell 100% 兼容）")
    print("="*60)
    for fam in thresholds:
        print(f"  {fam:12s} | pChEMBL={thresholds[fam]:.2f} | nM={threshold_nm[fam]:>8.1f}")
    print("-"*60)
    print(f"ATC 分类总数: {len(nested_dict)}")
    print(f"药物总数: {sum(len(d) for d in nested_dict.values())}")
    print(f"输出文件: {output_pkl}")
    print("="*60)
    print("✅ 运行完成！")

    conn.close()

def run_parse_chembl36(input_db: str, output_dir: str, thresholds_dict=None):
    os.makedirs(output_dir, exist_ok=True)
    output_pkl = os.path.join(output_dir, "ChEMBL_parseChEMBLfiltered-chembl_36_nested_drug_dict.pkl")
    parse_chembl36_official_pipeline(input_db, output_pkl, thresholds_dict)
    return output_pkl
