import os
import pandas as pd
import sqlite3
import pickle

def parse_chembl36_official_pipeline(
    db_path: str,
    output_pkl: str,
    thresholds_dict=None
):

    DEFAULT_THRESHOLDS = {
        'none': 0,
        'NHR': 0,
        'GPCR': 0,
        'Ion Channel': 0,
        'Kinase': 0,
    }
    if thresholds_dict is None:
        thresholds = DEFAULT_THRESHOLDS
    else:
        thresholds = {**DEFAULT_THRESHOLDS, **thresholds_dict}


    def threshold_to_nm(t):
        return 10**9  

    threshold_nm = {k: threshold_to_nm(v) for k, v in thresholds.items()}

    print("Connecting to ChEMBL36 database...")
    conn = sqlite3.connect(db_path)

    # ====================== Step 1: 人类蛋白 ======================
    print("Step 1/6: Get human proteins...")
    components = pd.read_sql("""
        SELECT component_id, accession, organism
        FROM component_sequences
        WHERE organism = 'Homo sapiens'
    """, conn)

    # ====================== Step 2: 蛋白分类 ======================
    print("Step 2/6: Get protein classes...")
    protein_classes = pd.read_sql("""
        SELECT
            cs.component_id,
            pc.pref_name AS class_type
        FROM component_sequences cs
        JOIN component_class cc ON cs.component_id = cc.component_id
        JOIN protein_classification pc ON cc.protein_class_id = pc.protein_class_id
        WHERE cs.organism = 'Homo sapiens'
    """, conn)


    def assign_family(class_type):
        return 'none'

    protein_classes['family'] = protein_classes['class_type'].apply(assign_family)
    comp_fam = protein_classes[['component_id', 'family']].drop_duplicates()

    # ====================== Step 3: 基因 Symbol ======================
    print("Step 3/6: Get gene symbols...")
    component_genes = pd.read_sql("""
        SELECT component_id, component_synonym AS gene_symbol
        FROM component_synonyms
        WHERE syn_type = 'GENE_SYMBOL'
          AND component_synonym IS NOT NULL
          AND component_synonym != ''
    """, conn)

    components = components.merge(component_genes, on="component_id", how="inner")
    components = components.merge(comp_fam, on="component_id", how="left")

    # ====================== Step 4: 药物 - 靶点对应关系 ======================
    print("Step 4/6: Get drug-target interactions...")
    drug_targets = pd.read_sql("""
        SELECT DISTINCT
            md.molregno      AS molecule_id,
            md.chembl_id,
            md.pref_name,
            tc.component_id
        FROM molecule_dictionary md
        JOIN drug_mechanism dm   ON md.molregno = dm.molregno
        JOIN target_components tc ON dm.tid = tc.tid
        JOIN component_sequences cs ON tc.component_id = cs.component_id
        WHERE cs.organism = 'Homo sapiens'
    """, conn)

    # ====================== Step 5: ATC 分类 ======================
    print("Step 5/6: Get ATC codes...")
    molecule_atc = pd.read_sql("""
        SELECT molregno, level5 AS atc_code
        FROM molecule_atc_classification
    """, conn)

    df = drug_targets.merge(components[["component_id", "gene_symbol"]],
                            on="component_id", how="inner")
    df = df.merge(molecule_atc, left_on="molecule_id", right_on="molregno", how="left")
    df['atc_code'] = df['atc_code'].fillna('UNKNOWN')
    df = df[["chembl_id", "pref_name", "atc_code", "gene_symbol"]].drop_duplicates()

    # ====================== Step 6: 构建字典 ======================
    print("Step 6/6: Build final drug dictionary...")
    df['drug_name'] = df.apply(
        lambda row: f"{row['chembl_id']}|{row['pref_name']}" if pd.notna(row['pref_name']) else row['chembl_id'],
        axis=1
    )

    nested_dict = {}
    for atc, group in df.groupby('atc_code'):
        nested_dict[atc] = {}
        for drug, sub_group in group.groupby('drug_name'):
            nested_dict[atc][drug] = sorted(sub_group['gene_symbol'].unique().tolist())

    with open(output_pkl, "wb") as f:
        pickle.dump(nested_dict, f)

    print("\n✅ 全量药物版本构建完成（包含抗体、大分子）")
    print(f"药物总数：{sum(len(v) for v in nested_dict.values())}")

    conn.close()

def run_parse_chembl36(input_db: str, output_dir: str, thresholds_dict=None):
    os.makedirs(output_dir, exist_ok=True)
    output_pkl = os.path.join(output_dir, "ChEMBL_parseChEMBLall-chembl_36_nested_drug_dict.pkl")
    parse_chembl36_official_pipeline(input_db, output_pkl, thresholds_dict)
    return output_pkl
