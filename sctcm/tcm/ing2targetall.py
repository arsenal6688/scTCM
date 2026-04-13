import pandas as pd
import pickle
import os

def run_ing2targetall(herb_info_path, known_target_path, output_dir, ob_min=0.0, drug_likeness_min=0.0):
    os.makedirs(output_dir, exist_ok=True)

    # 读取靶点库
    known_df = pd.read_csv(known_target_path, sep="\t", compression="infer")
    known_df = known_df.dropna(subset=["PubChem_CID"])
    known_df["cid"] = known_df["PubChem_CID"].astype(str).str.strip()
    known_df["cid"] = known_df["cid"].str.replace(".0", "", regex=False)
    known_df = known_df[known_df["cid"].str.isdigit()]

    template_df = known_df[["cid"]].drop_duplicates()
    print(f"从靶点库提取成分总数：{len(template_df)}")

    # 读取 HERB
    herb_df = pd.read_csv(herb_info_path, sep="\t")
    herb_df = herb_df.dropna(subset=["PubChem_id"])
    herb_df["cid"] = herb_df["PubChem_id"].astype(str).str.strip()
    herb_df["cid"] = herb_df["cid"].str.replace(".0", "", regex=False)
    herb_df = herb_df[herb_df["cid"].str.isdigit()]

    herb_df["OB_score"] = pd.to_numeric(herb_df["OB_score"], errors="coerce")
    herb_df["Drug_likeness"] = pd.to_numeric(herb_df["Drug_likeness"], errors="coerce")
    herb_df = herb_df.dropna(subset=["OB_score", "Drug_likeness"])

    # 匹配
    merge_df = template_df.merge(herb_df[["cid", "OB_score", "Drug_likeness"]], on="cid", how="inner")
    print(f"匹配成功数量：{len(merge_df)}")

    # 过滤
    keep_df = merge_df[(merge_df["OB_score"] >= ob_min) & (merge_df["Drug_likeness"] >= drug_likeness_min)].copy()
    print(f"满足条件数量：{len(keep_df)}")

    # 提取靶点
    cid2targets = {}
    for _, row in keep_df.iterrows():
        cid = row["cid"]
        hit = known_df[known_df["cid"] == cid]
        if hit.empty:
            continue
        ts = hit["known_target_proteins"].iloc[0]
        if pd.isna(ts):
            continue
        targets = [t.strip() for t in str(ts).split("|") if t.strip()]
        if targets:
            cid2targets[cid] = sorted(list(set(targets)))

    print(f"最终可输出成分数：{len(cid2targets)}")

    out = os.path.join(output_dir, "tcm_ing2targetall.pkl")
    with open(out, "wb") as f:
        pickle.dump(cid2targets, f)
    print(f"文件已保存：{out}")
    return cid2targets
