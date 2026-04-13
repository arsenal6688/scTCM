import pandas as pd
import pickle
from pathlib import Path
from typing import Dict, List, Union

# ====================== 步骤1：CID匹配HERB + 参数化筛选 ======================
def match_and_filter_herb(
    template_path: Union[str, Path],
    herb_info_path: Union[str, Path],
    drug_likeness_min: float = 0.1,
    ob_score_min: float = 30.0,
    sep: str = "\t"
) -> pd.DataFrame:
    print("===== 步骤1：HERB数据库匹配与成分筛选 =====")


    template_df = pd.read_csv(template_path, sep=sep, na_values=["NA", "na", "NaN", "nan"])
    template_df = template_df.drop_duplicates(subset=["PubChem_CID"], keep="first")
    total_original = len(template_df)
    print(f"原始输入成分总数：{total_original}")

    herb_df = pd.read_csv(herb_info_path, sep=sep, na_values=["NA", "na", "NaN", "nan"])
    herb_df.rename(columns={"PubChem_id": "PubChem_CID"}, inplace=True)
    herb_df = herb_df.drop_duplicates(subset=["PubChem_CID"], keep="first")

    herb_df["Drug_likeness"] = pd.to_numeric(herb_df["Drug_likeness"], errors="coerce")
    herb_df["OB_score"] = pd.to_numeric(herb_df["OB_score"], errors="coerce")

    merge_df = template_df.merge(herb_df[["PubChem_CID", "Drug_likeness", "OB_score"]],
                                 on="PubChem_CID", how="left")

    matched = merge_df["Drug_likeness"].notna().sum()
    unmatched = total_original - matched
    print(f"在HERB库中匹配成功：{matched} 个")

    if unmatched > 0:
        unmatched_cids = merge_df.loc[merge_df["Drug_likeness"].isna(), "PubChem_CID"].tolist()
        print(f"在HERB库中无记录：{unmatched} 个 → 标记【无法筛选】")
        print(f"无法筛选的CID：{unmatched_cids}")

    merge_df["filter_status"] = "无法筛选"
    pass_mask = (merge_df["Drug_likeness"].notna()) & \
                (merge_df["Drug_likeness"] >= drug_likeness_min) & \
                (merge_df["OB_score"] >= ob_score_min)
    merge_df.loc[pass_mask, "filter_status"] = "筛选通过"
    pass_num = pass_mask.sum()
    
    print(f"\n过滤条件：Drug_likeness ≥ {drug_likeness_min} 且 OB_score ≥ {ob_score_min}")
    print(f"满足筛选条件通过：{pass_num} 个")
    print(f"未通过筛选（含无HERB记录）：{total_original - pass_num} 个\n")

    return merge_df

# ====================== 步骤2：靶点提取 ======================
def extract_targets(
    filtered_df: pd.DataFrame,
    known_target_path: Union[str, Path],
    use_filtered_only: bool = False,
    sep: str = "\t"
) -> Dict[str, List[str]]:
    print("===== 步骤2：成分-靶点映射提取 =====")

    if use_filtered_only:
        df_use = filtered_df[filtered_df["filter_status"] == "筛选通过"].copy()
        print(f"📌 当前模式：仅对【筛选通过】的成分提取靶点 | 总数：{len(df_use)} 个")
    else:
        df_use = filtered_df.copy()
        print(f"📌 当前模式：对【全部原始成分】提取靶点 | 总数：{len(df_use)} 个")

    known_df = pd.read_csv(known_target_path, sep=sep, compression="infer", na_values=["NA", "na", "NaN", "nan"])
    known_df = known_df.drop_duplicates(subset=["PubChem_CID"], keep="first")

    result_dict = {}
    no_target_cids = []

    for _, row in df_use.iterrows():
        cid = row["PubChem_CID"]
        custom = row["Self_defined_targets"]

        # --------------- 过滤掉 NA ----------------
        if pd.notna(custom) and str(custom).strip() != "":
            targets = [t.strip() for t in str(custom).split("|") if t.strip() and t.strip() != "NA"]
            result_dict[cid] = targets
            continue

        hit = known_df[known_df["PubChem_CID"] == cid]
        if not hit.empty:
            ts = hit["known_target_proteins"].iloc[0]
            if pd.notna(ts) and str(ts).strip() != "":
                # --------------- 过滤掉 NA ----------------
                targets = [t.strip() for t in str(ts).split("|") if t.strip() and t.strip() != "NA"]
                result_dict[cid] = targets
                continue

        no_target_cids.append(cid)

    print(f"成功获取靶点：{len(result_dict)} 个")
    print(f"未找到任何靶点：{len(no_target_cids)} 个")

    if len(no_target_cids) > 0:
        print(f"无靶点CID列表：{no_target_cids}\n")

    return result_dict

# ====================== 主函数 ======================
def run_ingredient2target(
    template_path: Union[str, Path],
    herb_info_path: Union[str, Path],
    known_target_path: Union[str, Path],
    output_pkl: str = "ingredient2target_result.pkl",
    drug_likeness_min: float = 0.1,
    ob_score_min: float = 30.0,
    use_filtered_only: bool = False
) -> Dict[str, List[str]]:

    filtered_df = match_and_filter_herb(
        template_path=template_path,
        herb_info_path=herb_info_path,
        drug_likeness_min=drug_likeness_min,
        ob_score_min=ob_score_min
    )

    cid2target = extract_targets(
        filtered_df=filtered_df,
        known_target_path=known_target_path,
        use_filtered_only=use_filtered_only
    )

    with open(output_pkl, "wb") as f:
        pickle.dump(cid2target, f)

    print("===== 任务完成 =====")
    print(f"结果已保存至：{output_pkl}")
    print(f"最终成分-靶点映射总数：{len(cid2target)}\n")

    return cid2target

# ====================== 直接运行 ======================
if __name__ == "__main__":
    run_ingredient2target(
        template_path="tcm_ingredient2target-template.txt",
        herb_info_path="HERB_ingredient_info_v2_tmp.txt",
        known_target_path="known_browse_by_ingredients_head100.txt.gz",
        output_pkl="ingredient2target.pkl",
        drug_likeness_min=0.1,
        ob_score_min=30.0,
        use_filtered_only=False
    )
