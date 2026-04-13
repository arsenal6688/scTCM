import pandas as pd
import pickle
import tempfile
from pathlib import Path
from typing import Dict, List, Union, Optional
from .ingredient2target import match_and_filter_herb
from .ingredient2target import extract_targets

def parse_ingredients(ingredient_str: Optional[str]) -> List[str]:
    cids = []
    if pd.isna(ingredient_str):
        return cids
    for part in str(ingredient_str).split("|"):
        part = part.strip()
        if "(" in part and ")" in part:
            cid = part[part.find("(")+1:part.find(")")].strip()
            if cid.isdigit():
                cids.append(cid)
    return cids

def run_herb2target(
    template_path: Union[str, Path],
    herb_browse_path: Union[str, Path],
    herb_info_path: Union[str, Path],
    known_target_path: Union[str, Path],
    output_pkl: Union[str, Path] = "herb2target_result.pkl",
    drug_likeness_min: float = 0.1,
    ob_score_min: float = 10.0,
    use_filtered_only: bool = False,
    sep: str = "\t"
) -> Dict[str, Dict[str, List[str]]]:

    print("=" * 70)
    print("📌 开始中药-靶点映射流程（herb2target）")
    print("=" * 70)

    print("\n【步骤1：读取输入文件】")
    template_df = pd.read_csv(template_path, sep=sep, na_values=["NA", "na", "NaN", "nan"])
    template_df = template_df.drop_duplicates("Name_pinyin")
    print(f"模板文件中药总数：{len(template_df)}")
    print(f"中药拼音名列表：{template_df['Name_pinyin'].tolist()}")

    browse_df = pd.read_csv(herb_browse_path, sep=sep, na_values=["NA", "na", "NaN", "nan"])
    browse_df = browse_df.drop_duplicates("Pinyin.Name")
    print(f"herb_browse文件中药总数：{len(browse_df)}")

    print("\n【步骤2：中药与成分匹配】")
    df = template_df.merge(
        browse_df[["Pinyin.Name", "Ingredients"]],
        left_on="Name_pinyin",
        right_on="Pinyin.Name",
        how="left"
    )
    matched = df["Ingredients"].notna().sum()
    print(f"成功匹配成分的中药数：{matched} 个")
    print(f"未匹配成分的中药数：{len(template_df)-matched} 个")

    print("\n【步骤3：提取成分PubChem_CID】")
    herb_cid_map = {}
    for _, row in df.iterrows():
        name = row["Name_pinyin"]
        self_add = row["Self_added_ingredients"]
        ingredients = row["Ingredients"]

        if pd.notna(self_add) and str(self_add).strip() != "":
            cids = [c.strip() for c in str(self_add).split("|") if c.strip().isdigit()]
            print(f"中药[{name}]：使用自定义CID → {cids[:5]}...")
        else:
            cids = parse_ingredients(ingredients)
            print(f"中药[{name}]：从browse提取CID → {cids[:5]}...")

        herb_cid_map[name] = list(set(cids))

    total_cids = sum(len(cs) for cs in herb_cid_map.values())
    print(f"\n提取结果：")
    print(f"有有效成分的中药数：{len(herb_cid_map)} 个")
    print(f"提取的总成分CID数：{total_cids} 个")

    print("\n" + "=" * 70)
    print("【步骤4：成分筛选 + 靶点提取】")
    print("=" * 70)

    temp_rows = []
    for herb, cids in herb_cid_map.items():
        for cid in cids:
            temp_rows.append({
                "PubChem_CID": cid,
                "Name_IUPAC": f"Herb_{herb}_{cid}",
                "Self_defined_targets": ""
            })

    temp_df = pd.DataFrame(temp_rows)
    print(f"临时成分模板总数：{len(temp_df)} 个")

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
        temp_df.to_csv(f, sep=sep, index=False, na_rep="NA")
        temp_file = f.name

    filtered_df = match_and_filter_herb(
        template_path=temp_file,
        herb_info_path=herb_info_path,
        drug_likeness_min=drug_likeness_min,
        ob_score_min=ob_score_min,
        sep=sep
    )
    Path(temp_file).unlink()

    cid2target = extract_targets(
        filtered_df=filtered_df,
        known_target_path=known_target_path,
        use_filtered_only=use_filtered_only,
        sep=sep
    )

    # 统一CID为字符串
    cid2target_str = {str(k): v for k, v in cid2target.items()}

    print("\n【步骤5：构建中药-靶点最终映射】")
    result = {}
    for herb, cids in herb_cid_map.items():
        herb_dict = {}
        for cid in cids:
            ts = cid2target_str.get(str(cid), [])
            herb_dict[cid] = ts
        
        # ====================== 只保留有靶点的成分 ======================
        herb_dict_filtered = {cid: ts for cid, ts in herb_dict.items() if len(ts) > 0}
        result[herb] = herb_dict_filtered

    print("\n" + "="*70)
    print("📊 中药-靶点映射统计结果")
    print("="*70)
    print(f"总中药数：{len(result)} 个")

    for herb, cid_dict in result.items():
        total_c = len(cid_dict)
        has_t = sum(1 for t in cid_dict.values() if len(t) > 0)
        all_ts = list({t for ts in cid_dict.values() for t in ts})

        print(f"\n中药[{herb}]：")
        print(f"  - 成分CID数：{total_c} 个")
        print(f"  - 有靶点的成分数：{has_t} 个")
        print(f"  - 总靶点数：{sum(len(ts) for ts in cid_dict.values())} 个")
        if all_ts:
            print(f"  - 去重后靶点列表：{all_ts[:10]}...")
        else:
            print(f"  - 去重后靶点列表：[]")

    with open(output_pkl, "wb") as f:
        pickle.dump(result, f)

    print("\n" + "="*70)
    print(f"✅ 中药-靶点映射完成！")
    print(f"📁 结果已保存至：{output_pkl}")
    print("="*70)

    return result


if __name__ == "__main__":
    run_herb2target(
        template_path="tcm_herb2target-template.txt",
        herb_browse_path="herb_browse.txt",
        herb_info_path="HERB_ingredient_info_v2_tmp.txt",
        known_target_path="known_browse_by_ingredients.txt.gz",
        drug_likeness_min=0.1,
        ob_score_min=10.0,
        use_filtered_only=False
    )
