import pandas as pd
import pickle
from pathlib import Path
from typing import Dict, List, Union, Optional
from .herb2target import parse_ingredients  
from .herb2target import run_herb2target  
from .ingredient2target import match_and_filter_herb, extract_targets  

def run_formula2target(
    template_path: Union[str, Path],
    formula_browse_path: Union[str, Path],
    herb_browse_path: Union[str, Path],
    herb_info_path: Union[str, Path],
    known_target_path: Union[str, Path],
    output_pkl: Union[str, Path] = "formula2target_result.pkl",
    drug_likeness_min: float = 0.1,
    ob_score_min: float = 10.0,
    use_filtered_only: bool = False,
    sep: str = "\t"
) -> Dict[str, Dict[str, Dict[str, List[str]]]]:
    """
    方剂-靶点映射主函数：复用已有核心逻辑，实现 方剂→中药→成分→靶点 的完整流程
    
    参数
    -----
    template_path: 方剂模板文件 tcm_formula2target-template.txt
    formula_browse_path: 方剂-中药映射文件 formula_browse_100.txt
    herb_browse_path: 中药-成分映射文件 herb_browse.txt
    herb_info_path: HERB成分信息文件 HERB_ingredient_info_v2.txt
    known_target_path: 已知靶点文件 known_browse_by_ingredients.txt(.gz)
    output_pkl: 输出pkl文件路径
    drug_likeness_min: 类药性最小阈值（复用）
    ob_score_min: 口服生物利用度最小阈值（复用）
    use_filtered_only: 仅对筛选通过的成分提取靶点（复用）
    sep: 文件分隔符
    """
    print("=" * 80)
    print("📌 开始方剂-靶点映射流程（formula2target）")
    print("=" * 80)

    # ====================== 步骤1：读取并校验输入文件 ======================
    print("\n【步骤1：读取输入文件】")
    # 读取方剂模板
    template_df = pd.read_csv(template_path, sep=sep, na_values=["NA", "na", "NaN", "nan", ""])
    required_template_cols = ["Name_pinyin", "Composition_pinyin", "Self_added_herbs"]
    for col in required_template_cols:
        if col not in template_df.columns:
            raise ValueError(f"方剂模板文件缺少必需列：{col}")
    template_df = template_df.drop_duplicates(subset=["Name_pinyin"], keep="first")
    print(f"方剂模板总数：{len(template_df)} 个")
    print(f"方剂拼音名列表：{template_df['Name_pinyin'].tolist()}")

    # 读取方剂-中药映射文件
    formula_browse_df = pd.read_csv(formula_browse_path, sep=sep, na_values=["NA", "na", "NaN", "nan", ""])
    required_formula_cols = ["Pinyin.Name", "Pinyin.composition"]
    for col in required_formula_cols:
        if col not in formula_browse_df.columns:
            raise ValueError(f"formula_browse文件缺少必需列：{col}")
    formula_browse_df = formula_browse_df.drop_duplicates(subset=["Pinyin.Name"], keep="first")
    print(f"formula_browse文件方剂总数：{len(formula_browse_df)} 个")

    # ====================== 步骤2：匹配方剂与中药（Name_pinyin ↔ Pinyin.Name） ======================
    print("\n【步骤2：方剂与中药匹配】")
    # 合并模板与formula_browse文件（按拼音名匹配）
    merge_df = template_df.merge(
        formula_browse_df[["Pinyin.Name", "Pinyin.composition"]],
        left_on="Name_pinyin",
        right_on="Pinyin.Name",
        how="left"
    )
    merge_df = merge_df.drop(columns=["Pinyin.Name"])  # 删除重复列

    # 统计匹配情况
    matched_formula = merge_df["Pinyin.composition"].notna().sum()
    unmatched_formula = len(merge_df) - matched_formula
    print(f"成功匹配中药的方剂数：{matched_formula} 个")
    print(f"未匹配中药的方剂数：{unmatched_formula} 个")
    if unmatched_formula > 0:
        unmatched_names = merge_df.loc[merge_df["Pinyin.composition"].isna(), "Name_pinyin"].tolist()
        print(f"未匹配中药的方剂：{unmatched_names}")

    # ====================== 步骤3：提取方剂包含的中药（区分自定义/自动提取） ======================
    print("\n【步骤3：提取方剂包含的中药】")
    formula_herb_map = {}  # 存储每个方剂的中药列表：{方剂拼音: [中药拼音1, 中药拼音2, ...]}
    
    for _, row in merge_df.iterrows():
        formula_name = row["Name_pinyin"]
        self_added_herbs = row["Self_added_herbs"]
        browse_composition = row["Pinyin.composition"]
        template_composition = row["Composition_pinyin"]
        
        herbs = []
        # 情况1：自定义中药不为NA → 按 | 分割提取中药拼音
        if pd.notna(self_added_herbs) and str(self_added_herbs).strip() != "":
            self_herbs = [herb.strip() for herb in str(self_added_herbs).split(",") if herb.strip()]
            herbs.extend(self_herbs)
            print(f"方剂[{formula_name}]：使用自定义中药 → {self_herbs}")
        
        # 情况2：自定义中药为NA → 优先使用模板的Composition_pinyin，无则用browse文件
        else:
            if pd.notna(template_composition) and str(template_composition).strip() != "":
                # 从模板提取中药（按逗号分割）
                template_herbs = [herb.strip() for herb in str(template_composition).split(",") if herb.strip()]
                herbs.extend(template_herbs)
                print(f"方剂[{formula_name}]：从模板提取中药 → {template_herbs}")
            elif pd.notna(browse_composition) and str(browse_composition).strip() != "":
                # 从browse文件提取中药（按逗号分割）
                browse_herbs = [herb.strip() for herb in str(browse_composition).split(",") if herb.strip()]
                herbs.extend(browse_herbs)
                print(f"方剂[{formula_name}]：从formula_browse提取中药 → {browse_herbs[:10]}...")  # 超长列表截断
        
        # 去重并存储（过滤空字符串）
        herbs = list(set([h for h in herbs if h]))
        formula_herb_map[formula_name] = herbs

    # 统计中药提取情况
    total_herbs = sum(len(herbs) for herbs in formula_herb_map.values())
    formula_with_herb = sum(1 for herbs in formula_herb_map.values() if len(herbs) > 0)
    print(f"\n提取结果：")
    print(f"有有效中药的方剂数：{formula_with_herb} 个")
    print(f"提取的总中药数（去重前）：{total_herbs} 个")
    all_herbs = list(set([h for herbs in formula_herb_map.values() for h in herbs]))
    print(f"提取的唯一中药数：{len(all_herbs)} 个")
    print(f"唯一中药列表：{all_herbs[:10]}...")

    # ====================== 步骤4：复用herb2target，批量获取中药-靶点映射 ======================
    print("\n" + "=" * 80)
    print("【步骤4：批量获取中药-靶点映射（复用herb2target逻辑）】")
    print("=" * 80)

    # 构建临时中药模板（适配herb2target的输入格式）
    temp_herb_template_data = []
    for formula_name, herbs in formula_herb_map.items():
        for herb in herbs:
            temp_herb_template_data.append({
                "Name_pinyin": herb,
                "Self_added_ingredients": "NA"  # 中药成分从herb_browse提取
            })
    temp_herb_template_df = pd.DataFrame(temp_herb_template_data).drop_duplicates(subset=["Name_pinyin"])
    print(f"临时中药模板总数（去重后）：{len(temp_herb_template_df)} 个")

    # 保存临时中药模板文件（供herb2target调用）
    with pd.option_context('mode.chained_assignment', None):
        temp_herb_template_path = Path("./temp_herb_template.txt")
        temp_herb_template_df.to_csv(temp_herb_template_path, sep=sep, index=False, na_rep="NA")

    # 调用herb2target，批量获取所有中药的靶点映射
    herb2target_result = run_herb2target(
        template_path=temp_herb_template_path,
        herb_browse_path=herb_browse_path,
        herb_info_path=herb_info_path,
        known_target_path=known_target_path,
        output_pkl="./temp_herb2target_result.pkl",
        drug_likeness_min=drug_likeness_min,
        ob_score_min=ob_score_min,
        use_filtered_only=use_filtered_only,
        sep=sep
    )

    # 删除临时文件
    temp_herb_template_path.unlink()
    Path("./temp_herb2target_result.pkl").unlink(missing_ok=True)

    # ====================== 步骤5：构建方剂-靶点最终映射（汇总结果） ======================
    print("\n【步骤5：构建方剂-靶点最终映射】")
    formula2target_result = {}
    for formula_name, herbs in formula_herb_map.items():
        formula_herb_targets = {}
        for herb in herbs:
            # 获取该中药的靶点映射（无则为空字典）
            herb_targets = herb2target_result.get(herb, {})
            if herb_targets:  # 只保留有靶点的中药
                formula_herb_targets[herb] = herb_targets
        
        formula2target_result[formula_name] = formula_herb_targets

    # ====================== 步骤6：统计并保存结果 ======================
    print("\n" + "=" * 80)
    print("📊 方剂-靶点映射统计结果")
    print("=" * 80)
    print(f"总方剂数：{len(formula2target_result)} 个")
    
    for formula_name, herb_targets in formula2target_result.items():
        total_herb = len(herb_targets)
        total_ingredient = sum(len(cid_targets) for cid_targets in herb_targets.values())
        total_target = sum(len(targets) for cid_targets in herb_targets.values() for targets in cid_targets.values())
        all_targets = list(set([t for cid_targets in herb_targets.values() for targets in cid_targets.values() for t in targets]))

        print(f"\n方剂[{formula_name}]：")
        print(f"  - 有靶点的中药数：{total_herb} 个")
        print(f"  - 有靶点的成分数：{total_ingredient} 个")
        print(f"  - 总靶点数：{total_target} 个")
        print(f"  - 去重后靶点总数：{len(all_targets)} 个")
        if all_targets:
            print(f"  - 去重后靶点列表：{all_targets[:15]}...")  # 只打印前15个
        else:
            print(f"  - 去重后靶点列表：[]")

        # 打印每个中药的详细统计
        for herb, cid_targets in herb_targets.items():
            herb_ingredient_num = len(cid_targets)
            herb_target_num = sum(len(targets) for targets in cid_targets.values())
            print(f"    → 中药[{herb}]：{herb_ingredient_num}个成分 → {herb_target_num}个靶点")

    # 保存为PKL
    with open(output_pkl, "wb") as f:
        pickle.dump(formula2target_result, f)
    
    print("\n" + "=" * 80)
    print(f"✅ 方剂-靶点映射完成！")
    print(f"📁 结果已保存至：{output_pkl}")
    print(f"📋 结果格式：Dict[方剂拼音名: Dict[中药拼音名: Dict[PubChem_CID: List[靶点]]]]")
    print("=" * 80)

    return formula2target_result

# ====================== 直接运行 ======================
if __name__ == "__main__":
    run_formula2target(
        template_path="tcm_formula2target-template.txt",
        formula_browse_path="formula_browse_100.txt",
        herb_browse_path="herb_browse.txt",
        herb_info_path="HERB_ingredient_info_v2.txt",
        known_target_path="known_browse_by_ingredients.txt.gz",
        output_pkl="formula2target.pkl",
        drug_likeness_min=0.1,
        ob_score_min=10.0,
        use_filtered_only=False
    )
