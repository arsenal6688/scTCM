#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
form2targetall.py
功能：从 formula_browse_100.txt 提取所有方剂 → 中药 → 成分 → 靶点映射
步骤：
1) 读取 formula_browse_100.txt 所有方剂及对应中药（拼音名）
2) 复用 herb2targetall.py 逻辑，批量获取中药-成分-靶点映射
3) 按方剂聚合所有中药的靶点数据
4) 输出嵌套式方剂-靶点字典（PKL格式）
"""

import os
import pandas as pd
import pickle
import warnings
warnings.filterwarnings("ignore")

__all__ = ["run_form2targetall"]

def parse_formula_herbs(formula_browse_path: str) -> dict:
    """解析 formula_browse_100.txt，提取方剂-中药映射"""
    print("读取并解析 formula_browse_100.txt ...")
    formula_df = pd.read_csv(formula_browse_path, sep="\t", encoding="utf-8")
    
    # 列校验
    required_cols = ["Pinyin.Name", "Pinyin.composition"]
    for col in required_cols:
        if col not in formula_df.columns:
            raise ValueError(f"formula_browse_100.txt 缺少必需列：{col}")
    
    # 过滤无效行
    formula_df = formula_df.dropna(subset=required_cols).copy()
    formula_df = formula_df[formula_df["Pinyin.composition"] != ""].copy()
    formula_df = formula_df.drop_duplicates(subset=["Pinyin.Name"], keep="first")
    
    # 解析中药：按逗号分割多个中药拼音
    formula_herb_map = {}
    for _, row in formula_df.iterrows():
        formula_pinyin = str(row["Pinyin.Name"]).strip()
        herbs_str = str(row["Pinyin.composition"]).strip()
        # 分割中药并去重、过滤空值
        herbs = list(set([h.strip() for h in herbs_str.split(",") if h.strip()]))
        if herbs:
            formula_herb_map[formula_pinyin] = herbs
    
    print(f"有效方剂数量：{len(formula_herb_map)}")
    print(f"涉及中药总数（去重后）：{len(set([h for herbs in formula_herb_map.values() for h in herbs]))}")
    return formula_herb_map

def load_herb2target_dict(herb2target_pkl_path: str) -> dict:
    """加载中药-成分-靶点字典（herb2targetall.py生成的PKL）"""
    print(f"\n加载中药-靶点字典：{herb2target_pkl_path}")
    if not os.path.exists(herb2target_pkl_path):
        raise FileNotFoundError(f"中药-靶点字典文件不存在，请先运行herb2targetall.py生成")
    
    with open(herb2target_pkl_path, "rb") as f:
        herb2target = pickle.load(f)
    
    print(f"中药-靶点字典包含中药数：{len(herb2target)}")
    return herb2target

def build_form2target(formula_herb_map: dict, herb2target: dict) -> dict:
    """构建方剂-靶点映射（嵌套结构：方剂→中药→成分→靶点）"""
    print("\n开始构建方剂-靶点映射...")
    form2target = {}
    
    for formula_pinyin, herbs in formula_herb_map.items():
        formula_herb_targets = {}
        for herb_pinyin in herbs:
            # 获取该中药的成分-靶点映射（无则跳过）
            herb_target_data = herb2target.get(herb_pinyin, {})
            if herb_target_data:  # 只保留有靶点的中药
                formula_herb_targets[herb_pinyin] = herb_target_data
        
        # 只保留有有效中药靶点的方剂
        if formula_herb_targets:
            form2target[formula_pinyin] = formula_herb_targets
    
    print(f"成功关联靶点的方剂数量：{len(form2target)}")
    return form2target

def run_form2targetall(
    formula_browse_path: str,
    herb2target_pkl_path: str,
    output_dir: str,
    output_filename: str = "tcm_form2targetall.pkl"
) -> dict:
    """
    一键执行全库方剂-靶点映射流程
    :param formula_browse_path: formula_browse_100.txt路径
    :param herb2target_pkl_path: 中药-靶点字典PKL路径（herb2targetall.py输出）
    :param output_dir: 输出目录
    :param output_filename: 输出PKL文件名
    :return: 方剂-靶点映射字典
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 步骤1：解析方剂-中药映射
    formula_herb_map = parse_formula_herbs(formula_browse_path)
    
    # 步骤2：加载中药-靶点字典
    herb2target = load_herb2target_dict(herb2target_pkl_path)
    
    # 步骤3：构建方剂-靶点映射
    form2target = build_form2target(formula_herb_map, herb2target)
    
    # 步骤4：保存结果（仅PKL格式）
    output_path = os.path.join(output_dir, output_filename)
    with open(output_path, "wb") as f:
        pickle.dump(form2target, f)
    
    # 输出统计摘要
    print(f"\n✅ 方剂-靶点字典已保存：{output_path}")
    print("\n📊 映射统计摘要：")
    print(f"  - 输入方剂总数：{len(formula_herb_map)}")
    print(f"  - 成功关联靶点的方剂数：{len(form2target)}")
    total_herbs = sum(len(herbs) for herbs in form2target.values())
    print(f"  - 涉及有靶点的中药总数：{total_herbs}")
    total_components = sum(len(cids) for herbs in form2target.values() for cids in herbs.values())
    print(f"  - 涉及有靶点的成分总数：{total_components}")
    all_targets = set()
    for herbs in form2target.values():
        for cids in herbs.values():
            all_targets.update([t for targets in cids.values() for t in targets])
    print(f"  - 去重后总靶点数：{len(all_targets)}")
    
    print("\n🎉 全库方剂-靶点映射构建完成！")
    return form2target

if __name__ == "__main__":
    run_form2targetall(
        formula_browse_path="./formula_browse_100.txt",
        herb2target_pkl_path="./tcm_herb2targetall.pkl",
        output_dir="./form_target_results",
        output_filename="tcm_form2targetall.pkl"
    )
