#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, "/scTCM")

from sctcm.tcm.form2targetall import run_form2targetall

# ===================== 【路径配置】 =====================
# formula_browse_100.txt 路径
FORMULA_BROWSE_PATH = "/scTCM/sctcm/data/BATMAN-TCM/4-Formulas_information_and_their_herbal_components/formula_browse.txt"
# herb2targetall.py 生成的中药-靶点字典路径
HERB2TARGET_PKL_PATH = "/scTCM/sctcm/tests/result/tcm_herb2targetall.pkl"
# 输出目录
OUTPUT_DIR = "/scTCM/sctcm/tests/result"

# ===================== 【执行测试】 =====================
if __name__ == "__main__":
    print("="*70)
    print("          开始构建全库方剂→中药→成分→靶点映射")
    print("="*70)
    print(f"📋 配置信息：")
    print(f"  - 方剂数据文件：{FORMULA_BROWSE_PATH}")
    print(f"  - 中药-靶点字典：{HERB2TARGET_PKL_PATH}")
    print(f"  - 输出目录：{OUTPUT_DIR}")
    print(f"  - 输出文件：tcm_form2targetall.pkl")
    print("="*70)
    
    # 核心调用
    form2target = run_form2targetall(
        formula_browse_path=FORMULA_BROWSE_PATH,
        herb2target_pkl_path=HERB2TARGET_PKL_PATH,
        output_dir=OUTPUT_DIR
    )
    
    # 输出结果统计与示例
    print("\n" + "="*70)
    print("                     运行完成")
    print(f"📊 最终结果：")
    print(f"  - 成功关联靶点的方剂数：{len(form2target)}")
    
    # 输出前3个方剂示例
    print(f"\n📌 前3个方剂示例：")
    count = 0
    for formula, herbs in form2target.items():
        if count >= 3:
            break
        print(f"\n  方剂：{formula}")
        print(f"    包含有靶点的中药数：{len(herbs)}")
        for herb, cids in herbs.items():
            print(f"      → 中药[{herb}]：{len(cids)}个成分，{sum(len(t) for t in cids.values())}个靶点")
        count += 1
    
    print("\n" + "="*70)
    print(f"✅ 输出格式：Dict[方剂拼音: Dict[中药拼音: Dict[CID: List[靶点]]]]")
    print("="*70)
