#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, "/scTCM")

from sctcm.tcm.herb2targetall import run_herb2targetall

# ===================== 路径配置 =====================
HERB_BROWSE    = "/scTCM/sctcm/data/BATMAN-TCM/3-Herbs_information_and_their_ingredient_compounds/herb_browse.txt"
ING_TARGET_PKL = "/scTCM/sctcm/tests/result/tcm_ing2targetall-ing2targetall.pkl"
OUTPUT_DIR     = "/scTCM//sctcm/tests/result"

# ===================== 运行 =====================
if __name__ == "__main__":
    print("="*60)
    print("        全库中药 → 靶点映射")
    print("="*60)

    result = run_herb2targetall(
        herb_browse_path=HERB_BROWSE,
        ing_target_pkl=ING_TARGET_PKL,
        output_dir=OUTPUT_DIR
    )

    print(f"\n最终成功：{len(result)} 个中药")
    print("输出文件：tcm_herb2targetall.pkl")
    print("="*60)
