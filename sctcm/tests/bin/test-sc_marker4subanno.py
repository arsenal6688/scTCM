#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from sctcm.sc.marker4subanno import subtype_anno_pipeline

# ===================== 路径配置 =====================
# 提取后的主细胞类型h5ad（巨噬细胞）
ADATA_PATH = "/scTCM/sctcm/tests/result/sc_subtypeProc-macrophage.h5ad"
# 自定义亚型marker字典（如巨噬细胞的Mac-LYVE1/Mac-IL1β/Mac-undefined）
SUBTYPE_MARKER_PATH = "/scTCM/docs/templates/sc_marker4anno-subtype_marker_dict.txt"
# CellMarker数据库文件
CELLMARKER_DB_PATH = "/scTCM/sctcm/data/CellMarker_Augmented_2021.txt"
# 输出目录（结果会在该目录下生成macrophage子文件夹）
ROOT_OUTPUT_DIR = "/scTCM/sctcm/tests/result"

# ===================== 【参数配置】 =====================
MAJOR_CELL_TYPE = "Macrophage"  # 主细胞类型名
GROUPBY = "subtype_leiden_res1_0"       # 亚型聚类列名
ORGANISM = "Human"              # 物种
NUM_WORKERS = 1                 # AUCell线程数
DEG_METHOD = "wilcoxon"         # 差异分析方法
TOP_N_DEGS = 10                 # 每个亚型保存的top差异基因数
DPI = 300                       # 图片分辨率

# ===================== 【一键执行亚型注释流程】 =====================
if __name__ == "__main__":
    print(f"🚀 开始{MAJOR_CELL_TYPE}亚型注释分析...")
    adata_marker, adata_aucs = subtype_anno_pipeline(
        adata_path=ADATA_PATH,
        subtype_marker_path=SUBTYPE_MARKER_PATH,
        cellmarker_db_path=CELLMARKER_DB_PATH,
        root_output_dir=ROOT_OUTPUT_DIR,
        major_cell_type=MAJOR_CELL_TYPE,
        groupby=GROUPBY,
        organism=ORGANISM,
        num_workers=NUM_WORKERS,
        deg_method=DEG_METHOD,
        top_n_degs=TOP_N_DEGS,
        dpi=DPI
    )
    print(f"✅ {MAJOR_CELL_TYPE}亚型注释测试流程执行完成！")
