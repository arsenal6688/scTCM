import os
import omicverse as ov
import scanpy as sc
from sctcm.sc.AUCell4anno import aucell_core_pipeline

# ===================== 路径配置 =====================
ADATA_PATH = "/scTCM/sctcm/tests/result/sc_cluster.h5ad"
OUTPUT_DIR = "/scTCM/sctcm/tests/result"
CELLMARKER_DB = "/scTCM/sctcm/data/CellMarker_Augmented_2021.txt"

adata = ov.read(ADATA_PATH)

# ===================== 运行 AUCell 注释 =====================
aucell_core_pipeline(
    adata_path= ADATA_PATH,
    cellmarker_db_path= CELLMARKER_DB,
    save_dir= OUTPUT_DIR,
    groupby="leiden_res1_0"
)
