import omicverse as ov
import scanpy as sc
import os
import matplotlib.pyplot as plt

# -------------------------- 1. 配置路径 --------------------------
INPUT_H5AD = '/scTCM/sctcm/tests/result/sc_dimred-scale_HVGs_PCA_harmony.h5ad'
OUTPUT_PNG = '/scTCM/sctcm/tests/result/sc_cluster.leiden_resolution.png'
OUTPUT_H5AD = '/scTCM/sctcm/tests/result/sc_cluster.h5ad'

# -------------------------- 2. 导入封装的聚类函数 --------------------------

from sctcm.sc.cluster import (
    build_neighbors,
    leiden_cluster,
    run_mde,
    plot_cluster,
    run_cluster_pipeline
)

# -------------------------- 3. 核心测试逻辑 --------------------------
# 3.1 读取原始数据
adata = ov.read(INPUT_H5AD)

# 3.2 执行一站式聚类流程(复用封装函数)
run_cluster_pipeline(
    adata=adata,
    neighbors_params={
        "n_neighbors": 15,
        "n_pcs": 30,
        "use_rep": "X_harmony"
    },
    leiden_params={
        "resolutions": [0.5, 1.0],  
        "key_prefix": "leiden_res"
    },
    mde_params={
        "use_rep": "X_harmony",
        "mde_key": "X_mde"
    },
    plot_params={
        "basis": "X_mde",

        "cluster_keys": ["leiden_res0_5", "leiden_res1_0"],

        "titles": ["Resolution:0.5", "Resolution:1.0"],
        "save_path": OUTPUT_PNG,
        "dpi": 300,
        "frameon": "small"
    },
    save_adata_path=OUTPUT_H5AD
)

print("✅ 聚类流程执行完成")

