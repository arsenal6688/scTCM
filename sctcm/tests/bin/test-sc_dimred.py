from sctcm import dimred

# 核心配置
adata_path = "/scTCM/sctcm/tests/result/sc_io-merged_samples_qc.h5ad"  # 输入：已QC的h5ad文件路径
save_dir = "/scTCM/sctcm/tests/result"  
custom_config = {
    "preprocess_params": {"n_hvgs": 2000},
    "pca_params": {"n_pcs": 30},  
    "batch_correction_params": {
        "method": "harmony",
        "n_pcs": 30  
    }
}


dimred.preprocess_dimred(adata_path, save_dir, custom_config)

