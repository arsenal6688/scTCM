import omicverse as ov
import scanpy as sc
import os
import matplotlib.pyplot as plt
import numpy as np
from typing import Optional, Literal, List
from anndata import AnnData

# ========== 配置项（统一管理参数） ==========
DEFAULT_CONFIG = {
    "preprocess_params": {
        "mode": "shiftlog|pearson",
        "n_hvgs": 3000,
        "batch_key": "Batch"
    },
    "pca_params": {
        "layer": "scaled",
        "n_pcs": 50
    },
    "batch_correction_params": {
        "method": "harmony",  # 支持后续扩展：combat、scanorama等
        "n_pcs": 50
    },
    "embedding_params": {
        "basis": "mde",
        "plot_color": ["Batch"]
    }
}

# ========== 函数 ==========
def preprocess_scale_hvgs(adata, preprocess_params: dict) -> AnnData:
    """预处理：仅负责计数保存、归一化、HVGs筛选、缩放"""
    # 保存原始计数
    ov.utils.store_layers(adata, layers='counts')
    # 预处理（归一化+HVGs筛选）
    adata = ov.pp.preprocess(adata, 
                             mode=preprocess_params["mode"],
                             n_HVGs=preprocess_params["n_hvgs"],
                             batch_key=preprocess_params["batch_key"])
    # 筛选HVGs+缩放
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable_features].copy()
    ov.pp.scale(adata)
    return adata

def run_pca(adata, pca_params: dict) -> AnnData:
    """PCA降维：仅负责PCA计算"""
    ov.pp.pca(adata, layer=pca_params["layer"], n_pcs=pca_params["n_pcs"])
    return adata

def run_batch_correction(adata, batch_correction_params: dict, batch_key: str) -> AnnData:
    """批次校正：支持扩展多种方法"""
    method = batch_correction_params["method"]
    n_pcs = batch_correction_params["n_pcs"]
    if method == "harmony":
        ov.single.batch_correction(adata, batch_key=batch_key, methods='harmony', n_pcs=n_pcs)
    elif method == "combat":  # 扩展示例：新增Combat方法
        ov.single.batch_correction(adata, batch_key=batch_key, methods='combat', n_pcs=n_pcs)
    else:
        raise ValueError(f"不支持的批次校正方法：{method}")
    return adata

def plot_embedding(adata, basis: str, save_dir: str, filename: str, plot_color: List[str]):
    """可视化：仅负责降维图绘制与保存（直接拼接sc_dimred-前缀）"""
    # 硬编码前缀：sc_dimred-
    full_filename = f"sc_dimred-{filename}"
    save_path = os.path.join(save_dir, full_filename)
    ov.utils.embedding(adata, basis=basis, frameon='small', color=plot_color, show=False)
    plt.savefig(save_path, bbox_inches="tight", dpi=300)
    plt.close()

# ========== 顶层调度函数 ==========
def preprocess_dimred(
    adata_path: str,
    save_dir: str,
    config: Optional[dict] = None
) -> AnnData:
    """
    顶层调度函数：整合所有细粒度步骤（基于已完成QC的数据集）
    :param adata_path: 输入h5ad路径（已完成QC）
    :param save_dir: 结果保存目录
    :param config: 自定义配置（覆盖默认值）
    :return: 处理后的AnnData
    """
    # 合并默认配置与自定义配置
    final_config = DEFAULT_CONFIG.copy()
    if config is not None:
        for k, v in config.items():
            if isinstance(v, dict) and k in final_config:
                final_config[k].update(v)
            else:
                final_config[k] = v
    
    # 1. 读取数据（已完成QC）
    adata = ov.read(adata_path)
    
    # 2. 预处理（scale+HVGs）
    adata = preprocess_scale_hvgs(adata, final_config["preprocess_params"])
    print(f"预处理后保留细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
    
    # 3. PCA降维
    adata = run_pca(adata, final_config["pca_params"])
    
    # 4. 批次校正前可视化（硬编码前缀）
    adata.obsm["X_mde_pca"] = ov.utils.mde(adata.obsm["scaled|original|X_pca"])
    plot_embedding(
        adata,
        basis="X_mde_pca",
        save_dir=save_dir,
        filename="before_batch_correction.png",
        plot_color=final_config["embedding_params"]["plot_color"]
    )
    
    # 5. 批次校正
    adata = run_batch_correction(
        adata,
        final_config["batch_correction_params"],
        final_config["preprocess_params"]["batch_key"]
    )
    
    # 6. 批次校正后可视化
    correction_method = final_config["batch_correction_params"]["method"]
    adata.obsm[f"X_mde_{correction_method}"] = ov.utils.mde(adata.obsm[f"X_{correction_method}"])
    plot_embedding(
        adata,
        basis=f"X_mde_{correction_method}",
        save_dir=save_dir,
        filename=f"after_{correction_method}_correction.png",
        plot_color=final_config["embedding_params"]["plot_color"]
    )
    
    # 7. 保存最终结果
    h5ad_filename = f"sc_dimred-scale_HVGs_PCA_{correction_method}.h5ad"
    adata.write(os.path.join(save_dir, h5ad_filename))
    return adata
