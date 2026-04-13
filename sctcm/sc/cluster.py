import os
import omicverse as ov
import scanpy as sc
import matplotlib.pyplot as plt
from typing import Optional, List, Union


def build_neighbors(
    adata,
    n_neighbors: int = 15,  # 根据细胞数调整（10-30）
    n_pcs: int = 30,        # 需≤前期PCA保留的主成分数
    use_rep: str = "X_harmony",  # 替换为X_pca/X_scaled等降维矩阵
    key_added: Optional[str] = None,
    **kwargs
) -> None:
    """
    基于降维结果构建细胞邻域图（为聚类做准备）
    
    Parameters
    ----------
    adata : anndata.AnnData
        经过降维的单细胞数据对象
    n_neighbors : int, default=15
        构建邻域图的近邻数（可自行设置，常用10-30）
    n_pcs : int, default=30
        使用的主成分数量（可自行设置，不超过PCA维度）
    use_rep : str, default="X_harmony"
        用于构建邻域图的降维结果矩阵（可自行替换：X_pca/X_harmony/X_scaled）
    key_added : str, optional
        邻域图结果的存储键（默认由scanpy自动命名）
    **kwargs
        传递给sc.pp.neighbors的其他参数（如metric等）
    
    Returns
    -------
    None
        直接修改adata，添加邻域图相关数据
    """

    print(f"开始构建邻域图：n_neighbors={n_neighbors}, n_pcs={n_pcs}, use_rep={use_rep}")
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=use_rep,
        key_added=key_added,
        **kwargs
    )
    print("邻域图构建完成")


def leiden_cluster(
    adata,
    resolutions: Union[List[float], float] = 1.0,  # 默认仅保留1.0
    key_prefix: str = "leiden_res",
    neighbors_key: Optional[str] = None,
    **kwargs
) -> None:
    """
    执行Leiden聚类（支持多分辨率，默认仅1.0）
    
    Parameters
    ----------
    adata : anndata.AnnData
        已构建邻域图的单细胞数据对象
    resolutions : list or float, default=1.0
        聚类分辨率（可自行设置：如0.5/2.0，值越大聚类越细）
    key_prefix : str, default="leiden_res"
        聚类结果在adata.obs中的存储键前缀（如resolution=1.0则存储为leiden_res1_0）
    neighbors_key : str, optional
        邻域图的存储键（默认由scanpy自动命名）
    **kwargs
        传递给sc.tl.leiden的其他参数（如random_state等）
    
    Returns
    -------
    None
        直接修改adata，在obs中添加各分辨率的聚类结果
    """

    # 标准化resolutions为列表
    if isinstance(resolutions, float) or isinstance(resolutions, int):
        resolutions = [float(resolutions)]
    
    # 移除邻域图存在性判断逻辑
    neighbors_key = neighbors_key or "neighbors"
    
    # 执行聚类（默认1.0分辨率）
    for res in resolutions:
        res_key = f"{key_prefix}{str(res).replace('.', '_')}"
        print(f"\n执行Leiden聚类：resolution={res}, 结果存储键={res_key}")
        
        sc.tl.leiden(
            adata,
            resolution=res,
            key_added=res_key,
            neighbors_key=neighbors_key,
            **kwargs
        )
        
        # ===== 统计并打印聚类结果 =====
        cluster_counts = adata.obs[res_key].value_counts().sort_index()
        total_clusters = len(cluster_counts)
        total_cells = cluster_counts.sum()
        
        print(f"├─ 分辨率 {res} 聚类结果统计 ─────────────────")
        print(f"│  总cluster数量：{total_clusters}")
        print(f"│  总细胞数量：{total_cells}")
        print(f"│  各cluster细胞分布：")
        for cluster_id, count in cluster_counts.items():
            percentage = (count / total_cells) * 100
            print(f"│    Cluster {cluster_id}: {count} 个细胞 ({percentage:.2f}%)")
        print(f"└────────────────────────────────────────────")
    
    print(f"\n所有分辨率聚类完成：{resolutions}")


def run_mde(
    adata,
    use_rep: str = "X_harmony",  # 可自定义：X_pca/X_scaled等
    mde_key: str = "X_mde",
    **kwargs
) -> None:
    """
    执行MDE降维（用于聚类结果可视化）
    
    Parameters
    ----------
    adata : anndata.AnnData
        单细胞数据对象
    use_rep : str, default="X_harmony"
        用于MDE降维的输入矩阵（可自行替换：X_pca/X_harmony等）
    mde_key : str, default="X_mde"
        MDE结果在adata.obsm中的存储键
    **kwargs
        传递给ov.utils.mde的其他参数
    
    Returns
    -------
    None
        直接修改adata，在obsm中添加MDE降维结果
    """

    # 移除输入矩阵存在性判断
    print(f"\n执行MDE降维：输入矩阵={use_rep}, 结果存储键={mde_key}")
    adata.obsm[mde_key] = ov.utils.mde(adata.obsm[use_rep], **kwargs)
    print("MDE降维完成")


def plot_cluster(
    adata,
    basis: str = "X_mde",
    cluster_keys: Union[List[str], str] = "leiden_res1_0",  # 匹配1.0分辨率的存储键
    titles: Optional[List[str]] = None,
    save_path: Optional[str] = None,
    dpi: int = 300,
    palette: Optional[List[str]] = None,
    frameon: str = "small",
    **kwargs
) -> None:
    """
    可视化聚类结果（基于MDE/UMAP等降维空间）
    
    Parameters
    ----------
    adata : anndata.AnnData
        包含聚类结果和降维数据的单细胞数据对象
    basis : str, default="X_mde"
        可视化的降维空间（可自行设置：X_mde/X_umap/X_tsne）
    cluster_keys : list or str, default="leiden_res1_0"
        要可视化的聚类结果键（对应leiden_cluster的resolutions=1.0）
    titles : list, optional
        每个子图的标题（需与cluster_keys长度一致）
    save_path : str, optional
        图片保存路径（如指定则自动保存并关闭画布）
    dpi : int, default=300
        图片分辨率（可自行调整）
    palette : list, optional
        配色方案（默认使用omicverse默认配色）
    frameon : str, default="small"
        边框样式（参考omicverse.embedding参数）
    **kwargs
        传递给ov.utils.embedding的其他参数
    
    Returns
    -------
    None
        显示或保存可视化结果
    """

    # 标准化参数
    if isinstance(cluster_keys, str):
        cluster_keys = [cluster_keys]
    if titles is None:
        titles = [f"Resolution: {k.replace('leiden_res', '').replace('_', '.')}" for k in cluster_keys]
    if len(titles) != len(cluster_keys):
        raise ValueError("titles长度必须与cluster_keys一致！")
    if palette is None:
        palette = ov.palette()[:]
    
    
    print(f"\n可视化聚类结果：降维空间={basis}, 聚类键={cluster_keys}")
    ov.utils.embedding(
        adata,
        basis=basis,
        color=cluster_keys,
        title=titles,
        palette=palette,
        show=False,
        frameon=frameon,
        **kwargs
    )
    
    # 保存或显示图片
    if save_path is not None:
        save_dir = os.path.dirname(save_path)
        if save_dir and not os.path.exists(save_dir):
            os.makedirs(save_dir, exist_ok=True)
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        print(f"聚类图已保存至：{save_path}")
        plt.close()
    else:
        plt.show()


def run_cluster_pipeline(
    adata,
    neighbors_params: dict = None,
    leiden_params: dict = None,
    mde_params: dict = None,
    plot_params: dict = None,
    save_adata_path: Optional[str] = None
) -> None:
    """
    一站式执行聚类全流程：构建邻域图 → Leiden聚类 → MDE降维 → 可视化
    
    Parameters
    ----------
    adata : anndata.AnnData
        经过降维的单细胞数据对象
    neighbors_params : dict, optional
        传递给build_neighbors的参数（可自定义n_neighbors/n_pcs/use_rep等）
    leiden_params : dict, optional
        传递给leiden_cluster的参数（默认resolutions=1.0）
    mde_params : dict, optional
        传递给run_mde的参数（可自定义use_rep/mde_key等）
    plot_params : dict, optional
        传递给plot_cluster的参数（可自定义basis/cluster_keys等）
    save_adata_path : str, optional
        聚类完成后adata的保存路径（h5ad格式）
    
    Returns
    -------
    None
        完成全流程并可选保存adata
    """
    # 默认参数（resolutions固定为1.0，其他参数均可自定义）
    neighbors_params = neighbors_params or {"n_neighbors": 15, "n_pcs": 30, "use_rep": "X_harmony"}
    leiden_params = leiden_params or {"resolutions": 1.0, "key_prefix": "leiden_res"}
    mde_params = mde_params or {"use_rep": "X_harmony", "mde_key": "X_mde"}
    plot_params = plot_params or {"basis": "X_mde", "cluster_keys": "leiden_res1_0"}
    
    # 执行流程
    print("="*60)
    print("开始执行聚类全流程")
    print("="*60)
    build_neighbors(adata, **neighbors_params)
    leiden_cluster(adata, **leiden_params)
    run_mde(adata, **mde_params)
    plot_cluster(adata, **plot_params)
    
    # 保存adata
    if save_adata_path is not None:
        save_dir = os.path.dirname(save_adata_path)
        if save_dir and not os.path.exists(save_dir):
            os.makedirs(save_dir, exist_ok=True)
        adata.write(save_adata_path)
        print(f"\n聚类后的数据已保存至：{save_adata_path}")
    
    print("\n" + "="*60)
    print("聚类全流程执行完成")
    print("="*60)
