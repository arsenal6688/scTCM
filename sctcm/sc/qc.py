# sctcm/sc/qc.py
import anndata as ad
import numpy as np

def calculate_qc_metrics(adata: ad.AnnData, mito_pattern: str = "^MT-") -> ad.AnnData:
    """
    计算单细胞QC核心指标
    
    参数:
        adata: AnnData对象
        mito_pattern: 线粒体基因正则匹配模式（人：^MT-，小鼠：^mt-）
    
    返回:
        追加QC指标后的AnnData对象：
        - obs: n_genes(细胞表达基因数)、n_umis(细胞UMI总数)、mito_percent(线粒体比例%)
        - var: n_cells(基因表达的细胞数)
    """
    # 细胞维度：统计每个细胞表达的基因数
    adata.obs["n_genes"] = np.sum(adata.X > 0, axis=1).A1
    # 细胞维度：统计每个细胞的UMI总数
    adata.obs["n_umis"] = np.sum(adata.X, axis=1).A1
    
    # 基因维度：统计每个基因表达的细胞数
    adata.var["n_cells"] = np.sum(adata.X > 0, axis=0).A1
    
    # 计算线粒体基因占比（基于UMI数，结果为百分比）
    mito_genes = adata.var_names.str.match(mito_pattern)
    adata.obs["mito_percent"] = (np.sum(adata.X[:, mito_genes], axis=1).A1 / adata.obs["n_umis"]) * 100
    
    return adata

def filter_cells(
    adata: ad.AnnData,
    min_genes: int = None,
    max_genes: int = None,
    max_mito: float = None,
    min_umis: int = None,
    max_umis: int = None
) -> tuple[ad.AnnData, dict]:
    """
    过滤低质量细胞（参数为None时跳过对应过滤条件）
    
    参数:
        adata: 已计算QC指标的AnnData对象（需先运行calculate_qc_metrics）
        min_genes: 细胞最少表达基因数（None=不过滤）
        max_genes: 细胞最多表达基因数（None=不过滤）
        max_mito: 线粒体比例上限（百分比，None=不过滤）
        min_umis: 细胞最少UMI数（None=不过滤）
        max_umis: 细胞最多UMI数（None=不过滤）
    
    返回:
        tuple:
            - ad.AnnData: 过滤后的AnnData对象
            - dict: 过滤统计信息（包含过滤前后数量、移除比例等）
    """
    # 记录过滤前细胞总数
    cell_count_before = adata.n_obs
    
    # 初始化过滤掩码：默认保留所有细胞
    keep_mask = np.ones(adata.n_obs, dtype=bool)
    
    # 按条件更新过滤掩码（仅参数非None时执行）
    if min_genes is not None:
        keep_mask &= adata.obs["n_genes"] >= min_genes
    if max_genes is not None:
        keep_mask &= adata.obs["n_genes"] <= max_genes
    if max_mito is not None:
        keep_mask &= adata.obs["mito_percent"] <= max_mito
    if min_umis is not None:
        keep_mask &= adata.obs["n_umis"] >= min_umis
    if max_umis is not None:
        keep_mask &= adata.obs["n_umis"] <= max_umis
    
    # 执行过滤并统计结果
    adata_filtered = adata[keep_mask, :].copy()
    cell_count_after = adata_filtered.n_obs
    cell_count_removed = cell_count_before - cell_count_after
    cell_percent_removed = (cell_count_removed / cell_count_before) * 100 if cell_count_before > 0 else 0.0
    
    # 构建过滤统计字典
    filter_stats = {
        "step": "cell_filtering",
        "count_before": cell_count_before,
        "count_after": cell_count_after,
        "count_removed": cell_count_removed,
        "percent_removed": round(cell_percent_removed, 2)
    }
    
    # 打印易读的过滤结果
    print("=== 细胞过滤统计 ===")
    print(f"过滤前细胞数: {cell_count_before}")
    print(f"过滤后细胞数: {cell_count_after}")
    print(f"移除细胞数: {cell_count_removed} ({filter_stats['percent_removed']}%)")
    print(f"过滤参数: min_genes={min_genes}, max_genes={max_genes}, "
          f"max_mito={max_mito}, min_umis={min_umis}, max_umis={max_umis}\n")
    
    return adata_filtered, filter_stats

def filter_genes(adata: ad.AnnData, min_cells: int = None) -> tuple[ad.AnnData, dict]:
    """
    过滤低表达基因（参数为None时跳过过滤）
    
    参数:
        adata: AnnData对象
        min_cells: 基因最少表达的细胞数（None=不过滤）
    
    返回:
        tuple:
            - ad.AnnData: 过滤后的AnnData对象
            - dict: 过滤统计信息（包含过滤前后数量、移除比例等）
    """
    # 记录过滤前基因总数
    gene_count_before = adata.n_vars
    
    # 初始化过滤掩码：默认保留所有基因
    keep_mask = np.ones(adata.n_vars, dtype=bool)
    if min_cells is not None:
        keep_mask = adata.var["n_cells"] >= min_cells
    
    # 执行过滤并统计结果
    adata_filtered = adata[:, keep_mask].copy()
    gene_count_after = adata_filtered.n_vars
    gene_count_removed = gene_count_before - gene_count_after
    gene_percent_removed = (gene_count_removed / gene_count_before) * 100 if gene_count_before > 0 else 0.0
    
    # 构建过滤统计字典
    filter_stats = {
        "step": "gene_filtering",
        "count_before": gene_count_before,
        "count_after": gene_count_after,
        "count_removed": gene_count_removed,
        "percent_removed": round(gene_percent_removed, 2)
    }
    
    # 打印易读的过滤结果
    print("=== 基因过滤统计 ===")
    print(f"过滤前基因数: {gene_count_before}")
    print(f"过滤后基因数: {gene_count_after}")
    print(f"移除基因数: {gene_count_removed} ({filter_stats['percent_removed']}%)")
    print(f"过滤参数: min_cells={min_cells}\n")
    return adata_filtered, filter_stats
