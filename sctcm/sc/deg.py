# sctcm/sc/deg.py
import anndata as ad
import scanpy as sc
import pandas as pd

def calculate_deg(adata: ad.AnnData, groupby: str = "celltype", reference: str = "rest", method: str = "wilcoxon") -> pd.DataFrame:
    """
    计算差异表达基因
    :param adata: AnnData对象（含细胞类型注释）
    :param groupby: 分组列名（如celltype）
    :param reference: 参考组（rest表示与其他所有组比较）
    :param method: 差异分析方法（wilcoxon/t-test）
    :return: DEG结果DataFrame（含log2FC、pvalue、padj等）
    """
    sc.tl.rank_genes_groups(adata, groupby=groupby, reference=reference, method=method)
    # 转换为DataFrame
    deg_df = sc.get.rank_genes_groups_df(adata, group=None)
    return deg_df

def get_deg_genes(deg_df: pd.DataFrame, group: str, padj_cutoff: float = 0.05, log2fc_cutoff: float = 1.0) -> list:
    """
    筛选显著差异基因
    :param deg_df: calculate_deg返回的DataFrame
    :param group: 目标细胞类型组
    :param padj_cutoff: 校正后p值阈值
    :param log2fc_cutoff: log2倍变化阈值
    :return: 显著差异基因列表
    """
    deg_group = deg_df[deg_df["group"] == group]
    sig_genes = deg_group[(deg_group["padj"] < padj_cutoff) & (abs(deg_group["log2fc"]) > log2fc_cutoff)]["names"].tolist()
    return sig_genes
