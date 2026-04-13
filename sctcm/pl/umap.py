# sctcm/pl/umap.py
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt

def plot_umap_celltype(adata: ad.AnnData, celltype_key: str = "celltype", figsize: tuple = (8, 6), palette: str = "tab20") -> plt.Figure:
    """
    绘制细胞类型UMAP图
    :param adata: AnnData对象
    :param celltype_key: 细胞类型列名
    :param figsize: 图大小
    :param palette: 配色方案
    :return: matplotlib Figure对象
    """
    if "X_umap" not in adata.obsm:
        raise ValueError("UMAP结果不存在，请先运行sctcm.sc.run_umap")
    fig, ax = plt.subplots(figsize=figsize)
    sc.pl.umap(adata, color=celltype_key, ax=ax, palette=palette, show=False)
    ax.set_title("UMAP of Cell Types")
    plt.tight_layout()
    return fig

def plot_umap_gene_expr(adata: ad.AnnData, gene: str, figsize: tuple = (8, 6), cmap: str = "Reds") -> plt.Figure:
    """
    绘制基因表达UMAP图
    """
    if gene not in adata.var_names:
        raise ValueError(f"基因{gene}不在数据中")
    if "X_umap" not in adata.obsm:
        raise ValueError("UMAP结果不存在")
    fig, ax = plt.subplots(figsize=figsize)
    sc.pl.umap(adata, color=gene, ax=ax, cmap=cmap, show=False)
    ax.set_title(f"UMAP of {gene} Expression")
    plt.tight_layout()
    return fig
