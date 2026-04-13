# sctcm/pl/dotplot.py
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_celltype_target_dotplot(adata: ad.AnnData, target_genes: list, celltype_key: str = "celltype", figsize: tuple = (10, 6), cmap: str = "viridis") -> plt.Figure:
    """
    绘制细胞类型-靶点点图（点大小=表达细胞比例，颜色=平均表达）
    :param adata: AnnData对象
    :param target_genes: 目标靶点列表
    :param celltype_key: 细胞类型列名
    :param figsize: 图大小
    :param cmap: 颜色映射
    :return: matplotlib Figure对象
    """
    
    target_genes_in_adata = [g for g in target_genes if g in adata.var_names]
    if len(target_genes_in_adata) == 0:
        raise ValueError("无靶点在单细胞数据中")
    
    # 计算每个细胞类型的靶点平均表达和表达比例
    expr_mean = []
    expr_frac = []
    celltypes = []
    genes = []
    for ct in adata.obs[celltype_key].unique():
        ct_adata = adata[adata.obs[celltype_key] == ct, target_genes_in_adata]
        for i, gene in enumerate(target_genes_in_adata):
            # 平均表达
            mean = ct_adata[:, gene].X.mean()
            # 表达比例（有表达的细胞数/总细胞数）
            frac = np.sum(ct_adata[:, gene].X > 0) / ct_adata.n_obs
            expr_mean.append(mean)
            expr_frac.append(frac)
            celltypes.append(ct)
            genes.append(gene)
    # 构建DataFrame
    plot_df = pd.DataFrame({
        "celltype": celltypes,
        "target": genes,
        "expr_mean": expr_mean,
        "expr_frac": expr_frac
    })
    # 绘图
    fig, ax = plt.subplots(figsize=figsize)
    sns.scatterplot(
        data=plot_df,
        x="target",
        y="celltype",
        size="expr_frac",
        hue="expr_mean",
        sizes=(50, 500),
        cmap=cmap,
        ax=ax
    )
    ax.set_xlabel("Target Genes")
    ax.set_ylabel("Cell Type")
    ax.set_title("Cell Type - Target Expression Dotplot")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    return fig
