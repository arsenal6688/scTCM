#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import omicverse as ov
import scanpy as sc
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

# 保留核心绘图样式设置
ov.ov_plot_set()


def load_cellmarker_geneset(db_path: str, organism: str = "Human") -> dict:
    """
    加载CellMarker数据库基因集
    :param db_path: CellMarker_Augmented_2021.txt路径
    :param organism: 物种（默认Human）
    :return: 细胞类型-基因集字典
    """
    print(f"\n[1/5] 加载CellMarker数据库: {db_path}")
    pathway_dict = ov.utils.geneset_prepare(
        db_path,
        organism=organism
    )
    # 打印基因集基础信息
    print(f"✅ 加载 {len(pathway_dict)} 个细胞类型基因集")
    return pathway_dict


def run_aucell_enrichment(
    adata_path: str,
    pathway_dict: dict,
    num_workers: int = 1
) -> sc.AnnData:
    """
    执行AUCell基因集富集
    :param adata_path: 原始h5ad文件路径
    :param pathway_dict: CellMarker基因集字典
    :param num_workers: 线程数
    :return: 含AUCell分数的adata_aucs
    """
    print(f"\n[2/5] 加载单细胞数据并执行AUCell富集")
    # 加载数据并提取raw矩阵
    adata = ov.read(adata_path)
    adata_raw = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    
    # 核心AUCell计算
    adata_aucs = ov.single.pathway_aucell_enrichment(
        adata_raw,
        pathways_dict=pathway_dict,
        num_workers=num_workers,
    )
    
    # 对齐元信息（保证分组信息一致）
    adata_aucs.obs = adata[adata_aucs.obs.index].obs
    adata_aucs.obsm = adata[adata_aucs.obs.index].obsm
    adata_aucs.obsp = adata[adata_aucs.obs.index].obsp
    
    print(f"✅ AUCell完成，结果维度: {adata_aucs.shape}")
    return adata_aucs


def run_aucell_diff_analysis(
    adata_aucs: sc.AnnData,
    groupby: str = "leiden_res1",
    method: str = "wilcoxon",
    key_added: str = "dea_leiden_aucs_res1"
) -> sc.AnnData:
    """
    AUCell结果差异分析（核心步骤）
    :param adata_aucs: AUCell富集结果
    :param groupby: 分组列名（默认leiden_res1）
    :param method: 差异分析方法（默认wilcoxon）
    :param key_added: 结果存储key
    :return: 含差异分析结果的adata_aucs
    """
    print(f"\n[3/5] 执行AUCell差异分析 (groupby={groupby})")
    # 修复log1p base问题
    if 'log1p' not in adata_aucs.uns:
        adata_aucs.uns['log1p'] = {'base': None}
    else:
        adata_aucs.uns['log1p']['base'] = None
    
    # 核心差异分析
    sc.tl.rank_genes_groups(
        adata_aucs,
        groupby=groupby,
        use_raw=False,
        method=method,
        key_added=key_added
    )
    print("✅ 差异分析完成")
    return adata_aucs


def plot_aucell_dotplot(
    adata_aucs: sc.AnnData,
    save_path: str,
    groupby: str = "leiden_res1",
    key: str = "dea_leiden_aucs_res1",
    n_genes: int = 3,
    cmap: str = "RdBu_r"
):
    """
    绘制AUCell dotplot（核心可视化）
    :param adata_aucs: 含差异分析结果的adata
    :param save_path: dotplot保存路径
    :param groupby: 分组列名
    :param key: 差异分析结果key
    :param n_genes: 每个分组展示基因数
    :param cmap: 配色
    """
    print(f"\n[4/5] 绘制AUCell dotplot")
    # 核心绘图逻辑
    sc.pl.rank_genes_groups_dotplot(
        adata_aucs,
        groupby=groupby,
        key=key,
        dendrogram=True,
        standard_scale="var",
        n_genes=n_genes,
        cmap=cmap,
        show=False
    )
    # 保存图片
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Dotplot已保存至: {save_path}")


def save_aucell_cluster_annotations(
    adata_aucs: sc.AnnData,
    save_path: str,
    groupby: str = "leiden_res1",
    key: str = "dea_leiden_aucs_res1"
):
    """
    保存AUCell分析后每个cluster的富集细胞类型
    :param adata_aucs: 含差异分析结果的adata
    :param save_path: 结果保存路径
    :param groupby: 分组列名
    :param key: 差异分析结果key
    """
    print(f"\n[5/5] 保存AUCell cluster注释结果")
    # 获取dendrogram排序的cluster（参考原代码）
    clusters = (
        adata_aucs.uns[f'dendrogram_{groupby}']['categories_ordered']
        if f'dendrogram_{groupby}' in adata_aucs.uns
        else sorted(adata_aucs.obs[groupby].cat.categories)
    )
    
    # 保存每个cluster的top2富集细胞类型
    with open(save_path, 'w') as f:
        for cluster in clusters:
            special_cluster = str(cluster)
            degs = sc.get.rank_genes_groups_df(
                adata_aucs,
                group=special_cluster,
                key=key
            ).dropna()
            f.write(f'{special_cluster}: {"|".join(degs.names[:2].tolist())}\n')
    print(f"✅ Cluster注释已保存至: {save_path}")


def aucell_core_pipeline(
    adata_path: str,
    cellmarker_db_path: str,
    save_dir: str,
    groupby: str = "leiden_res1",
    organism: str = "Human",
    num_workers: int = 1
):
    """
    :param adata_path: 原始h5ad文件路径
    :param cellmarker_db_path: CellMarker数据库路径
    :param save_dir: 结果保存目录
    :param groupby: 分组列名（默认leiden_res1）
    :param organism: 物种（默认Human）
    :param num_workers: AUCell线程数（默认1）
    """
    # 创建保存目录
    os.makedirs(save_dir, exist_ok=True)
    
    # 1. 加载CellMarker基因集
    pathway_dict = load_cellmarker_geneset(cellmarker_db_path, organism)
    
    # 2. 执行AUCell富集
    adata_aucs = run_aucell_enrichment(adata_path, pathway_dict, num_workers)
    
    # 3. AUCell差异分析
    adata_aucs = run_aucell_diff_analysis(adata_aucs, groupby)
    
    # 4. 绘制并保存AUCell dotplot
    dotplot_path = os.path.join(save_dir, "sc_AUCell4anno-aucell_dotplot.png")
    plot_aucell_dotplot(adata_aucs, dotplot_path, groupby)
    
    # 5. 保存cluster注释结果
    cluster_anno_path = os.path.join(save_dir, "sc_AUCell4anno-aucell_cluster_markers.txt")
    save_aucell_cluster_annotations(adata_aucs, cluster_anno_path, groupby)
    
    print("\n🎉 AUCell核心流程执行完成！")
    print(f"📁 结果保存目录: {save_dir}")


# 主函数
if __name__ == "__main__":
    aucell_core_pipeline(
        adata_path="cluster.h5ad",
        cellmarker_db_path="CellMarker_Augmented_2021.txt",
        save_dir="/result",
        groupby="leiden_res1",
        organism="Human",
        num_workers=1
    )
