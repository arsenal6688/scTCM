#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
单细胞亚型注释专用工具（marker基因+AUCell双流程）
核心功能：
1. 基于自定义亚型marker字典做已知marker分析
2. 基于CellMarker数据库做AUCell富集分析
3. 自动生成主细胞类型命名的结果文件夹，结果独立归档
适用场景：提取单一主细胞类型后的亚型分析（如巨噬细胞、T细胞亚型）
"""
import os
import pickle
import warnings
from typing import Dict, List, Optional, Union
import scanpy as sc
import omicverse as ov
import pandas as pd
import matplotlib.pyplot as plt
from anndata import AnnData

# 全局配置
warnings.filterwarnings('ignore')
ov.ov_plot_set()
__all__ = [
    "load_marker_dict", "check_marker_genes", "plot_marker_dotplot",
    "compute_cluster_degs", "plot_deg_dotplot", "save_cluster_top_markers",
    "load_cellmarker_geneset", "run_aucell_enrichment", "run_aucell_diff_analysis",
    "plot_aucell_dotplot", "save_aucell_cluster_annotations", "subtype_anno_pipeline"
]

# ===================== 【模块1：自定义亚型marker分析（源自marker4anno.py）】 =====================
def load_marker_dict(
    marker_file_path: str,
    cell_type_col: Union[int, str] = 0,
    genes_col: Union[int, str] = 1,
    gene_sep: str = "|",
    encoding: str = "utf-8"
) -> Dict[str, List[str]]:
    marker_dict = {}
    if not os.path.exists(marker_file_path):
        raise FileNotFoundError(f"Marker文件不存在: {marker_file_path}")
    with open(marker_file_path, 'r', encoding=encoding) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                print(f"⚠️ 第{line_num}行格式错误，跳过: {line}")
                continue
            cell_type = parts[cell_type_col].strip()
            genes_str = parts[genes_col].strip()
            genes = [gene.strip() for gene in genes_str.split(gene_sep) if gene.strip()]
            if cell_type and genes:
                marker_dict[cell_type] = genes
            else:
                print(f"⚠️ 第{line_num}行数据无效，跳过: {line}")
    if not marker_dict:
        raise ValueError("未从文件中加载到有效marker数据")
    print(f"✅ 成功加载{len(marker_dict)}种亚型的marker基因")
    return marker_dict

def check_marker_genes(
    adata: AnnData,
    marker_dict: Dict[str, List[str]],
    verbose: bool = True
) -> Dict[str, List[str]]:
    filtered_marker_dict = {}
    del_markers = []
    for ct, markers in marker_dict.items():
        markers_found = [marker for marker in markers if marker in adata.var.index]
        filtered_marker_dict[ct] = markers_found
        if verbose:
            missing = set(markers) - set(markers_found)
            if missing:
                print(f"⚠️ {ct}: 缺失基因 {missing}")
        if not markers_found:
            del_markers.append(ct)
    for ct in del_markers:
        del filtered_marker_dict[ct]
        if verbose:
            print(f"❌ 删除无有效marker的亚型: {ct}")
    return filtered_marker_dict

def plot_marker_dotplot(
    adata: AnnData,
    marker_dict: Dict[str, List[str]],
    groupby: str = "leiden_res1",
    save_path: Optional[str] = None,
    standard_scale: str = "var",
    dendrogram: bool = True,
    dpi: int = 300,
    bbox_inches: str = "tight"
) -> None:
    sc.pl.dotplot(
        adata,
        groupby=groupby,
        var_names=marker_dict,
        dendrogram=dendrogram,
        standard_scale=standard_scale,
        show=False
    )
    if save_path is not None:
        plt.savefig(save_path, dpi=dpi, bbox_inches=bbox_inches)
        plt.close()
        print(f"✅ 亚型marker dotplot已保存至: {save_path}")

def save_marker_genes(
    marker_dict: Dict[str, List[str]],
    save_path: str,
    sep: str = ","
) -> None:
    with open(save_path, 'w') as f:
        for cell_type, markers in marker_dict.items():
            f.write(f"{cell_type}: {sep.join(markers)}\n")
    print(f"✅ 过滤后亚型marker已保存至: {save_path}")

def compute_cluster_degs(
    adata: AnnData,
    groupby: str = "leiden_res1",
    method: str = "wilcoxon",
    key_added: str = "rank_genes_groups",
    use_raw: bool = False,
    min_in_group_fraction: float = 0.1,
    max_out_group_fraction: float = 0.2,
    filter_key_added: str = "rank_genes_groups_filtered"
) -> AnnData:
    if 'log1p' in adata.uns and adata.uns['log1p'].get('base') is not None:
        adata.uns['log1p']['base'] = None
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        key_added=key_added,
        use_raw=use_raw
    )
    sc.tl.filter_rank_genes_groups(
        adata,
        min_in_group_fraction=min_in_group_fraction,
        max_out_group_fraction=max_out_group_fraction,
        key=key_added,
        key_added=filter_key_added
    )
    print(f"✅ 亚型差异基因计算完成，结果存储在 adata.uns['{filter_key_added}']")
    return adata

def plot_deg_dotplot(
    adata: AnnData,
    groupby: str = "leiden_res1",
    key: str = "rank_genes_groups_filtered",
    n_genes: int = 3,
    save_path: Optional[str] = None,
    standard_scale: str = "var",
    dendrogram: bool = True,
    dpi: int = 300,
    bbox_inches: str = "tight"
) -> None:
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=groupby,
        key=key,
        n_genes=n_genes,
        dendrogram=dendrogram,
        standard_scale=standard_scale,
        show=False
    )
    if save_path is not None:
        plt.savefig(save_path, dpi=dpi, bbox_inches=bbox_inches)
        plt.close()
        print(f"✅ 亚型DEG dotplot已保存至: {save_path}")

def save_cluster_top_markers(
    adata: AnnData,
    save_path: str,
    groupby: str = "leiden_res1",
    key: str = "rank_genes_groups_filtered",
    top_n: int = 10
) -> None:
    clusters = sorted(adata.obs[groupby].cat.categories)
    with open(save_path, 'w') as f:
        for cluster in clusters:
            degs = sc.get.rank_genes_groups_df(
                adata,
                group=str(cluster),
                key=key
            ).dropna().head(top_n)
            f.write(f'=== Subtype {cluster} ===\n')
            f.write('\n'.join(degs['names'].astype(str)) + '\n\n')
    print(f"✅ 各亚型Top{top_n}标记基因已保存至: {save_path}")

# ===================== 【模块2：AUCell富集分析（源自AUCell4anno.py）】 =====================
def load_cellmarker_geneset(db_path: str, organism: str = "Human") -> dict:
    print(f"\n[AUCell] 加载CellMarker数据库: {db_path}")
    pathway_dict = ov.utils.geneset_prepare(
        db_path,
        organism=organism
    )
    print(f"✅ [AUCell] 加载 {len(pathway_dict)} 个细胞类型基因集")
    return pathway_dict

def run_aucell_enrichment(
    adata_path: str,
    pathway_dict: dict,
    num_workers: int = 1
) -> sc.AnnData:
    print(f"\n[AUCell] 加载亚型数据并执行AUCell富集")
    adata = ov.read(adata_path)
    adata_raw = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    adata_aucs = ov.single.pathway_aucell_enrichment(
        adata_raw,
        pathways_dict=pathway_dict,
        num_workers=num_workers,
    )
    adata_aucs.obs = adata[adata_aucs.obs.index].obs
    adata_aucs.obsm = adata[adata_aucs.obs.index].obsm
    adata_aucs.obsp = adata[adata_aucs.obs.index].obsp
    print(f"✅ [AUCell] 完成，结果维度: {adata_aucs.shape}")
    return adata_aucs

def run_aucell_diff_analysis(
    adata_aucs: sc.AnnData,
    groupby: str = "leiden_res1",
    method: str = "wilcoxon",
    key_added: str = "dea_leiden_aucs_res1"
) -> sc.AnnData:
    print(f"\n[AUCell] 执行亚型AUCell差异分析 (groupby={groupby})")
    if 'log1p' not in adata_aucs.uns:
        adata_aucs.uns['log1p'] = {'base': None}
    else:
        adata_aucs.uns['log1p']['base'] = None
    sc.tl.rank_genes_groups(
        adata_aucs,
        groupby=groupby,
        use_raw=False,
        method=method,
        key_added=key_added
    )
    print("✅ [AUCell] 差异分析完成")
    return adata_aucs

def plot_aucell_dotplot(
    adata_aucs: sc.AnnData,
    save_path: str,
    groupby: str = "leiden_res1",
    key: str = "dea_leiden_aucs_res1",
    n_genes: int = 3,
    cmap: str = "RdBu_r"
):
    print(f"\n[AUCell] 绘制亚型AUCell dotplot")
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
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ [AUCell] Dotplot已保存至: {save_path}")

def save_aucell_cluster_annotations(
    adata_aucs: sc.AnnData,
    save_path: str,
    groupby: str = "leiden_res1",
    key: str = "dea_leiden_aucs_res1"
):
    print(f"\n[AUCell] 保存亚型AUCell注释结果")
    clusters = (
        adata_aucs.uns[f'dendrogram_{groupby}']['categories_ordered']
        if f'dendrogram_{groupby}' in adata_aucs.uns
        else sorted(adata_aucs.obs[groupby].cat.categories)
    )
    with open(save_path, 'w') as f:
        for cluster in clusters:
            special_cluster = str(cluster)
            degs = sc.get.rank_genes_groups_df(
                adata_aucs,
                group=special_cluster,
                key=key
            ).dropna()
            f.write(f'{special_cluster}: {"|".join(degs.names[:2].tolist())}\n')
    print(f"✅ [AUCell] 亚型注释已保存至: {save_path}")

# ===================== 【核心主流程：亚型注释一键执行】 =====================
def subtype_anno_pipeline(
    adata_path: str,
    subtype_marker_path: str,
    cellmarker_db_path: str,
    root_output_dir: str,
    major_cell_type: str = "macrophage",
    groupby: str = "leiden_res1",
    organism: str = "Human",
    num_workers: int = 1,
    deg_method: str = "wilcoxon",
    top_n_degs: int = 10,
    dpi: int = 300
):
    """
    单细胞亚型注释核心流程（marker+AUCell双流程）
    :param adata_path: 提取后的主细胞类型h5ad路径
    :param subtype_marker_path: 自定义亚型marker字典路径
    :param cellmarker_db_path: CellMarker数据库路径
    :param root_output_dir: 根输出目录
    :param major_cell_type: 主细胞类型名（用于生成专属结果文件夹）
    :param groupby: 亚型聚类列名
    :param organism: 物种（Human/Mouse）
    :param num_workers: AUCell线程数
    :param deg_method: 差异分析方法
    :param top_n_degs: 每个亚型保存的top差异基因数
    :param dpi: 图片分辨率
    """
    # 1. 生成主细胞类型专属结果文件夹
    sub_output_dir = os.path.join(root_output_dir, major_cell_type)
    os.makedirs(sub_output_dir, exist_ok=True)
    print(f"📁 亚型分析结果文件夹已创建: {sub_output_dir}")

    # 2. 流程1：自定义亚型marker基因分析
    print("\n" + "="*50)
    print("【开始：自定义亚型marker基因分析】")
    print("="*50)
    adata = ov.read(adata_path)
    # 2.1 加载并过滤亚型marker
    marker_dict = load_marker_dict(subtype_marker_path)
    filtered_markers = check_marker_genes(adata, marker_dict)
    save_marker_genes(filtered_markers, os.path.join(sub_output_dir, f"{major_cell_type}-filtered_subtype_markers.txt"))
    # 2.2 绘制亚型marker dotplot
    plot_marker_dotplot(
        adata, filtered_markers, groupby=groupby,
        save_path=os.path.join(sub_output_dir, f"{major_cell_type}-subtype_marker_dotplot.png"),
        dpi=dpi
    )
    # 2.3 计算亚型差异基因
    adata = compute_cluster_degs(
        adata, groupby=groupby, method=deg_method,
        key_added=f"dea_{major_cell_type}_{groupby}",
        filter_key_added=f"dea_{major_cell_type}_{groupby}_filtered"
    )
    # 2.4 绘制亚型DEG dotplot
    plot_deg_dotplot(
        adata, groupby=groupby,
        key=f"dea_{major_cell_type}_{groupby}_filtered",
        save_path=os.path.join(sub_output_dir, f"{major_cell_type}-subtype_deg_dotplot.png"),
        dpi=dpi
    )
    # 2.5 保存亚型top差异基因
    save_cluster_top_markers(
        adata, groupby=groupby,
        key=f"dea_{major_cell_type}_{groupby}_filtered",
        top_n=top_n_degs,
        save_path=os.path.join(sub_output_dir, f"{major_cell_type}-subtype_top{top_n_degs}_markers.txt")
    )

 

    keys_to_del = [
        f"dea_{major_cell_type}_{groupby}",
        f"dea_{major_cell_type}_{groupby}_filtered"
    ]
    for k in keys_to_del:
        if k in adata.uns:
            del adata.uns[k]

    # 保存
    marker_adata_path = os.path.join(sub_output_dir, f"{major_cell_type}-marker_annotated_subtype.h5ad")
    try:
        adata.write(marker_adata_path)
        print(f"✅ marker分析后的adata已保存至: {marker_adata_path}")
    except Exception as e:
        print(f"⚠️ 跳过保存adata（低版本兼容）: {e}")

    # 3. 流程2：CellMarker数据库AUCell富集分析
    print("\n" + "="*50)
    print("【开始：AUCell数据库富集亚型分析】")
    print("="*50)
    # 3.1 加载CellMarker基因集
    pathway_dict = load_cellmarker_geneset(cellmarker_db_path, organism)
    # 3.2 执行AUCell富集
    adata_aucs = run_aucell_enrichment(adata_path, pathway_dict, num_workers)
    # 3.3 AUCell差异分析
    adata_aucs = run_aucell_diff_analysis(adata_aucs, groupby=groupby, key_added=f"dea_{major_cell_type}_aucs_{groupby}")
    # 3.4 绘制AUCell dotplot
    plot_aucell_dotplot(
        adata_aucs, groupby=groupby,
        key=f"dea_{major_cell_type}_aucs_{groupby}",
        save_path=os.path.join(sub_output_dir, f"{major_cell_type}-aucell_subtype_dotplot.png")
    )
    # 3.5 保存AUCell亚型注释结果
    save_aucell_cluster_annotations(
        adata_aucs, groupby=groupby,
        key=f"dea_{major_cell_type}_aucs_{groupby}",
        save_path=os.path.join(sub_output_dir, f"{major_cell_type}-aucell_subtype_annotations.txt")
    )

    # AUCell adata 保存
    aucell_adata_path = os.path.join(sub_output_dir, f"{major_cell_type}-aucell_annotated_subtype.h5ad")
    try:
        adata_aucs.write(aucell_adata_path)
        print(f"✅ AUCell分析后的adata已保存至: {aucell_adata_path}")
    except Exception as e:
        print(f"⚠️ 跳过保存AUCell adata: {e}")

    # 4. 流程结束
    print("\n🎉 【亚型注释全流程执行完成】")
    print(f"📁 所有结果已保存至: {sub_output_dir}")
    return adata, adata_aucs

# 主函数（命令行调用）
if __name__ == "__main__":
    subtype_anno_pipeline(
        adata_path="./sc_subtypeProc-macrophage.h5ad",
        subtype_marker_path="./sc_marker4anno-subtype_marker_dict.txt",
        cellmarker_db_path="./CellMarker_Augmented_2021.txt",
        root_output_dir="./subtype_annotation_results",
        major_cell_type="macrophage",
        groupby="leiden_res1_0",
        organism="Human",
        num_workers=1,
        deg_method="wilcoxon",
        top_n_degs=10,
        dpi=300
    )
