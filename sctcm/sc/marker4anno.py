"""
sctcm.marker4anno
-----------------
单细胞数据marker基因分析与注释工具
包含marker基因加载、过滤、差异基因计算、可视化等功能
"""

import os
import pickle
from typing import Dict, List, Optional, Union, Tuple

import scanpy as sc
import omicverse as ov
import pandas as pd
import matplotlib.pyplot as plt
from anndata import AnnData

# 设置全局绘图样式
ov.ov_plot_set()

__all__ = [
    "load_marker_dict",
    "check_marker_genes",
    "plot_marker_dotplot",
    "sdea_leidenave_marker_genes",
    "compute_cluster_degs",
    "plot_deg_dotplot",
    "save_cluster_top_markers",
    "marker_analysis_pipeline"
]


def load_marker_dict(
    marker_file_path: str,
    cell_type_col: Union[int, str] = 0,
    genes_col: Union[int, str] = 1,
    gene_sep: str = "|",
    encoding: str = "utf-8"
) -> Dict[str, List[str]]:
    """
    从txt文件加载marker字典（格式：CellType\tGene1|Gene2|Gene3）
    
    Parameters
    ----------
    marker_file_path
        marker文件路径
    cell_type_col
        细胞类型所在列索引/列名（默认0列）
    genes_col
        基因列表所在列索引/列名（默认1列）
    gene_sep
        基因间分隔符（默认|）
    encoding
        文件编码格式
    
    Returns
    -------
    marker_dict
        细胞类型-标记基因字典
    """
    marker_dict = {}
    if not os.path.exists(marker_file_path):
        raise FileNotFoundError(f"Marker文件不存在: {marker_file_path}")
    
    with open(marker_file_path, 'r', encoding=encoding) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue  # 跳过空行
            
            # 按制表符分割行
            parts = line.split('\t')
            if len(parts) < 2:
                print(f"⚠️ 第{line_num}行格式错误，跳过: {line}")
                continue
            
            cell_type = parts[cell_type_col].strip()
            genes_str = parts[genes_col].strip()
            
            # 分割基因并去重、过滤空值
            genes = [gene.strip() for gene in genes_str.split(gene_sep) if gene.strip()]
            
            if cell_type and genes:
                marker_dict[cell_type] = genes
            else:
                print(f"⚠️ 第{line_num}行数据无效，跳过: {line}")
    
    if not marker_dict:
        raise ValueError("未从文件中加载到有效marker数据")
    
    print(f"✅ 成功加载{len(marker_dict)}种细胞类型的marker基因")
    return marker_dict


def check_marker_genes(
    adata: AnnData,
    marker_dict: Dict[str, List[str]],
    verbose: bool = True
) -> Dict[str, List[str]]:
    """
    检查输入的marker基因是否存在于adata的var索引中，过滤不存在的基因
    
    Parameters
    ----------
    adata
        聚类后的AnnData对象
    marker_dict
        细胞类型-标记基因字典，格式: {cell_type: [gene1, gene2,...]}
    verbose
        是否打印过滤信息
    
    Returns
    -------
    filtered_marker_dict
        过滤后的标记基因字典（仅保留存在于adata中的基因）
    """
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
    
    # 删除无有效marker的细胞类型
    for ct in del_markers:
        del filtered_marker_dict[ct]
        if verbose:
            print(f"❌ 删除无有效marker的细胞类型: {ct}")
    
    return filtered_marker_dict


def plot_marker_dotplot(
    adata: AnnData,
    marker_dict: Dict[str, List[str]],
    groupby: str = "leiden",
    save_path: Optional[str] = None,
    standard_scale: str = "var",
    dendrogram: bool = True,
    dpi: int = 300,
    bbox_inches: str = "tight"
) -> None:
    """
    根据已知marker基因绘制dotplot
    
    Parameters
    ----------
    adata
        聚类后的AnnData对象
    marker_dict
        细胞类型-标记基因字典（建议先通过check_marker_genes过滤）
    groupby
        分组依据（如leiden聚类结果列名）
    save_path
        图片保存路径，None则不保存
    standard_scale
        标准化方式，"var"按基因标准化，"group"按分组标准化
    dendrogram
        是否显示聚类树
    dpi
        图片分辨率
    bbox_inches
        图片边界设置
    """
    # 绘制dotplot
    sc.pl.dotplot(
        adata,
        groupby=groupby,
        var_names=marker_dict,
        dendrogram=dendrogram,
        standard_scale=standard_scale,
        show=False
    )
    
    # 保存图片
    if save_path is not None:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=dpi, bbox_inches=bbox_inches)
        plt.close()
        print(f"✅ 已知marker dotplot已保存至: {save_path}")


def save_marker_genes(
    marker_dict: Dict[str, List[str]],
    save_path: str,
    sep: str = ","
) -> None:
    """
    保存过滤后的marker基因到文本文件
    
    Parameters
    ----------
    marker_dict
        过滤后的标记基因字典
    save_path
        保存路径
    sep
        基因间分隔符
    """
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    with open(save_path, 'w') as f:
        for cell_type, markers in marker_dict.items():
            f.write(f"{cell_type}: {sep.join(markers)}\n")
    print(f"✅ 标记基因列表已保存至: {save_path}")


def compute_cluster_degs(
    adata: AnnData,
    groupby: str = "leiden",
    method: str = "wilcoxon",
    key_added: str = "rank_genes_groups",
    use_raw: bool = False,
    min_in_group_fraction: float = 0.1,
    max_out_group_fraction: float = 0.2,
    filter_key_added: str = "rank_genes_groups_filtered"
) -> AnnData:
    """
    计算每个cluster的差异表达基因，并过滤低可信度的DEGs
    
    Parameters
    ----------
    adata
        聚类后的AnnData对象
    groupby
        分组依据（聚类结果列名）
    method
        DEG计算方法（wilcoxon/t-test等）
    key_added
        DEG结果存储的key
    use_raw
        是否使用原始表达矩阵
    min_in_group_fraction
        基因在组内表达的最小细胞比例
    max_out_group_fraction
        基因在组外表达的最大细胞比例
    filter_key_added
        过滤后DEG结果存储的key
    
    Returns
    -------
    adata
        更新后的AnnData对象（包含DEG计算结果）
    """
    # 重置log1p基数
    if 'log1p' in adata.uns and adata.uns['log1p'].get('base') is not None:
        adata.uns['log1p']['base'] = None
    
    # 计算差异基因
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        key_added=key_added,
        use_raw=use_raw
    )
    
    # 过滤差异基因
    sc.tl.filter_rank_genes_groups(
        adata,
        min_in_group_fraction=min_in_group_fraction,
        max_out_group_fraction=max_out_group_fraction,
        key=key_added,
        key_added=filter_key_added
    )
    
    print(f"✅ 差异基因计算完成，结果存储在 adata.uns['{filter_key_added}']")
    return adata


def plot_deg_dotplot(
    adata: AnnData,
    groupby: str = "leiden",
    key: str = "rank_genes_groups_filtered",
    n_genes: int = 3,
    save_path: Optional[str] = None,
    standard_scale: str = "var",
    dendrogram: bool = True,
    dpi: int = 300,
    bbox_inches: str = "tight"
) -> None:
    """
    绘制差异基因dotplot
    
    Parameters
    ----------
    adata
        包含DEG结果的AnnData对象
    groupby
        分组依据
    key
        DEG结果的key
    n_genes
        每个cluster展示的top基因数
    save_path
        图片保存路径
    standard_scale
        标准化方式
    dendrogram
        是否显示聚类树
    dpi
        图片分辨率
    bbox_inches
        图片边界设置
    """
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
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=dpi, bbox_inches=bbox_inches)
        plt.close()
        print(f"✅ 差异基因dotplot已保存至: {save_path}")


def save_cluster_top_markers(
    adata: AnnData,
    save_path: str,
    groupby: str = "leiden",
    key: str = "rank_genes_groups_filtered",
    top_n: int = 10
) -> None:
    """
    保存每个cluster的top N差异标记基因
    
    Parameters
    ----------
    adata
        包含DEG结果的AnnData对象
    groupby
        分组依据
    key
        DEG结果的key
    top_n
        每个cluster保存的top基因数
    save_path
        保存路径
    """
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    # 获取排序后的cluster列表
    clusters = sorted(adata.obs[groupby].cat.categories)
    
    with open(save_path, 'w') as f:
        for cluster in clusters:
            # 提取当前cluster的top N DEGs
            degs = sc.get.rank_genes_groups_df(
                adata,
                group=str(cluster),
                key=key
            ).dropna().head(top_n)
            
            # 写入文件
            f.write(f'=== Cluster {cluster} ===\n')
            f.write('\n'.join(degs['names'].astype(str)) + '\n\n')
    
    print(f"✅ 各Cluster Top{top_n}标记基因已保存至: {save_path}")


def marker_analysis_pipeline(
    adata_path: str,
    marker_file_path: str,
    output_dir: str,
    groupby: str = "leiden_res1",
    deg_method: str = "wilcoxon",
    top_n_degs: int = 10,
    dpi: int = 300
) -> AnnData:
    """
    整合已知marker分析 + 差异marker计算的完整流程（从txt加载marker）
    
    Parameters
    ----------
    adata_path
        聚类后的h5ad文件路径
    marker_file_path
        marker字典txt文件路径（格式：CellType\tGene1|Gene2|Gene3）
    output_dir
        结果输出目录
    groupby
        分组依据（聚类结果列名）
    deg_method
        DEG计算方法
    top_n_degs
        每个cluster保存的top DEG数
    dpi
        图片分辨率
    
    Returns
    -------
    adata
        处理后的AnnData对象
    """
    # 1. 读取数据
    adata = ov.read(adata_path)
    print(f"✅ 成功读取数据: {adata_path}")
    
    # 2. 加载并过滤marker基因
    # 2.1 从txt加载marker字典
    marker_dict = load_marker_dict(marker_file_path)
    # 2.2 过滤marker基因
    filtered_markers = check_marker_genes(adata, marker_dict)
    # 2.3 保存过滤后的marker
    marker_save_path = os.path.join(output_dir, "sc_marker4anno-filtered_marker_genes.txt")
    save_marker_genes(filtered_markers, marker_save_path)
    # 2.4 绘制已知marker dotplot
    marker_dotplot_path = os.path.join(output_dir, "sc_marker4anno-marker_dotplot.png")
    plot_marker_dotplot(
        adata,
        filtered_markers,
        groupby=groupby,
        save_path=marker_dotplot_path,
        dpi=dpi
    )
    
    # 3. 差异marker计算
    # 3.1 计算DEGs
    adata = compute_cluster_degs(
        adata,
        groupby=groupby,
        method=deg_method,
        key_added=f"dea_{groupby}",
        filter_key_added=f"dea_{groupby}_filtered"
    )
    # 3.2 绘制DEG dotplot
    deg_dotplot_path = os.path.join(output_dir, "sc_marker4anno-deg_dotplot.png")
    plot_deg_dotplot(
        adata,
        groupby=groupby,
        key=f"dea_{groupby}_filtered",
        save_path=deg_dotplot_path,
        dpi=dpi
    )
    # 3.3 保存top DEGs
    top_marker_save_path = os.path.join(output_dir, f"sc_marker4anno-all_clusters_top{top_n_degs}_markers.txt")
    save_cluster_top_markers(
        adata,
        save_path=top_marker_save_path,
        groupby=groupby,
        key=f"dea_{groupby}_filtered",
        top_n=top_n_degs
    )
    
    # 4. 保存DEG结果
    deg_key = f"dea_{groupby}_filtered"
    if deg_key in adata.uns:
        deg_pkl_path = os.path.join(output_dir, f"sc_marker4anno-{deg_key}.pkl")
        with open(deg_pkl_path, 'wb') as f:
            pickle.dump(adata.uns[deg_key], f)
        # 移除复杂结构以支持h5ad保存
        adata.uns.pop(deg_key)
        print(f"✅ DEG结果已保存至: {deg_pkl_path}")
    
    # 5. 保存处理后的adata
    adata_save_path = os.path.join(output_dir, "sc_marker4anno-marker_annotated.h5ad")
    adata.write(adata_save_path)
    print(f"✅ 处理后的数据已保存至: {adata_save_path}")

    return adata
