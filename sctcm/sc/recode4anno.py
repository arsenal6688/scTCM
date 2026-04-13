import scanpy as sc
import omicverse as ov
import matplotlib.pyplot as plt
import pickle
from typing import Dict, Optional, Union, List
from pathlib import Path

def read_recode_template(template_path: Union[str, Path]) -> Dict[str, str]:
    """
    读取cluster注释模板：每行是 数字[TAB]细胞类型
    跳过空行、兼容TAB/空格分隔
    """
    cluster2annotation = {}
    with open(template_path, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue

            # 按【制表符TAB】或【空格】分割
            parts = line.split('\t')  # 优先TAB
            if len(parts) == 1:
                parts = line.split()   # 兼容空格

            if len(parts) >= 2:
                cluster = parts[0].strip()
                celltype = parts[1].strip()
                cluster2annotation[cluster] = celltype
            else:
                print(f"⚠️  跳过无效行 {line_num}: {line}")

    return cluster2annotation

def recode_cluster_annotation(
    adata: sc.AnnData,
    cluster_key: str = "leiden",
    annotation_template_path: Union[str, Path] = None,
    new_anno_key: str = "celltype_manual"
) -> sc.AnnData:
    """
    核心：根据模板文件对 cluster 进行手动重新注释
    """
    if annotation_template_path is None:
        raise ValueError("必须提供注释模板路径")

    annotation_dict = read_recode_template(annotation_template_path)

    # 执行注释
    adata.obs[new_anno_key] = adata.obs[cluster_key].astype(str).map(annotation_dict)
    adata.obs[new_anno_key] = adata.obs[new_anno_key].astype('category')

    return adata

def plot_umap(
    adata: sc.AnnData,
    color_key: str = "celltype_manual",
    save_path: Union[str, Path] = None,
    title: str = "Manual Annotated Cell Types (UMAP)",
    palette: Optional[List[str]] = None,
    dpi: int = 300,
    frameon: bool = True,
    **kwargs
):
    if palette is None:
        palette = ov.palette()

    sc.pl.umap(
        adata,
        color=color_key,
        title=title,
        palette=palette,
        frameon=frameon,
        show=False,
        **kwargs
    )

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
    plt.close()

def plot_mde(
    adata: sc.AnnData,
    color_keys: List[str] = ["leiden", "celltype_manual"],
    titles: List[str] = ["Clusters", "Manual Annotated Cell Types"],
    save_path: Union[str, Path] = None,
    palette: Optional[List[str]] = None,
    dpi: int = 300,
    **kwargs
):
    if palette is None:
        palette = ov.palette()

    ov.utils.embedding(
        adata,
        basis="X_mde",
        color=color_keys,
        title=titles,
        palette=palette,
        show=False,
        **kwargs
    )

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
    plt.close()

def save_adata_safely(
    adata: sc.AnnData,
    save_path: Union[str, Path],
    pkl_backup_path: Optional[Union[str, Path]] = None
):
    complex_keys = ["dea_leiden_filtered", "dea_results"]
    backup = {}

    for key in complex_keys:
        if key in adata.uns:
            backup[key] = adata.uns.pop(key)

    if pkl_backup_path and backup:
        with open(pkl_backup_path, 'wb') as f:
            pickle.dump(backup, f)

    adata.write(save_path)

def run_annotation_pipeline(
    adata_path: Union[str, Path],
    template_path: Union[str, Path],
    output_dir: Union[str, Path],
    cluster_key: str = "leiden",
    new_anno_key: str = "celltype_manual",
    run_umap: bool = True,
    run_mde: bool = True
):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    adata = sc.read_h5ad(adata_path)

    adata = recode_cluster_annotation(
        adata=adata,
        cluster_key=cluster_key,
        annotation_template_path=template_path,
        new_anno_key=new_anno_key
    )

    if run_mde and "X_mde" in adata.obsm:
        mde_path = output_dir / "sc_recode4anno-manual_annotation_mde.png"
        plot_mde(
            adata=adata,
            color_keys=[cluster_key, new_anno_key],
            save_path=mde_path
        )

    if run_umap:
        if "X_umap" not in adata.obsm:
            sc.tl.umap(adata)
        umap_path = output_dir / "sc_recode4anno-manual_annotation_umap.svg"
        plot_umap(
            adata=adata,
            color_key=new_anno_key,
            save_path=umap_path
        )

    adata_save = output_dir / "sc_recode4anno-celltype_manual_annotated.h5ad"
    pkl_save = output_dir / "sc_recode4anno-backup_uns.pkl"
    save_adata_safely(adata, save_path=adata_save, pkl_backup_path=pkl_save)

    return adata

__all__ = [
    "read_recode_template",
    "recode_cluster_annotation",
    "plot_umap",
    "plot_mde",
    "save_adata_safely",
    "run_annotation_pipeline"
]
