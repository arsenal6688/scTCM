# sctcm/sc/io.py
import os
import anndata as ad
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")

# ----------------------
# 1. 读取 10x .h5 文件
# ----------------------
def read_10x_h5(file_path: str) -> ad.AnnData:
    """读取 10x Genomics .h5 格式数据"""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"文件不存在: {file_path}")
    adata = sc.read_10x_h5(file_path)
    adata.var_names_make_unique()
    return adata

# ----------------------
# 2. 读取 10x mtx 文件夹
# ----------------------
def read_10x_mtx(folder_path: str) -> ad.AnnData:
    """读取 10x Genomics mtx 文件夹"""
    if not os.path.isdir(folder_path):
        raise NotADirectoryError(f"文件夹不存在: {folder_path}")
    adata = sc.read_10x_mtx(folder_path, var_names="gene_symbols")
    adata.var_names_make_unique()
    return adata

# ----------------------
# 3. 读取模板 + 分组 + 整合
# ----------------------
def read_and_merge_samples(
    sample_table: str,
    sep: str = "\t"
) -> ad.AnnData:
    """
    输入模板格式：Sample_name\tPath\tType\tGroup
    Type 支持：10x_h5, 10x_mtx
    整合方式完全对齐你提供的 sc.concat 标准写法
    """
    df = pd.read_csv(sample_table, sep=sep, dtype=str)
    required = ["Sample_name", "Path", "Type", "Group"]
    if not all(c in df.columns for c in required):
        raise ValueError(f"表格必须包含: {required}")

    adata_list = []
    sample_names = []

    # 逐个读取样本
    for _, row in df.iterrows():
        name = row["Sample_name"].strip()
        path = row["Path"].strip()
        typ = row["Type"].strip().lower()
        group = row["Group"].strip()

        print(f"▶ 读取样本: {name} | 类型: {typ} | 分组: {group}")

        # 读取数据
        if typ == "10x_h5":
            adata = read_10x_h5(path)
        elif typ == "10x_mtx":
            adata = read_10x_mtx(path)
        else:
            raise ValueError("仅支持 10x_h5 / 10x_mtx")

        # 标注分组
        adata.obs["Group"] = group
        adata.obs["Batch"] = name
        adata_list.append(adata)
        sample_names.append(name)

    print("\n✅ 开始合并所有样本...")
    adata_merged = sc.concat(
        adata_list,
        axis=0,
        join='outer',
        merge='same',
        index_unique="-",    # 细胞唯一标签（你要的！）
        keys=sample_names    # 样本名作为key
    )

    # 去重
    adata_merged.var_names_make_unique()
    adata_merged.obs_names_make_unique()

    # 输出统计
    print_stats(adata_merged)

    return adata_merged

# ----------------------
# 4. 统计信息（总细胞、总基因、分组）
# ----------------------
def print_stats(adata: ad.AnnData):
    print("\n" + "=" * 55)
    print("📊 数据整合完成（细胞ID已自动添加唯一标签）")
    print(f"总细胞数: {adata.n_obs}")
    print(f"总基因数: {adata.n_vars}")

    if "Group" in adata.obs:
        print("\n📦 分组统计:")
        print(adata.obs["Group"].value_counts())

    if "Batch" in adata.obs:
        print("\n📦 样本统计:")
        print(adata.obs["Batch"].value_counts())

    print("=" * 55 + "\n")

# ----------------------
# 5. 保存 h5ad
# ----------------------
def write_h5ad(adata: ad.AnnData, out_path: str):
    adata.write_h5ad(out_path, compression="gzip")
    print(f"✅ 数据已保存: {out_path}")


