import scanpy as sc
import pandas as pd
from sctcm.sc.io import read_and_merge_samples

def load_homolog_mapping(homolog_file: str) -> dict:
    """
    加载小鼠-人类基因同源映射表
    须包含列：MGI.symbol(小鼠基因), HGNC.symbol(人类基因)
    """
    homolog_df = pd.read_csv(homolog_file)
    mouse_to_human = homolog_df.drop_duplicates("MGI.symbol").set_index("MGI.symbol")["HGNC.symbol"].to_dict()
    print(f"✅ 同源映射表加载完成：{len(mouse_to_human)} 个基因对应关系")
    return mouse_to_human

def mouse2human_convert(
    adata: sc.AnnData,
    mouse_to_human: dict,
    min_cells: int = 1,
    keep_duplicates: str = "first"
) -> sc.AnnData:
    """
    小鼠基因 → 人类基因 转换主函数
    """
    print(f"原始数据：{adata.shape[0]} 个细胞，{adata.shape[1]} 个小鼠基因")

    # 1. 基因名映射
    adata.var["human_gene"] = [mouse_to_human.get(gene, None) for gene in adata.var_names]
    
    # 2. 过滤无同源基因
    adata_human = adata[:, ~adata.var["human_gene"].isna()].copy()
    
    # 3. 替换为人类基因名
    adata_human.var_names = adata_human.var["human_gene"].astype(str)
    del adata_human.var["human_gene"]
    adata_human.var.index.name = None

    # 4. 去重基因名
    adata_human = adata_human[:, ~adata_human.var_names.duplicated(keep=keep_duplicates)].copy()

    # 5. 过滤低表达基因
    sc.pp.filter_genes(adata_human, min_cells=min_cells)

    print(f"转换完成：{adata_human.shape[1]} 个人类基因保留")
    return adata_human

def process_mouse_sc_data(
    sample_table: str,
    homolog_file: str,
    output_h5ad: str,
    sep: str = "\t",
    min_cells: int = 1
) -> None:
    """
    一站式：读取合并样本 → 小鼠转人基因 → 保存h5ad
    """
    # 读取合并样本
    print("🔽 开始读取并合并样本...")
    adata_merged = read_and_merge_samples(sample_table=sample_table, sep=sep)
    
    # 加载映射
    mouse_to_human = load_homolog_mapping(homolog_file)
    
    # 基因转换
    print("🔁 开始小鼠基因 → 人类基因转换...")
    adata_human = mouse2human_convert(adata_merged, mouse_to_human, min_cells)
    
    # 保存结果
    print(f"💾 保存结果到：{output_h5ad}")
    sc.write(output_h5ad, adata_human)
    print("🎉 全部任务完成！")

if __name__ == "__main__":
    import os
    # 默认测试路径
    SAMPLE_TEMPLATE = "/scTCM/docs/templates/mouse_sc_input_template.txt"
    HOMOLOG_FILE = "/scTCM/docs/templates/mouse_human_homologs_dec2021.csv"
    OUTPUT_DIR = "/scTCM/sctcm/tests/result"
    OUTPUT_H5AD = os.path.join(OUTPUT_DIR, "mouse2human_converted.h5ad")

    process_mouse_sc_data(SAMPLE_TEMPLATE, HOMOLOG_FILE, OUTPUT_H5AD)
