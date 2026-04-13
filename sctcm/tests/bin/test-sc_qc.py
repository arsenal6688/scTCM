import anndata as ad
from sctcm.sc.qc import calculate_qc_metrics, filter_cells, filter_genes

# 配置参数
INPUT_H5AD = "/scTCM/sctcm/tests/result/sc_io-merged_samples.h5ad"
OUTPUT_FILE = "/scTCM/sctcm/tests/result/sc_io-merged_samples_qc.h5ad"
MITO_PATTERN = "^MT-"

# 加载数据
adata = ad.read_h5ad(INPUT_H5AD)

# 计算QC指标
adata_qc = calculate_qc_metrics(adata, mito_pattern=MITO_PATTERN)

# 细胞过滤
adata_cell_filtered, cell_filter_stats = filter_cells(
    adata_qc,
    min_genes=200,
    max_genes=5000,
    max_mito=10.0,
    min_umis=500
)

# 基因过滤
adata_final, gene_filter_stats = filter_genes(adata_cell_filtered, min_cells=3)

# 保存结果
adata_final.write_h5ad(OUTPUT_FILE)
