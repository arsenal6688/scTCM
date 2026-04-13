# -*- coding: utf-8 -*-
import os
import pickle
import scanpy as sc
import drug2cell as d2c
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")


def run_chembl_drug2cell(h5ad_path, drug_dict_pkl, output_dir, groupby="cell_type_plot"):
    """
    运行 ChEMBL36 + drug2cell 细胞类型药物打分
    原版流程：score → rank_genes_groups → dotplot
    """
    os.makedirs(output_dir, exist_ok=True)

    # -------------------- 1. 读取数据 --------------------
    print("[1/5] 读取 AnnData ...")
    adata = sc.read_h5ad(h5ad_path)

    # -------------------- 2. 标准化 --------------------
    print("[2/5] 数据标准化 ...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    # -------------------- 3. 加载药物 --------------------
    print("[3/5] 加载 ChEMBL 药物集 ...")
    with open(drug_dict_pkl, "rb") as f:
        drug_dict = pickle.load(f)

    # -------------------- 4. drug2cell 打分 --------------------
    print("[4/5] 运行 drug2cell.score ...")
    d2c.score(
        adata,
        targets=drug_dict,
        nested=True,
        use_raw=True
    )

    # -------------------- 5. 差异分析 --------------------
    print("[5/5] 差异分析 & 绘图 ...")
    sc.tl.rank_genes_groups(
        adata.uns["drug2cell"],
        groupby=groupby,
        method="wilcoxon"
    )

    # -------------------- 绘图：通用版 --------------------
    sc.pl.rank_genes_groups_dotplot(
        adata.uns["drug2cell"],
        groupby=groupby,
        swap_axes=True,
        dendrogram=False,
        n_genes=10,        # TOP10
        cmap="Blues",
        vmax=0.6,
        standard_scale="var",  
        show=False

    )
    
    plt.savefig(
        os.path.join(output_dir, "dotplot_fig2k.pdf"),
        bbox_inches="tight"
    )
    plt.close()

    # -------------------- 保存差异结果 --------------------
    res = adata.uns["drug2cell"].uns["rank_genes_groups"]
    groups = res["names"].dtype.names
    deg_df = pd.DataFrame({
        f"{g}_{k[:1]}": res[k][g]
        for g in groups
        for k in ["names", "scores", "pvals_adj", "logfoldchanges"]
    })
    deg_df.to_csv(os.path.join(output_dir, "chembl_degs.csv"), index=False)

    print(f"\n✅ 完成！结果已保存至：{output_dir}")
    return adata
