# -*- coding: utf-8 -*-
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from pathlib import Path
import matplotlib.pyplot as plt
import sys
import warnings
warnings.filterwarnings('ignore')


def TCMscoreScan_case_high(
    adata_path,
    output_dir,
    group_col="Group",
    case_label="CASE",
    ctrl_label="CTRL",
    celltype_col="celltype_manual",
    pval_cutoff=0.05,
    fc_cutoff=1.5
):
    """
    仅绘制 CASE 显著高于 CTRL 的 TCMscore 细胞类型点图
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 读取数据
    adata = sc.read_h5ad(adata_path)
    print(f"✅ 成功读取 AnnData：{adata_path}")
    print(f"✅ 细胞数量：{adata.n_obs} | 基因数量：{adata.n_vars}")

    # 自动识别所有 TCMscore
    tcm_scores = [col for col in adata.obs.columns if "TCMscore_" in col]
    if len(tcm_scores) == 0:
        print("❌ 未找到任何 TCMscore 列")
        return

    print(f"✅ 共找到 {len(tcm_scores)} 个 TCMscore")

    # 差异分析
    res = []
    for score in tcm_scores:
        try:
            case = adata[adata.obs[group_col] == case_label].obs[score].values
            ctrl = adata[adata.obs[group_col] == ctrl_label].obs[score].values

            stat, pval = mannwhitneyu(case, ctrl, alternative="greater")
            mean_case = np.mean(case)
            mean_ctrl = np.mean(ctrl)
            fc = mean_case / (mean_ctrl + 1e-8)

            res.append({
                "score": score,
                "mean_case": mean_case,
                "mean_ctrl": mean_ctrl,
                "fc": fc,
                "pval": pval
            })
        except Exception as e:
            continue

    # 筛选
    df = pd.DataFrame(res)
    df = df[(df["pval"] < pval_cutoff) & (df["fc"] > fc_cutoff)]
    keep = df["score"].tolist()

    print(f"✅ 筛选后 CASE 显著高的条目：{len(keep)} 个")

    if len(keep) == 0:
        print("❌ 无满足条件的条目")
        return

    # 保存筛选结果
    df.to_csv(out_dir / "d4c_TCMscoreScan_case_high_significant.csv", index=False)

    # 绘图
    sc.pl.dotplot(
        adata,
        var_names=keep,
        groupby=celltype_col,
        show=False,
        standard_scale="var",
        title="TCMscore (CASE > CTRL)"
    )
    plt.savefig(
        out_dir / "d4c_TCMscoreScan-case_high_only_dotplot.pdf",
        bbox_inches="tight"
    )
    plt.close()

    print(f"🎉 绘图完成，已保存至：{output_dir}")
    return df
