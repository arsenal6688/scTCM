# -*- coding: utf-8 -*-
import scanpy as sc
import drug2cell as d2c
import pandas as pd
import pickle
import warnings
from pathlib import Path
import numpy as np
from scipy.sparse import issparse
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sc.settings.set_figure_params(dpi=300, fontsize=10, facecolor="white", figsize=(6, 4))


def load_formula_targets(pkl_path):
    """
    加载四层嵌套复方PKL：
    格式：{
        "复方名": {
            "草药名": {
                "CID": [靶点1, 靶点2...],
                ...
            },
            ...
        },
        ...
    }
    合并为：复方 → 总靶点列表（去重）
    """
    with open(pkl_path, 'rb') as f:
        formula_dict = pickle.load(f)

    formula2targets = {}
    for formula_name, herb_dict in formula_dict.items():
        formula_name = str(formula_name).strip()
        all_targets = set()

        # 遍历草药
        for herb_name, cid_dict in herb_dict.items():
            # 遍历成分CID
            for cid, targets in cid_dict.items():
                if isinstance(targets, list) and len(targets) > 0:
                    all_targets.update([str(t).strip() for t in targets])

        target_list = sorted(list(all_targets))
        if len(target_list) > 0:
            formula2targets[formula_name] = target_list

    print(f"✅ 加载全量复方靶点：{len(formula2targets)} 个有效复方")
    return formula2targets


def formulaScan(
    sc_template,
    adata_path,
    formula_pkl,
    output_dir,
    min_pct: float = 0.25,
    logfc_threshold: float = 2.0,
    use_filtered_only: bool = False,
    celltype_col: str = "celltype_manual",
    group_col: str = "Group",
    case_name: str = "CASE",
    ctrl_name: str = "CTRL"
):
    """
    全复方批量drug2cell分析
    筛选标准：
    1. logFC > 阈值
    2. 得分在【细胞类型内部】至少 min_pct 细胞表达
    """
    print("\n" + "=" * 80)
    print("🚀 formulaScan - 中药全复方drug2cell批量扫描")
    print("=" * 80)

    print("\n📋 运行参数")
    print(f"• 分组文件: {sc_template}")
    print(f"• 单细胞数据: {adata_path}")
    print(f"• 全复方靶点: {formula_pkl}")
    print(f"• 输出目录: {output_dir}")
    print(f"• min_pct(细胞类型内): {min_pct} | logfc_threshold: {logfc_threshold}")
    print(f"• 严格过滤绘图: {use_filtered_only}")
    print(f"• 分组: {case_name} vs {ctrl_name}")

    # 读取分组
    df_temp = pd.read_csv(sc_template, sep="\t")
    print(f"\n📄 样本分组: {df_temp[group_col].value_counts().to_dict()}")

    # 读取h5ad
    print(f"\n📂 读取h5ad: {adata_path}")
    adata = sc.read_h5ad(adata_path)
    print(f"✅ 数据: {adata.n_obs} 细胞 | {adata.n_vars} 基因")
    print(f"✅ 细胞类型数: {len(adata.obs[celltype_col].unique())}")

    # 加载复方靶点
    formula2target = load_formula_targets(formula_pkl)
    if len(formula2target) == 0:
        print("❌ 错误：未加载到任何有效复方！")
        return
    print(f"• 平均每复方靶点: {np.mean([len(v) for v in formula2target.values()]):.2f}")

    Path(output_dir).mkdir(exist_ok=True, parents=True)

    # ====================== drug2cell 打分 ======================
    print("\n🧪 运行 drug2cell.score ...")
    d2c.score(adata, targets=formula2target, use_raw=True)
    d2c_adata = adata.uns["drug2cell"]

    # ====================== 差异分析 CASE vs CTRL ======================
    print(f"\n📈 差异分析: {case_name} vs {ctrl_name}")
    sc.tl.rank_genes_groups(
        d2c_adata,
        groupby=group_col,
        reference=ctrl_name,
        method="wilcoxon"
    )

    # 解析结果
    rank_res = d2c_adata.uns['rank_genes_groups']
    all_formulas = list(formula2target.keys())
    diff_metrics = {}

    for formula in all_formulas:
        if formula in rank_res['names'][case_name]:
            idx = list(rank_res['names'][case_name]).index(formula)
            diff_metrics[formula] = {
                'logfoldchange': rank_res['logfoldchanges'][case_name][idx],
                'pval': rank_res['pvals'][case_name][idx],
                'pval_adj': rank_res['pvals_adj'][case_name][idx],
                'score': rank_res['scores'][case_name][idx]
            }
        else:
            diff_metrics[formula] = {
                'logfoldchange': np.nan,
                'pval': np.nan,
                'pval_adj': np.nan,
                'score': np.nan
            }

    # logFC 过滤
    satisfy_lfc = []
    for formula, m in diff_metrics.items():
        lfc = m['logfoldchange']
        if not np.isnan(lfc) and abs(lfc) >= logfc_threshold:
            satisfy_lfc.append((formula, m))

    print("\n📊 logFC过滤结果")
    print(f"• 总复方: {len(all_formulas)}")
    print(f"• 满足|logFC|≥{logfc_threshold}: {len(satisfy_lfc)}")

    # ====================== 双条件过滤：细胞类型内表达占比 ======================
    pass_both = []
    cell_types = d2c_adata.obs[celltype_col].unique()

    for formula, _ in satisfy_lfc:
        if formula in d2c_adata.var_names:
            max_pct = 0.0
            for ct in cell_types:
                ct_mask = d2c_adata.obs[celltype_col] == ct
                adata_ct = d2c_adata[ct_mask, formula]
                x = adata_ct.X
                total = x.shape[0]
                if total == 0:
                    continue
                nonzero = x.nnz if issparse(x) else np.count_nonzero(x)
                pct = nonzero / total
                if pct > max_pct:
                    max_pct = pct

            if max_pct >= min_pct:
                pass_both.append(formula)

    # 绘图
    plot_list = pass_both.copy() if use_filtered_only else all_formulas.copy()
    print(f"\n🔍 严格模式：仅绘制双条件满足复方（{len(plot_list)}个）")

    # 输出表格
    df_diff = pd.DataFrame.from_dict(diff_metrics, orient='index')
    df_diff.reset_index(names='formula', inplace=True)
    df_diff['pass_logFC'] = df_diff['formula'].isin([f for f, _ in satisfy_lfc])
    df_diff['pass_both'] = df_diff['formula'].isin(pass_both)
    df_diff.to_csv(Path(output_dir)/"d2c_formulaScan-all_diff.csv", index=False)

    if pass_both:
        df_pass = df_diff[df_diff['pass_both']]
        df_pass.to_csv(Path(output_dir)/"d2c_formulaScan-pass_both.csv", index=False)

    # 绘图
    print("\n📊 生成图表")
    if len(plot_list) > 0:
        # 差异点图
        sc.pl.rank_genes_groups_dotplot(
            d2c_adata, groupby=group_col, var_names=plot_list[:10],
            swap_axes=True, dendrogram=False, show=False
        )
        plt.savefig(Path(output_dir)/"d2c_formulaScan-diff_dotplot.pdf", bbox_inches="tight")
        plt.close()

        # 细胞类型点图
        sc.pl.dotplot(
            d2c_adata, var_names=plot_list[:10], groupby=celltype_col,
            swap_axes=True, show=False
        )
        plt.savefig(Path(output_dir)/"d2c_formulaScan-celltype_dotplot.pdf", bbox_inches="tight")
        plt.close()

        # UMAP
        top10 = plot_list[:10]
        n = len(top10)
        if n > 0:
            n_rows = (n + 3) // 4
            n_cols = min(4, n)
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows), dpi=300)
            axes = axes.flatten() if n > 1 else [axes]
            for i, formula in enumerate(top10):
                sc.pl.umap(d2c_adata, color=formula, cmap="OrRd", ax=axes[i], show=False, title=formula)
            for j in range(n, len(axes)):
                axes[j].set_visible(False)
            plt.tight_layout()
            plt.savefig(Path(output_dir)/"d2c_formulaScan-umap_top10.pdf", bbox_inches="tight")
            plt.close()

        print("✅ 所有图表已生成")

    print("\n📈 最终统计")
    print(f"• 总复方: {len(formula2target)}")
    print(f"• 满足logFC: {len(satisfy_lfc)}")
    print(f"• 双条件通过: {len(pass_both)}")
    print(f"• 绘图复方: {len(plot_list)}")
    print(f"\n🎉 完成！结果: {output_dir}")


# ====================== 测试 ======================
if __name__ == "__main__":
    SC_TEMPLATE = "sc_io-input_template.txt"
    ADATA_PATH = "sc_recode4anno-celltype_manual_annotated.h5ad"
    FORMULA_PKL = "tcm_form2targetall.pkl"
    OUTPUT_DIR = "./formulaScan_result"

    formulaScan(
        sc_template=SC_TEMPLATE,
        adata_path=ADATA_PATH,
        formula_pkl=FORMULA_PKL,
        output_dir=OUTPUT_DIR,
        min_pct=0.01,
        logfc_threshold=2.0,
        use_filtered_only=True
    )
