import scanpy as sc
import drug2cell as d2c
import pandas as pd
import pickle
import warnings
from pathlib import Path
import numpy as np
from scipy.sparse import issparse

warnings.filterwarnings("ignore")
sc.settings.set_figure_params(dpi=300, fontsize=10, facecolor="white", figsize=(6,4))

def load_formula_as_whole(pkl_path):
    with open(pkl_path, 'rb') as f:
        formula2herb2cid2target = pickle.load(f)

    formula2target = {}
    for formula_name, herb_dict in formula2herb2cid2target.items():
        all_targets = []
        for herb_name, cid_dict in herb_dict.items():
            for cid, targets in cid_dict.items():
                if isinstance(targets, list) and len(targets) > 0:
                    all_targets.extend(targets)

        unique_targets = sorted(list(set([t.strip() for t in all_targets if t.strip()])))
        if len(unique_targets) > 0:
            formula2target[formula_name] = unique_targets
        else:
            print(f"⚠️ 无有效靶点，跳过：{formula_name}")

    print(f"✅ 配方整体靶点处理完成：")
    print(f"   - 有效配方数：{len(formula2target)} 个")
    for formula, tgts in formula2target.items():
        print(f"   → {formula}：{len(tgts)} 个去重靶点")
    return formula2target

def run_formula2cell(
    sc_template,
    adata_path,
    formula_pkl,
    output_dir,
    min_pct=0.25,
    logfc_threshold=2.0,
    use_filtered_only=False
):
    print("\n" + "=" * 80)
    print("🚀 运行 tcm_formula2cell（配方整体 → 细胞分析）")
    print("📌 模式：Formula-Target（配方整体，靶点去重）")
    print("=" * 80)

    print("\n📋 当前运行参数：")
    print(f"   - 分组模板: {sc_template}")
    print(f"   - 单细胞数据: {adata_path}")
    print(f"   - 配方靶点文件: {formula_pkl}")
    print(f"   - 输出目录: {output_dir}")
    print(f"   - 最小细胞表达占比 min_pct: {min_pct}")
    print(f"   - logFC 阈值: {logfc_threshold}")
    print(f"   - 仅绘制满足双条件的配方: {use_filtered_only}")

    df_temp = pd.read_csv(sc_template, sep="\t")
    print(f"\n📄 分组信息：\n{df_temp['Group'].value_counts().to_dict()}")
    ctrl_num = df_temp['Group'].value_counts().get('CTRL', 0)
    case_num = df_temp['Group'].value_counts().get('CASE', 0)
    print(f"   - CTRL={ctrl_num} | CASE={case_num}")

    print(f"\n📂 读取h5ad：{adata_path}")
    adata = sc.read_h5ad(adata_path)
    print(f"✅ 数据：{adata.n_obs} 细胞 | {adata.n_vars} 基因")
    celltype_col = "celltype_manual"
    group_col = "Group"
    print(f"   - 细胞类型数量：{len(adata.obs[celltype_col].unique())}")

    formula2target = load_formula_as_whole(formula_pkl)
    Path(output_dir).mkdir(exist_ok=True, parents=True)

    if not formula2target:
        print("❌ 无有效配方，退出")
        return

    print("\n🧪 运行 d2c.score（配方整体打分）...")
    d2c.score(adata, targets=formula2target, use_raw=True)
    d2c_adata = adata.uns["drug2cell"]

    print("\n📈 差异分析：CASE vs CTRL（配方整体）")
    sc.tl.rank_genes_groups(
        d2c_adata,
        groupby=group_col,
        reference="CTRL",
        method="wilcoxon"
    )

    rank_results = d2c_adata.uns['rank_genes_groups']
    formula_names = list(formula2target.keys())
    diff_metrics = {}

    for formula in formula_names:
        if formula in rank_results['names']['CASE']:
            idx = list(rank_results['names']['CASE']).index(formula)
            diff_metrics[formula] = {
                'logfoldchange': rank_results['logfoldchanges']['CASE'][idx],
                'pval': rank_results['pvals']['CASE'][idx],
                'pval_adj': rank_results['pvals_adj']['CASE'][idx],
                'score': rank_results['scores']['CASE'][idx],
            }
        else:
            diff_metrics[formula] = {
                'logfoldchange': np.nan,
                'pval': np.nan,
                'pval_adj': np.nan,
                'score': np.nan,
            }

    satisfy_logfc = []
    unsatisfy_logfc = []
    for formula, m in diff_metrics.items():
        lfc = m['logfoldchange']
        if not np.isnan(lfc) and abs(lfc) >= logfc_threshold:
            satisfy_logfc.append(formula)
        else:
            unsatisfy_logfc.append(formula)

    print("\n📊 logFC 统计：")
    print(f"   - 满足 logFC ≥ {logfc_threshold}：{len(satisfy_logfc)} 个")
    for f in satisfy_logfc:
        print(f"     ✅ {f} | logFC={diff_metrics[f]['logfoldchange']:.4f}")
    print(f"   - 不满足 logFC：{len(unsatisfy_logfc)} 个")
    for f in unsatisfy_logfc:
        lfc = diff_metrics[f]['logfoldchange']
        lfc_str = f"{lfc:.4f}" if not np.isnan(lfc) else "N/A"
        print(f"     ❌ {f} | logFC={lfc_str}")

    print("\n📊 细胞表达占比统计（cell percentage）：")
    cell_pct_stats = {}
    satisfy_cell_pct = []
    unsatisfy_cell_pct = []

    for formula in formula_names:
        if formula in d2c_adata.var_names:
            expr = d2c_adata[:, formula].X
            n_total = expr.shape[0]
            n_nonzero = expr.nnz if issparse(expr) else np.count_nonzero(expr)
            pct = n_nonzero / n_total
            cell_pct_stats[formula] = pct

            if pct >= min_pct:
                satisfy_cell_pct.append(formula)
                mark = "✅"
            else:
                unsatisfy_cell_pct.append(formula)
                mark = "❌"
            print(f"   {mark} {formula} | 细胞占比 = {pct:.1%}")
        else:
            cell_pct_stats[formula] = 0.0
            print(f"   ❌ {formula} | 细胞占比 = 0%")

    print(f"\n📌 细胞占比 ≥ {min_pct:.0%} 的配方：{len(satisfy_cell_pct)} 个")
    print(f"📌 细胞占比 ＜ {min_pct:.0%} 的配方：{len(unsatisfy_cell_pct)} 个")

    both_ok = [f for f in satisfy_logfc if f in satisfy_cell_pct]
    print("\n📊 双条件同时满足（logFC + 细胞占比）：")
    print(f"   - 符合配方数：{len(both_ok)}")
    for f in both_ok:
        print(f"     ✅ {f}")

    if use_filtered_only:
        print(f"\n🎯 use_filtered_only = {use_filtered_only} → 仅绘制双条件满足的配方")
        final_plot = both_ok.copy()
    else:
        print(f"\n🎯 use_filtered_only = {use_filtered_only} → 不过滤，绘制全部配方")
        final_plot = formula_names.copy()

    print(f"\n🎯 最终绘图配方数：{len(final_plot)}")

    print("\n📊 开始绘图...")
    import matplotlib.pyplot as plt

    if final_plot:
        sc.pl.rank_genes_groups_dotplot(
            d2c_adata, groupby=group_col, swap_axes=True, dendrogram=False,
            n_genes=min(5, len(final_plot)), show=False
        )
        plt.savefig(f"{output_dir}/d2c_formula2cell-dotplot_top.pdf", bbox_inches="tight")
        plt.close()

    if final_plot:
        plot_list = final_plot[:10]
        sc.pl.dotplot(
            d2c_adata, var_names=plot_list, groupby=celltype_col,
            swap_axes=True, show=False
        )
        plt.savefig(f"{output_dir}/d2c_formula2cell-dotplot_celltypes.pdf", bbox_inches="tight")
        plt.close()

    if final_plot:
        plot_list = final_plot[:10]
        n = len(plot_list)
        n_rows = (n + 3) // 4
        n_cols = min(4, n)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows), dpi=300)
        axes = axes.flatten() if n > 1 else [axes]
        for i, formula in enumerate(plot_list):
            sc.pl.umap(d2c_adata, color=formula, cmap="OrRd", ax=axes[i], show=False, title=formula)
        for i in range(n, len(axes)):
            axes[i].set_visible(False)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/d2c_formula2cell-umap_top10.pdf", bbox_inches="tight")
        plt.close()

    print("\n📈 最终汇总：")
    print(f"   - 总配方：{len(formula2target)}")
    print(f"   - 满足 logFC：{len(satisfy_logfc)}")
    print(f"   - 满足细胞占比：{len(satisfy_cell_pct)}")
    print(f"   - 双条件都满足：{len(both_ok)}")
    print(f"   - 最终绘图：{len(final_plot)}")
    print(f"\n🎉 formula2cell 运行完成！")

if __name__ == "__main__":
    run_formula2cell(
        sc_template="/scTCM/bin/v0.1.0/docs/templates/sc_io-input_template.txt",
        adata_path="/scTCM/bin/v0.1.0/sctcm/tests/result/sc_recode4anno-celltype_manual_annotated.h5ad",
        formula_pkl="/scTCM/bin/v0.1.0/sctcm/tests/result/tcm_formula2target-formula2target.pkl",
        output_dir="/scTCM/bin/v0.1.0/sctcm/tests/result",
        min_pct=0.25,
        logfc_threshold=2.0,
        use_filtered_only=False
    )
