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

def load_herb_as_whole(pkl_path):
    """
    【核心】加载草药 → 把所有成分靶点合并、去重 → 输出 {草药名: [去重靶点列表]}
    完全按 Herb 整体分析
    """
    with open(pkl_path, 'rb') as f:
        herb2cid2target = pickle.load(f)

    herb2target = {}
    total_herbs = 0

    for herb_name, cid_dict in herb2cid2target.items():
        # 收集该草药下所有成分的所有靶点
        all_targets = []
        for cid, targets in cid_dict.items():
            if isinstance(targets, list) and len(targets) > 0:
                all_targets.extend(targets)

        # 去重 + 去空
        unique_targets = sorted(list(set([t.strip() for t in all_targets if t.strip()])))

        if len(unique_targets) > 0:
            herb2target[herb_name] = unique_targets
            total_herbs += 1

    print(f"✅ 草药整体靶点处理完成（Herb-Target 模式）：")
    print(f"   - 有效草药数：{len(herb2target)} 个")
    print(f"   - 全部按【草药整体】分析，靶点已去重")

    for herb, tgts in herb2target.items():
        print(f"   → {herb}：{len(tgts)} 个去重靶点")

    return herb2target

def run_herb2cell(
    sc_template,
    adata_path,
    herb_pkl,
    output_dir,
    min_pct=0.25,
    logfc_threshold=2.0,
    use_filtered_only=False
):
    print("\n" + "=" * 80)
    print("🚀 运行 tcm_herb2cell（Herb 整体 → 细胞分析）")
    print("📌 模式：Herb-Target（草药整体，靶点去重）")
    print("=" * 80)

    # ========== 打印参数 ==========
    print("\n📋 当前运行参数：")
    print(f"   - 分组模板: {sc_template}")
    print(f"   - 单细胞数据: {adata_path}")
    print(f"   - 草药靶点文件: {herb_pkl}")
    print(f"   - 输出目录: {output_dir}")
    print(f"   - 最小细胞表达占比 min_pct: {min_pct}")
    print(f"   - logFC 阈值: {logfc_threshold}")
    print(f"   - 仅绘制满足双条件的成分: {use_filtered_only}")

    # 读取分组
    df_temp = pd.read_csv(sc_template, sep="\t")
    print(f"\n📄 分组信息：\n{df_temp['Group'].value_counts().to_dict()}")
    ctrl_num = df_temp['Group'].value_counts().get('CTRL', 0)
    case_num = df_temp['Group'].value_counts().get('CASE', 0)
    print(f"   - CTRL={ctrl_num} | CASE={case_num}")

    # 读取单细胞数据
    print(f"\n📂 读取h5ad：{adata_path}")
    adata = sc.read_h5ad(adata_path)
    print(f"✅ 数据：{adata.n_obs} 细胞 | {adata.n_vars} 基因")
    celltype_col = "major_celltype"
    group_col = "group"
    print(f"   - 细胞类型数量：{len(adata.obs[celltype_col].unique())}")

    # ========== 核心：加载草药整体靶点 ==========
    herb2target = load_herb_as_whole(herb_pkl)
    Path(output_dir).mkdir(exist_ok=True, parents=True)

    if not herb2target:
        print("❌ 无有效草药，退出")
        return

    # ========== drug2cell 打分（Herb 整体）==========
    print("\n🧪 运行 d2c.score（草药整体打分）...")
    d2c.score(adata, targets=herb2target, use_raw=True)
    d2c_adata = adata.uns["drug2cell"]

    # ========== 差异分析 CASE vs CTRL ==========
    print("\n📈 差异分析：CASE vs CTRL（草药整体）")
    sc.tl.rank_genes_groups(
        d2c_adata,
        groupby=group_col,
        reference="CTRL",
        method="wilcoxon"
    )

    # ========== 解析差异结果 ==========
    rank_results = d2c_adata.uns['rank_genes_groups']
    herb_names = list(herb2target.keys())
    diff_metrics = {}

    for herb in herb_names:
        if herb in rank_results['names']['CASE']:
            idx = list(rank_results['names']['CASE']).index(herb)
            diff_metrics[herb] = {
                'logfoldchange': rank_results['logfoldchanges']['CASE'][idx],
                'pval': rank_results['pvals']['CASE'][idx],
                'pval_adj': rank_results['pvals_adj']['CASE'][idx],
                'score': rank_results['scores']['CASE'][idx],
            }
        else:
            diff_metrics[herb] = {
                'logfoldchange': np.nan,
                'pval': np.nan,
                'pval_adj': np.nan,
                'score': np.nan,
            }

    # ========== 筛选满足 logFC 的草药 ==========
    satisfy_herbs = []
    unsatisfy_herbs = []
    for herb, m in diff_metrics.items():
        lfc = m['logfoldchange']
        if not np.isnan(lfc) and abs(lfc) >= logfc_threshold:
            satisfy_herbs.append((herb, m))
        else:
            unsatisfy_herbs.append((herb, m))

    # ========== 打印统计 ==========
    print("\n📊 草药差异统计（CASE vs CTRL）：")
    total = len(herb_names)
    satisfy_cnt = len(satisfy_herbs)
    print(f"   - 总草药数：{total}")
    print(f"   - 满足 logFC≥{logfc_threshold}：{satisfy_cnt} ({satisfy_cnt/total*100:.1f}%)")
    print(f"   - 不满足：{len(unsatisfy_herbs)} ({len(unsatisfy_herbs)/total*100:.1f}%)")

    if satisfy_herbs:
        print("\n✅ 满足 logFC 的草药：")
        for herb, m in satisfy_herbs:
            print(f"   {herb:<12} | logFC={m['logfoldchange']:.4f}")

    if unsatisfy_herbs:
        print("\n❌ 不满足 logFC 的草药：")
        for herb, m in unsatisfy_herbs:
            lfc_str = f"{m['logfoldchange']:.4f}" if not np.isnan(m['logfoldchange']) else "N/A"
            print(f"   {herb:<12} | logFC={lfc_str}")

    # ========== 双条件过滤：logFC + 细胞表达占比 ==========
    true_both_satisfy = []
    for herb, _ in satisfy_herbs:
        if herb in d2c_adata.var_names:
            expr = d2c_adata[:, herb].X
            n_total = expr.shape[0]
            n_nonzero = expr.nnz if issparse(expr) else np.count_nonzero(expr)
            pct = n_nonzero / n_total
            if pct >= min_pct:
                true_both_satisfy.append(herb)

    # 绘图列表
    if use_filtered_only:
        print(f"\n🔍 严格过滤：仅保留 logFC≥{logfc_threshold} + 细胞占比≥{min_pct}")
        final_plot = true_both_satisfy.copy()
    else:
        print(f"\n📌 不过滤：绘制全部草药")
        final_plot = herb_names

    n_final = len(final_plot)
    print(f"\n🎯 最终绘图草药数：{n_final}")
    if n_final == 0:
        print("⚠️ 无符合条件草药，跳过绘图")

    # ========== 绘图 ==========
    print("\n📊 开始绘图（草药整体）...")
    import matplotlib.pyplot as plt

    # 1 差异点图
    if final_plot:
        sc.pl.rank_genes_groups_dotplot(
            d2c_adata, groupby=group_col, swap_axes=True, dendrogram=False,
            n_genes=min(5, len(final_plot)), show=False
        )
        plt.savefig(f"{output_dir}/d2c_herb2cell-dotplot_top.pdf", bbox_inches="tight")
        plt.close()
        print("   - ✅ 草药差异点图已保存")

    # 2 细胞类型点图
    if final_plot:
        plot_list = final_plot[:10]
        sc.pl.dotplot(
            d2c_adata, var_names=plot_list, groupby=celltype_col,
            swap_axes=True, show=False
        )
        plt.savefig(f"{output_dir}/d2c_herb2cell-dotplot_celltypes.pdf", bbox_inches="tight")
        plt.close()
        print(f"   - ✅ 细胞类型点图（{len(plot_list)} 个草药）")

    # 3 UMAP TOP10 + OrRd
    if final_plot:
        plot_list = final_plot[:10]
        n = len(plot_list)
        print(f"   - 生成草药 UMAP 组图 TOP{n}（配色 OrRd）")

        n_rows = (n + 3) // 4
        n_cols = min(4, n)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows), dpi=300)
        axes = axes.flatten() if n > 1 else [axes]

        for i, herb in enumerate(plot_list):
            sc.pl.umap(d2c_adata, color=herb, cmap="OrRd", ax=axes[i], show=False, title=herb)

        for i in range(n, len(axes)):
            axes[i].set_visible(False)

        plt.tight_layout()
        plt.savefig(f"{output_dir}/d2c_herb2cell-umap_top10.pdf", bbox_inches="tight")
        plt.close()
        print("   - ✅ 草药 UMAP 组图已保存")

    # ========== 最终汇总 ==========
    print("\n📈 最终统计（Herb 整体）：")
    print(f"   - 总草药：{len(herb2target)}")
    print(f"   - 满足 logFC：{len(satisfy_herbs)}")
    print(f"   - 满足双条件：{len(true_both_satisfy)}")
    print(f"   - 实际绘图：{len(final_plot)}")

    print(f"\n🎉 herb2cell（草药整体分析）完成！")
