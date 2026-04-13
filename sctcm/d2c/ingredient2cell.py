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

def load_ingredient_targets(pkl_path):
    with open(pkl_path, 'rb') as f:
        cid2target = pickle.load(f)
    cid2target = {str(k).strip(): v for k, v in cid2target.items()}
    cid2target = {cid: ts for cid, ts in cid2target.items() if isinstance(ts, list) and len(ts) > 0}
    print(f"✅ 加载成分靶点：{len(cid2target)} 个有效成分")
    return cid2target

def run_ingredient2cell(
    sc_template,
    adata_path,
    ingredient_pkl,
    output_dir,
    min_pct=0.25,
    logfc_threshold=2.0,
    use_filtered_only=False
):
    print("\n" + "=" * 80)
    print("🚀 运行 tcm_ingredient2cell")
    print("=" * 80)

    # ========== 打印运行参数 ==========
    print("\n📋 当前运行参数：")
    print(f"   - 分组文件路径: {sc_template}")
    print(f"   - 单细胞数据路径: {adata_path}")
    print(f"   - 成分靶点文件路径: {ingredient_pkl}")
    print(f"   - 输出目录: {output_dir}")
    print(f"   - 最小细胞表达占比 (min_pct): {min_pct}")
    print(f"   - logFC阈值 (logfc_threshold): {logfc_threshold}")
    print(f"   - 启用双条件过滤绘图 (use_filtered_only): {use_filtered_only}")

    # 读取分组
    df_temp = pd.read_csv(sc_template, sep="\t")
    print(f"\n📄 分组信息：\n{df_temp['Group'].value_counts().to_dict()}")
    print(f"   - CTRL={df_temp['Group'].value_counts().get('CTRL',0)} | CASE={df_temp['Group'].value_counts().get('CASE',0)}")

    # 读取单细胞数据
    print(f"\n📂 读取h5ad：{adata_path}")
    adata = sc.read_h5ad(adata_path)
    print(f"✅ 数据维度：{adata.n_obs} 细胞 | {adata.n_vars} 基因")
    print(f"   - 细胞类型数量: {len(adata.obs['celltype_manual'].unique())}")

    # 加载成分靶点
    targets = load_ingredient_targets(ingredient_pkl)
    print(f"   - 平均每个成分靶点数: {np.mean([len(v) for v in targets.values()]):.2f}")
    Path(output_dir).mkdir(exist_ok=True, parents=True)

    celltype_col = "celltype_manual"
    group_col = "Group"

    # drug2cell 打分
    print("\n🧪 运行 d2c.score ...")
    d2c.score(adata, targets=targets, use_raw=True)
    d2c_adata = adata.uns["drug2cell"]

    # 差异分析 CASE vs CTRL
    print("\n📈 差异分析：CASE vs CTRL")
    sc.tl.rank_genes_groups(d2c_adata, groupby=group_col, reference="CTRL", method="wilcoxon")

    # 解析差异结果
    rank_results = d2c_adata.uns['rank_genes_groups']
    ingredient_names = list(targets.keys())
    diff_metrics = {}
    
    for ing in ingredient_names:
        if ing in rank_results['names']['CASE']:
            idx = list(rank_results['names']['CASE']).index(ing)
            diff_metrics[ing] = {
                'logfoldchange': rank_results['logfoldchanges']['CASE'][idx],
                'pval': rank_results['pvals']['CASE'][idx],
                'pval_adj': rank_results['pvals_adj']['CASE'][idx],
                'score': rank_results['scores']['CASE'][idx]
            }
        else:
            diff_metrics[ing] = {
                'logfoldchange': np.nan,
                'pval': np.nan,
                'pval_adj': np.nan,
                'score': np.nan
            }
    
    # 筛选满足 logFC 的成分
    satisfy_ingredients = []
    unsatisfy_ingredients = []
    for ing, m in diff_metrics.items():
        lfc = m['logfoldchange']
        if not np.isnan(lfc) and abs(lfc) >= logfc_threshold:
            satisfy_ingredients.append((ing, m))
        else:
            unsatisfy_ingredients.append((ing, m))
    
    # 打印统计
    print("\n📊 差异统计（CASE vs CTRL）：")
    print(f"   - 总成分数: {len(ingredient_names)}")
    print(f"   - 满足logFC: {len(satisfy_ingredients)} ({len(satisfy_ingredients)/len(ingredient_names)*100:.1f}%)")
    print(f"   - 不满足logFC: {len(unsatisfy_ingredients)} ({len(unsatisfy_ingredients)/len(ingredient_names)*100:.1f}%)")

    if satisfy_ingredients:
        print("\n✅ 满足 logFC 阈值的成分：")
        for ing, m in satisfy_ingredients:
            print(f"   {ing}: logFC = {m['logfoldchange']:.4f}")

    if unsatisfy_ingredients:
        print("\n❌ 不满足 logFC 阈值的成分：")
        for ing, m in unsatisfy_ingredients:
            lfc_str = f"{m['logfoldchange']:.4f}" if not np.isnan(m['logfoldchange']) else "N/A"
            print(f"   {ing}: logFC = {lfc_str}")

    # ====================== 核心：双条件过滤（logFC + 细胞表达占比）======================

    true_both_satisfy = []
    for ing, _ in satisfy_ingredients:
        if ing in d2c_adata.var_names:
            expr = d2c_adata[:, ing].X
            n_total = expr.shape[0]
            n_nonzero = expr.nnz if issparse(expr) else np.count_nonzero(expr)
            pct = n_nonzero / n_total
            if pct >= min_pct:
                true_both_satisfy.append(ing)

    # 用于绘图的成分列表
    final_plot_ingredients = []
    if use_filtered_only:
        print(f"\n🔍 严格过滤模式：仅绘制【logFC≥{logfc_threshold}】+【细胞表达占比≥{min_pct}】的成分")
        final_plot_ingredients = true_both_satisfy.copy()
    else:
        print(f"\n📌 不过滤模式：绘制所有成分")
        final_plot_ingredients = ingredient_names

    n_final = len(final_plot_ingredients)
    print(f"\n🎯 最终绘图成分数量：{n_final} 个")
    if n_final == 0:
        print("⚠️  无符合条件的成分 → 跳过所有绘图")

    # ====================== 绘图 ======================
    print("\n📊 开始绘图...")
    import matplotlib.pyplot as plt

    # 1. 差异点图
    if final_plot_ingredients:
        sc.pl.rank_genes_groups_dotplot(
            d2c_adata,
            groupby=group_col,
            swap_axes=True,
            dendrogram=False,
            n_genes=min(5, len(final_plot_ingredients)),
            show=False
        )
        plt.savefig(f"{output_dir}/d2c_ingredient2cell-dotplot_top.pdf", bbox_inches="tight")
        plt.close()
        print("   - ✅ 差异点图已保存")
    else:
        print("   - ⚠️  跳过差异点图")

    # 2. 细胞类型点图
    if final_plot_ingredients:
        plot_list = final_plot_ingredients[:10]
        sc.pl.dotplot(
            d2c_adata,
            var_names=plot_list,
            groupby=celltype_col,
            swap_axes=True,
            show=False
        )
        plt.savefig(f"{output_dir}/d2c_ingredient2cell-dotplot_celltypes.pdf", bbox_inches="tight")
        plt.close()
        print(f"   - ✅ 细胞类型点图已保存（{len(plot_list)} 个成分）")
    else:
        print("   - ⚠️  跳过细胞类型点图")

    # 3. UMAP 组图 TOP10 + OrRd
    if final_plot_ingredients:
        plot_list = final_plot_ingredients[:10]
        n = len(plot_list)
        print(f"   - 生成 UMAP 组图（TOP{n}，配色 OrRd）")

        n_rows = (n + 3) // 4
        n_cols = min(4, n)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows), dpi=300)
        axes = axes.flatten() if n > 1 else [axes]

        for i, ing in enumerate(plot_list):
            sc.pl.umap(d2c_adata, color=ing, cmap="OrRd", ax=axes[i], show=False, title=ing)

        for i in range(n, len(axes)):
            axes[i].set_visible(False)

        plt.tight_layout()
        plt.savefig(f"{output_dir}/d2c_ingredient2cell-umap_top10.pdf", bbox_inches="tight")
        plt.close()
        print("   - ✅ UMAP 组图已保存")
    else:
        print("   - ⚠️  跳过 UMAP 组图")

    # 最终汇总（逻辑完全正确）
    print("\n📈 最终真实统计：")
    print(f"   - 总成分: {len(targets)}")
    print(f"   - 满足logFC: {len(satisfy_ingredients)}")
    print(f"   - 同时满足logFC+细胞占比: {len(true_both_satisfy)}")
    print(f"   - 本次实际绘图成分: {len(final_plot_ingredients)}")

    print(f"\n🎉 全部完成！结果：{output_dir}")
