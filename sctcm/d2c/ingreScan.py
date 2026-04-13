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


def load_ingredient_targets(pkl_path):
    """加载成分-靶点字典（tcm_ing2targetall.pkl）"""
    with open(pkl_path, 'rb') as f:
        cid2target = pickle.load(f)
    cid2target = {str(k).strip(): v for k, v in cid2target.items()}
    cid2target = {cid: ts for cid, ts in cid2target.items()
                  if isinstance(ts, list) and len(ts) > 0}
    print(f"✅ 加载全量成分靶点：{len(cid2target)} 个有效成分")
    return cid2target


def ingreScan(
    sc_template,
    adata_path,
    ingredient_pkl,
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
    全成分批量drug2cell分析（CASE vs CTRL）
    筛选标准：
    1. logFC > 阈值 (vs 其他细胞)
    2. 靶基因在【感兴趣细胞类型内部】至少 min_pct 细胞表达

    Parameters
    ----------
    sc_template: str
        样本分组文件路径 sc_io-input_template.txt
    adata_path: str
        已注释单细胞h5ad路径
    ingredient_pkl: str
        全成分靶点文件 tcm_ing2targetall.pkl
    output_dir: str
        输出目录
    min_pct: float
        最小细胞表达占比（细胞类型内部）
    logfc_threshold: float
        logFC绝对值阈值
    use_filtered_only: bool
        True=仅绘制双条件满足成分；False=绘制所有成分
    celltype_col: str
        细胞类型列名
    group_col: str
        分组列名
    case_name: str
        病例组名称
    ctrl_name: str
        对照组名称
    """
    print("\n" + "=" * 80)
    print("🚀 ingreScan - 中药全成分drug2cell批量扫描")
    print("=" * 80)

    # 1. 打印运行参数
    print("\n📋 运行参数")
    print(f"• 分组文件: {sc_template}")
    print(f"• 单细胞数据: {adata_path}")
    print(f"• 全成分靶点: {ingredient_pkl}")
    print(f"• 输出目录: {output_dir}")
    print(f"• min_pct(细胞类型内): {min_pct} | logfc_threshold: {logfc_threshold}")
    print(f"• 严格过滤绘图: {use_filtered_only}")
    print(f"• 分组: {case_name} vs {ctrl_name}")

    # 2. 读取样本分组
    df_temp = pd.read_csv(sc_template, sep="\t")
    print(f"\n📄 样本分组: {df_temp[group_col].value_counts().to_dict()}")

    # 3. 读取注释后单细胞数据
    print(f"\n📂 读取h5ad: {adata_path}")
    adata = sc.read_h5ad(adata_path)
    print(f"✅ 数据: {adata.n_obs} 细胞 | {adata.n_vars} 基因")
    print(f"✅ 细胞类型数: {len(adata.obs[celltype_col].unique())}")

    # 4. 加载全成分靶点
    targets = load_ingredient_targets(ingredient_pkl)
    print(f"• 平均每成分靶点: {np.mean([len(v) for v in targets.values()]):.2f}")

    # 创建输出目录
    Path(output_dir).mkdir(exist_ok=True, parents=True)

    # ---------------------- 核心drug2cell流程 ----------------------
    print("\n🧪 运行 drug2cell.score ...")
    d2c.score(adata, targets=targets, use_raw=True)
    d2c_adata = adata.uns["drug2cell"]

    # 差异分析 CASE vs CTRL
    print(f"\n📈 差异分析: {case_name} vs {ctrl_name}")
    sc.tl.rank_genes_groups(
        d2c_adata,
        groupby=group_col,
        reference=ctrl_name,
        method="wilcoxon"
    )

    # 解析所有成分差异结果
    rank_res = d2c_adata.uns['rank_genes_groups']
    all_ingredients = list(targets.keys())
    diff_metrics = {}

    for ing in all_ingredients:
        if ing in rank_res['names'][case_name]:
            idx = list(rank_res['names'][case_name]).index(ing)
            diff_metrics[ing] = {
                'logfoldchange': rank_res['logfoldchanges'][case_name][idx],
                'pval': rank_res['pvals'][case_name][idx],
                'pval_adj': rank_res['pvals_adj'][case_name][idx],
                'score': rank_res['scores'][case_name][idx]
            }
        else:
            diff_metrics[ing] = {
                'logfoldchange': np.nan,
                'pval': np.nan,
                'pval_adj': np.nan,
                'score': np.nan
            }

    # ---------------------- logFC过滤 ----------------------
    satisfy_lfc = []
    unsatisfy_lfc = []
    for ing, m in diff_metrics.items():
        lfc = m['logfoldchange']
        if not np.isnan(lfc) and abs(lfc) >= logfc_threshold:
            satisfy_lfc.append((ing, m))
        else:
            unsatisfy_lfc.append((ing, m))

    print("\n📊 logFC过滤结果")
    print(f"• 总成分: {len(all_ingredients)}")
    print(f"• 满足|logFC|≥{logfc_threshold}: {len(satisfy_lfc)}")
    print(f"• 不满足: {len(unsatisfy_lfc)}")

    # ====================== 双条件过滤：logFC + 细胞类型内表达占比 ======================
    # 文献：expressed in at least X% of cells within the cell type
    pass_both = []
    cell_types = d2c_adata.obs[celltype_col].unique()

    for ing, _ in satisfy_lfc:
        if ing in d2c_adata.var_names:
            max_pct = 0.0
            # 遍历每一种细胞类型，计算该成分在【这个细胞类型内部】的表达占比
            for ct in cell_types:
                ct_mask = d2c_adata.obs[celltype_col] == ct
                adata_ct = d2c_adata[ct_mask, ing]
                
                x = adata_ct.X
                total_cells = x.shape[0]
                if total_cells == 0:
                    continue
                nonzero_cells = x.nnz if issparse(x) else np.count_nonzero(x)
                pct = nonzero_cells / total_cells
                if pct > max_pct:
                    max_pct = pct

            # 只要在任意一种细胞类型中表达占比 ≥ min_pct，就通过筛选
            if max_pct >= min_pct:
                pass_both.append(ing)

    # 绘图成分选择
    if use_filtered_only:
        plot_list = pass_both.copy()
        print(f"\n🔍 严格模式：仅绘制双条件满足成分（{len(plot_list)}个）")
    else:
        plot_list = all_ingredients.copy()
        print(f"\n📌 全量模式：绘制所有成分（{len(plot_list)}个）")

    # ---------------------- 结果输出 ----------------------

    df_diff = pd.DataFrame.from_dict(diff_metrics, orient='index')
    df_diff.reset_index(names='ingredient', inplace=True)
    df_diff['pass_logFC'] = df_diff['ingredient'].isin([i[0] for i in satisfy_lfc])
    df_diff['pass_both'] = df_diff['ingredient'].isin(pass_both)
    df_diff.to_csv(Path(output_dir)/"d2c_ingreScan-all_diff.csv", index=False)

    # 2. 保存双条件通过成分
    if pass_both:
        df_pass = df_diff[df_diff['pass_both']].copy()
        df_pass.to_csv(Path(output_dir)/"d2c_ingreScan-pass_both.csv", index=False)

    # ---------------------- 可视化绘图 ----------------------
    print("\n📊 生成图表")
    n_plot = len(plot_list)

    if n_plot == 0:
        print("⚠️ 无符合条件成分，跳过绘图")
    else:
        # 1. 分组差异点图
        sc.pl.rank_genes_groups_dotplot(
            d2c_adata,
            groupby=group_col,
            var_names=plot_list[:10],
            swap_axes=True,
            dendrogram=False,
            show=False
        )
        plt.savefig(Path(output_dir)/"d2c_ingreScan-diff_dotplot.pdf", bbox_inches="tight")
        plt.close()
        print("• 差异点图 ✅")

        # 2. 细胞类型点图
        sc.pl.dotplot(
            d2c_adata,
            var_names=plot_list[:10],
            groupby=celltype_col,
            swap_axes=True,
            show=False
        )
        plt.savefig(Path(output_dir)/"d2c_ingreScan-celltype_dotplot.pdf", bbox_inches="tight")
        plt.close()
        print("• 细胞类型点图 ✅")

        # 3. UMAP组图 TOP10
        top10 = plot_list[:10]
        n = len(top10)
        if n > 0:
            n_rows = (n + 3) // 4
            n_cols = min(4, n)
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows), dpi=300)
            axes = axes.flatten() if n > 1 else [axes]
            for i, ing in enumerate(top10):
                sc.pl.umap(d2c_adata, color=ing, cmap="OrRd", ax=axes[i], show=False, title=ing)
            for j in range(n, len(axes)):
                axes[j].set_visible(False)
            plt.tight_layout()
            plt.savefig(Path(output_dir)/"d2c_ingreScan-umap_top10.pdf", bbox_inches="tight")
            plt.close()
            print("• UMAP组图 ✅")

    # ---------------------- 最终统计 ----------------------
    print("\n📈 最终统计")
    print(f"• 总成分: {len(targets)}")
    print(f"• 满足logFC: {len(satisfy_lfc)}")
    print(f"• 双条件通过: {len(pass_both)}")
    print(f"• 绘图成分: {len(plot_list)}")
    print(f"\n🎉 完成！结果: {output_dir}")


# ====================== 测试入口 ======================
if __name__ == "__main__":
    # 路径配置
    SC_TEMPLATE = "sc_io-input_template.txt"
    ADATA_PATH = "sc_recode4anno-celltype_manual_annotated.h5ad"
    INGREDIENT_PKL = "tcm_ing2targetall.pkl"
    OUTPUT_DIR = "./ingreScan_result"

    # 运行
    ingreScan(
        sc_template=SC_TEMPLATE,
        adata_path=ADATA_PATH,
        ingredient_pkl=INGREDIENT_PKL,
        output_dir=OUTPUT_DIR,
        min_pct=0.01,
        logfc_threshold=2.0,
        use_filtered_only=True
    )
