# -*- coding: utf-8 -*-
import scanpy as sc
import pandas as pd
import pickle
import numpy as np
from pathlib import Path
import warnings
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sc.settings.set_figure_params(dpi=300, fontsize=10, facecolor="white")

# ====================== Jaccard 系数 ======================
def jaccard_coeff(set1, set2):
    inter = len(set1 & set2)
    union = len(set1 | set2)
    return inter / union if union > 0 else 0.0

def calculate_synergy(gene_dict, adata, synergy_scale):
    valid_sets = []
    for genes in gene_dict.values():
        s = set(g for g in genes if g in adata.var.index)
        if len(s) > 0:
            valid_sets.append(s)
    if len(valid_sets) < 2:
        return 0.0
    vals = []
    for i in range(len(valid_sets)):
        for j in range(i+1, len(valid_sets)):
            vals.append(jaccard_coeff(valid_sets[i], valid_sets[j]))
    return np.mean(vals) * synergy_scale

def get_basal_score(adata, gene_list):
    genes = [g for g in gene_list if g in adata.var.index]
    if not genes:
        return np.zeros(adata.n_obs)
    print(f"     - 基础分计算：输入基因数={len(gene_list)}, 有效基因数={len(genes)}")
    sc.tl.score_genes(adata, genes, score_name="tmp_score", copy=False)
    score = adata.obs["tmp_score"].values.copy()
    print(f"     - tmp_score统计：均值={score.mean():.4f}, 标准差={score.std():.4f}, 最小值={score.min():.4f}, 最大值={score.max():.4f}")
    del adata.obs["tmp_score"]
    return score

def load_pkl(pkl_path):
    try:
        with open(pkl_path, 'rb') as f:
            return pickle.load(f)
    except Exception as e:
        print(f"❌ PKL 加载失败: {e}")
        return {}

def print_adata_stats(adata):
    print("\n" + "-" * 60)
    print("📊 单细胞数据基础统计：")
    print(f"   - 总细胞数: {adata.n_obs:,}")
    print(f"   - 总基因数: {adata.n_vars:,}")
    if 'celltype_manual' in adata.obs.columns:
        print(f"   - 细胞类型数: {len(adata.obs['celltype_manual'].unique())}")
    else:
        print(f"   - 细胞类型数: N/A")
    if 'Group' in adata.obs.columns:
        print(f"   - 分组分布: {adata.obs.Group.value_counts().to_dict()}")
        cnt = adata.obs.Group.value_counts()
        for k, v in cnt.items():
            print(f"   - {k}={v:,}")
    print("-" * 60 + "\n")

def plot_tcmscore_results(adata, score_keys, output_dir):
    output_dir = Path(output_dir)
    valid_keys = [k for k in score_keys if k in adata.obs.columns]
    if not valid_keys:
        print("⚠️ 无有效 TCMscore 可绘制")
        return

    print(f"\n🎨 开始绘制 TCMscore 图表")
    try:
        sc.pl.dotplot(adata, valid_keys, groupby='celltype_manual', show=False)
        plt.savefig(output_dir / "d4c_TCMscore2cell-dotplot_celltype.pdf", bbox_inches='tight')
        plt.close()
        print("   - ✅ 细胞类型点图")
    except Exception as e:
        plt.close()
        print(f"   - ❌ 点图失败: {e}")

    try:
        sc.pl.violin(adata, valid_keys, groupby='Group', show=False)
        plt.savefig(output_dir / "d4c_TCMscore2cell-violin_group.pdf", bbox_inches='tight')
        plt.close()
        print("   - ✅ 分组小提琴图")
    except Exception as e:
        plt.close()
        print(f"   - ❌ 小提琴失败: {e}")

    try:
        sc.pl.umap(adata, color=valid_keys, cmap="OrRd", show=False)
        plt.savefig(output_dir / "d4c_TCMscore2cell-umap_TCMscore.pdf", bbox_inches='tight')
        plt.close()
        print("   - ✅ UMAP 得分图")
    except Exception as e:
        plt.close()
        print(f"   - ❌ UMAP 图失败: {e}")

def TCMscore2cell(
    input_txt, sc_meta_txt, anno_h5ad, formula_config_txt, output_dir, synergy_scale=0.3
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"📂 读取单细胞数据: {anno_h5ad}")
    adata = sc.read(anno_h5ad)
    print_adata_stats(adata)

    print(f"📄 读取输入配置: {input_txt}")
    df_input = pd.read_csv(input_txt, sep="\t", dtype=str)
    print(f"   - 输入样本数: {len(df_input)}")
    print(f"   - 样本列表: {df_input['Name'].tolist()}")

    print(f"📄 读取配方配置: {formula_config_txt}")
    df_form = pd.read_csv(formula_config_txt, sep="\t", dtype=str)
    df_form["JunChen_weight"] = df_form["JunChen_weight"].astype(float)
    df_form["Clinical_ratio"] = df_form["Clinical_ratio"].astype(float)
    print(f"   - 配方配置行数: {len(df_form)}")
    print(f"   - 涉及配方数: {df_form['Formula_name'].nunique()}")

    results = []
    score_keys = []

    for idx, row in df_input.iterrows():
        name = row["Name"].strip()
        typ = row["Type"].strip().lower()
        pkl_path = row["PKL_path"].strip()
        params = row["Parameters"].split(",") if pd.notna(row["Parameters"]) else []
        layer = float(params[0]) if len(params)>=1 and pd.notna(params[0]) else 1.0

        print(f"\n{'='*60}")
        print(f"📌 处理第{idx+1}/{len(df_input)}个样本: {name} | 类型: {typ}")
        print(f"{'='*60}")

        if typ == "ingredient":
            obj = load_pkl(pkl_path)
            print("   - ✅ 成功加载PKL")
            print(f"   - 层级权重: {layer}")
            target_genes = obj.get(name, [])
            valid = [g for g in target_genes if g in adata.var.index]
            print(f"   - 成分靶点数: {len(target_genes)}")
            print(f"   - 最终有效靶点数: {len(valid)}")
            if not valid:
                print("   - ⚠️  无有效基因，跳过")
                continue
            basal = get_basal_score(adata, valid)
            score = basal * (1 + layer)
            key = f"TCMscore_{name}"
            adata.obs[key] = score
            score_keys.append(key)
            results.append({
                "name": name, "type": "ingredient", "ngenes": len(valid),
                "layer_weight": layer, "synergy": 0.0, "final_weight": round(1+layer,2),
                "score_mean": round(score.mean(),4), "score_std": round(score.std(),4)
            })
            print(f"   - 成分得分计算：")
            print(f"     - 基础分均值={basal.mean():.4f}, 基础分标准差={basal.std():.4f}")
            print(f"     - 权重项=(1+{layer})={1+layer:.4f}")
            print(f"     - 最终得分均值={score.mean():.4f}, 标准差={score.std():.4f}")

        elif typ == "herb":
            obj = load_pkl(pkl_path)
            print("   - ✅ 成功加载PKL")
            print(f"   - 层级权重: {layer}")
            if name not in obj:
                print("❌ 草药不在 PKL 中")
                continue
            herb_data = obj[name]
            print(f"   - 草药成分数: {len(herb_data)}")

            all_genes = []
            for cid, tg in herb_data.items():
                all_genes.extend(tg)
            valid = [g for g in all_genes if g in adata.var.index]
            print(f"   - 原始靶点数: {len(all_genes)}, 最终有效靶点数: {len(valid)}")
            if not valid:
                print("⚠️ 无有效基因")
                continue

            basal = get_basal_score(adata, valid)
            syn_ing = calculate_synergy(herb_data, adata, synergy_scale)
            score = basal * (1 + syn_ing) * (1 + layer)
            key = f"TCMscore_{name}"
            adata.obs[key] = score
            score_keys.append(key)

            print(f"   - 单味药得分计算：")
            print(f"     - 基础分: 均值={basal.mean():.4f}, 标准差={basal.std():.4f}")
            print(f"     - 成分协同计算：")
            print(f"       - 有效成分集合数={len([s for s in herb_data.values() if len(set(s)&set(adata.var.index))>0])}")
            print(f"       - 协同系数(原始)={np.mean([jaccard_coeff(set(herb_data[k1]), set(herb_data[k2])) for k1 in herb_data for k2 in herb_data if k1<k2]) if len(herb_data)>=2 else 0:.4f}")
            print(f"       - 协同增益(×{synergy_scale})={syn_ing:.4f}")
            print(f"     - 权重项拆解：(1+{syn_ing:.4f}) × (1+{layer}) = {(1+syn_ing)*(1+layer):.4f}")
            print(f"     - 最终得分: 均值={score.mean():.4f}, 标准差={score.std():.4f}")
            print(f"   - TCMscore: {key} | 均值={score.mean():.3f}")

            results.append({
                "name": name, "type": "herb", "ngenes": len(valid),
                "layer_weight": layer, "synergy": round(syn_ing,4), "final_weight": round((1+syn_ing)*(1+layer),2),
                "score_mean": round(score.mean(),4), "score_std": round(score.std(),4)
            })

        elif typ == "formula":
            formula_obj = load_pkl(pkl_path)
            print("   - ✅ 成功加载复方PKL")
            print(f"   - 层级权重: {layer}")
            if name not in formula_obj:
                print("❌ 复方不在 PKL 中")
                continue

            formula_data = formula_obj[name]
            all_valid_genes = set()
            herb_score_dict = {}
            herb_weight_dict = {}
            herb_gene_sets = []
            herb_basal_stats = {}

            for herb_name, herb_content in formula_data.items():
                all_g = []
                for cid, targets in herb_content.items():
                    vs = [g for g in targets if g in adata.var.index]
                    all_g.extend(vs)
                    all_valid_genes.update(vs)
                if not all_g:
                    print(f"   - ⚠️ 草药 {herb_name} 无有效靶点，跳过")
                    continue

                gset = set(all_g)
                if len(gset) > 0:
                    herb_gene_sets.append(gset)

                print(f"   - 计算草药 {herb_name} 基础分：")
                basal = get_basal_score(adata, all_g)
                syn_ing = calculate_synergy(herb_content, adata, synergy_scale)
                h_basal_score = basal * (1 + syn_ing) * (1 + layer)
                herb_basal_stats[herb_name] = {
                    "basal_mean": basal.mean(),
                    "syn_ing": syn_ing,
                    "h_basal_mean": h_basal_score.mean(),
                    "valid_genes": len(gset)
                }

                sub = df_form[(df_form["Formula_name"] == name) & (df_form["Herb_name"] == herb_name)]
                jc = sub["JunChen_weight"].iloc[0] if not sub.empty else 1.0
                ratio = sub["Clinical_ratio"].iloc[0] if not sub.empty else 1.0
                weight = jc * ratio
                weight_rounded = round(weight, 4)

                herb_score_dict[herb_name] = h_basal_score
                herb_weight_dict[herb_name] = weight

                print(f"   - 草药 {herb_name}:")
                print(f"     - 有效靶点数: {len(gset)}")
                print(f"     - 君臣权重: {jc}, 临床配比: {ratio}, 配伍权重: {weight_rounded}")
                print(f"     - 单味药基础分：均值={basal.mean():.4f}, 成分协同={syn_ing:.4f}")
                print(f"     - 单味药加权分：{(1+syn_ing)*(1+layer):.4f} × 基础分 = {h_basal_score.mean():.4f}")

            if not herb_score_dict:
                print("⚠️ 复方中无有效草药，跳过")
                continue

            total_score = np.zeros(adata.n_obs)
            total_w = 0.0
            weight_sum = 0.0
            weighted_score_sum = 0.0
            print(f"   - 复方加权平均计算：")
            for hn, w in herb_weight_dict.items():
                h_score = herb_score_dict[hn]
                total_score += h_score * w
                total_w += w
                weight_sum += w
                weighted_score_sum += h_score.mean() * w
                print(f"     - {hn}: 得分均值={h_score.mean():.4f} × 权重={w:.4f} = {h_score.mean()*w:.4f}")
            if total_w > 0:
                weighted_avg = total_score / total_w
                weighted_avg_mean = weighted_score_sum / weight_sum
            else:
                weighted_avg = np.zeros(adata.n_obs)
                weighted_avg_mean = 0.0
            print(f"     - 加权总分和: {weighted_score_sum:.4f}, 总权重: {weight_sum:.4f}")
            print(f"     - 加权平均分均值: {weighted_avg_mean:.4f}")

            syn_form = 0.0
            if len(herb_gene_sets) >= 2:
                syn_vals = []
                print(f"   - 复方草药间协同计算（共{len(herb_gene_sets)}味药，{len(herb_gene_sets)*(len(herb_gene_sets)-1)//2}对）：")
                for i in range(len(herb_gene_sets)):
                    for j in range(i+1, len(herb_gene_sets)):
                        herb1 = list(formula_data.keys())[i]
                        herb2 = list(formula_data.keys())[j]
                        inter = len(herb_gene_sets[i] & herb_gene_sets[j])
                        union = len(herb_gene_sets[i] | herb_gene_sets[j])
                        if union > 0:
                            jaccard = inter / union
                            syn_vals.append(jaccard)
                            print(f"     - 草药对({herb1},{herb2}) Jaccard: 交集={inter}, 并集={union}, 系数={jaccard:.4f}")
                if syn_vals:
                    syn_form = np.mean(syn_vals)
                    print(f"     - 所有Jaccard值: {[round(v,4) for v in syn_vals]}")
                    print(f"     - 协同增益均值: 求和={sum(syn_vals):.4f}, 数量={len(syn_vals)}, 均值={syn_form:.4f}")

            final_weight = 1 + layer + syn_form
            final_score = weighted_avg * final_weight
            key = f"TCMscore_{name}"
            adata.obs[key] = final_score
            score_keys.append(key)

            print(f"   - 复方最终得分计算：")
            print(f"     - 加权平均分均值: {weighted_avg.mean():.4f}")
            print(f"     - 单味药间协同增益: {syn_form:.4f}")
            print(f"     - 最终权重项: 1+{layer}+{syn_form:.4f} = {final_weight:.4f}")
            print(f"     - 最终得分: {weighted_avg.mean():.4f} × {final_weight:.4f} = {final_score.mean():.4f}")
            print(f"   - 配方草药数: {len(formula_data)}")
            print(f"   - 最终有效靶点数(并集): {len(all_valid_genes)}")
            print(f"   - TCMscore: {key} | 均值={final_score.mean():.3f}, 标准差={final_score.std():.3f}")

            results.append({
                "name": name, "type": "formula", "ngenes": len(all_valid_genes),
                "layer_weight": layer, "synergy": round(syn_form,4), "final_weight": round(final_weight,2),
                "score_mean": round(final_score.mean(),4), "score_std": round(final_score.std(),4)
            })

    print("\n📊 最终结果：")
    res_df = pd.DataFrame(results)
    print(res_df.to_string(index=False))
    res_df.to_csv(output_dir / "d4c_TCMscore2cell-TCMscore_result.csv", index=False)
    adata.write(output_dir / "d4c_TCMscore2cell-TCMscore_annotated.h5ad")
    plot_tcmscore_results(adata, score_keys, output_dir)
    print(f"\n✅ 全部完成！输出目录: {output_dir}")
    return adata

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Calculate TCMscore for single cell data')
    parser.add_argument('--input_txt', required=True, help='Input configuration file')
    parser.add_argument('--sc_meta_txt', required=True, help='Single cell metadata file')
    parser.add_argument('--anno_h5ad', required=True, help='Annotated h5ad file')
    parser.add_argument('--formula_config_txt', required=True, help='Formula configuration file')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--synergy_scale', type=float, default=0.3, help='Synergy scale factor')
    args = parser.parse_args()

    TCMscore2cell(
        args.input_txt,
        args.sc_meta_txt,
        args.anno_h5ad,
        args.formula_config_txt,
        args.output_dir,
        args.synergy_scale
    )
