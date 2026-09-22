import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D

# ==================== 1. 全局配置与输入文件定义 ====================
WORK_DIR = "/mnt/chuand/jwchen/P_tree"

GENE_FILES = {
    "kinase_CDS": os.path.join(WORK_DIR, "kinase_CDS_genomic_synteny_ete3.csv"),
    "PepSY_25.6": os.path.join(WORK_DIR, "PepSY_25.6_genomic_synteny_ete3.csv"),
    "PepSY":      os.path.join(WORK_DIR, "PepSY_genomic_synteny_ete3.csv"),
    "RR_8.7":     os.path.join(WORK_DIR, "RR_8.7_genomic_synteny_ete3.csv"),
}

GENE_COLORS = {
    "kinase_CDS": "#E41A1C",  # 红色
    "PepSY_25.6": "#377EB8",  # 蓝色
    "PepSY":      "#4DAF4A",  # 绿色
    "RR_8.7":     "#984EA3",  # 紫色
}

# 操纵子间距判定阈值: <= 100 bp (允许负数重叠)
MAX_OPERON_DISTANCE = 100

OUTPUT_CSV = os.path.join(WORK_DIR, "operon_synteny_summary.csv")
OUTPUT_FIG_PNG = os.path.join(WORK_DIR, "operon_structures.png")
OUTPUT_FIG_PDF = os.path.join(WORK_DIR, "operon_structures.pdf")

# ==================== 2. 数据读取与清洗 ====================
print("[*] 正在载入各基因的定位与分类学数据...")
dfs = []
for gene_name, filepath in GENE_FILES.items():
    if not os.path.exists(filepath):
        print(f"[!] 警告: 未找到文件 {filepath}")
        continue
    df = pd.read_csv(filepath, low_memory=False)
    df["target_gene"] = gene_name
    dfs.append(df)

if not dfs:
    print("[x] 错误: 没有读取到任何 CSV 文件，请检查路径！")
    exit(1)

all_df = pd.concat(dfs, ignore_index=True)

all_df = all_df[all_df["assembly_accession"].notna() & (all_df["assembly_accession"] != "Unassembled/WGS")].copy()
all_df = all_df[all_df["cds_start"].notna() & all_df["cds_end"].notna()].copy()
all_df["cds_start"] = all_df["cds_start"].astype(int)
all_df["cds_end"] = all_df["cds_end"].astype(int)
all_df["bitscore"] = pd.to_numeric(all_df["bitscore"], errors="coerce").fillna(0)
all_df["pident"] = pd.to_numeric(all_df["pident"], errors="coerce").fillna(0)

print(f"[+] 成功整合数据记录: {len(all_df)} 条，涉及菌株 Assembly 数: {all_df['assembly_accession'].nunique()}")

# ==================== 3. 操纵子物理聚合分析 (<= 100 bp) ====================
print(f"[*] 正在聚类操纵子物理结构（严格阈值: 相邻基因间距 <= {MAX_OPERON_DISTANCE} bp）...")

strain_clusters = []

for (asm, contig), group in all_df.groupby(["assembly_accession", "replicon_contig"]):
    group = group.sort_values(by="cds_start")
    
    current_cluster = []
    for _, row in group.iterrows():
        if not current_cluster:
            current_cluster.append(row)
        else:
            prev_row = current_cluster[-1]
            distance = row["cds_start"] - prev_row["cds_end"]
            if distance <= MAX_OPERON_DISTANCE:
                current_cluster.append(row)
            else:
                strain_clusters.append(pd.DataFrame(current_cluster))
                current_cluster = [row]
    if current_cluster:
        strain_clusters.append(pd.DataFrame(current_cluster))

best_clusters = {}
for cluster_df in strain_clusters:
    asm = cluster_df["assembly_accession"].iloc[0]
    unique_genes = cluster_df["target_gene"].nunique()
    total_bitscore = cluster_df["bitscore"].sum()
    
    if asm not in best_clusters:
        best_clusters[asm] = cluster_df
    else:
        curr_best = best_clusters[asm]
        curr_unique = curr_best["target_gene"].nunique()
        curr_bitscore = curr_best["bitscore"].sum()
        if (unique_genes > curr_unique) or (unique_genes == curr_unique and total_bitscore > curr_bitscore):
            best_clusters[asm] = cluster_df

all_asms = set(all_df["assembly_accession"].unique())
found_asms = set(best_clusters.keys())
for missing_asm in all_asms - found_asms:
    best_clusters[missing_asm] = all_df[all_df["assembly_accession"] == missing_asm].head(1)

# ==================== 4. 构建汇总统计表格 ====================
print("[*] 正在构建操纵子分类统计表...")
summary_rows = []

for asm, c_df in best_clusters.items():
    c_df = c_df.sort_values(by="cds_start")
    
    org = c_df["organism"].dropna().iloc[0] if "organism" in c_df.columns and not c_df["organism"].dropna().empty else "Unknown"
    taxid = str(c_df["taxid"].dropna().iloc[0]) if "taxid" in c_df.columns and not c_df["taxid"].dropna().empty else "N/A"
    leaf_label = c_df["tree_leaf_label"].dropna().iloc[0] if "tree_leaf_label" in c_df.columns and not c_df["tree_leaf_label"].dropna().empty else f"{taxid}_{asm}"
    phylum = c_df["phylum"].dropna().iloc[0] if "phylum" in c_df.columns and not c_df["phylum"].dropna().empty else "Unclassified"
    
    genes_present = list(dict.fromkeys(c_df["target_gene"].tolist()))
    gene_count = len(genes_present)
    
    strands = c_df["strand"].tolist()
    is_strand_consistent = (len(set(strands)) == 1) if strands else False
    
    gene_order_str = " -> ".join([f"{row['target_gene']}({row['strand']})" for _, row in c_df.iterrows()])
    
    row_data = {
        "assembly_accession": asm,
        "organism": org,
        "taxid": taxid,
        "tree_leaf_label": leaf_label,
        "phylum": phylum,
        "gene_count": gene_count,
        "present_genes": ";".join(genes_present),
        "operon_contig": c_df["replicon_contig"].iloc[0],
        "operon_start": c_df["cds_start"].min(),
        "operon_end": c_df["cds_end"].max(),
        "operon_span_bp": c_df["cds_end"].max() - c_df["cds_start"].min(),
        "operon_strand_consistent": is_strand_consistent,
        "gene_order": gene_order_str
    }
    
    for gname in GENE_FILES.keys():
        g_hits = c_df[c_df["target_gene"] == gname]
        if not g_hits.empty:
            g_hit = g_hits.sort_values(by="bitscore", ascending=False).iloc[0]
            row_data[f"{gname}_protein_acc"] = g_hit["protein_accession"]
            row_data[f"{gname}_pident"] = g_hit["pident"]
            row_data[f"{gname}_evalue"] = g_hit["evalue"]
            row_data[f"{gname}_bitscore"] = g_hit["bitscore"]
            row_data[f"{gname}_start"] = g_hit["cds_start"]
            row_data[f"{gname}_end"] = g_hit["cds_end"]
            row_data[f"{gname}_strand"] = g_hit["strand"]
        else:
            row_data[f"{gname}_protein_acc"] = None
            row_data[f"{gname}_pident"] = np.nan
            row_data[f"{gname}_evalue"] = np.nan
            row_data[f"{gname}_bitscore"] = np.nan
            row_data[f"{gname}_start"] = np.nan
            row_data[f"{gname}_end"] = np.nan
            row_data[f"{gname}_strand"] = None

    summary_rows.append(row_data)

summary_df = pd.DataFrame(summary_rows)
summary_df.sort_values(by=["gene_count", "operon_span_bp"], ascending=[False, True], inplace=True)
summary_df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")

print(f"\n[✓] 操纵子共线性特征分析完成！结果已保存至: {OUTPUT_CSV}")
print("\n【基因数量分类统计结果 (<= 100 bp)】")
print(summary_df["gene_count"].value_counts().sort_index(ascending=False).rename(lambda x: f"含 {x} 个基因的菌株数"))

# ==================== 5. 绘制操纵子结构图 ====================
print("\n[*] 正在绘制基因排布示意图...")
SAMPLE_PER_CLASS = 6
plot_asms = []
for count in [4, 3, 2, 1]:
    sub = summary_df[summary_df["gene_count"] == count]
    if not sub.empty:
        plot_asms.extend(sub.head(SAMPLE_PER_CLASS)["assembly_accession"].tolist())

num_strains = len(plot_asms)
fig, ax = plt.subplots(figsize=(14, max(6, num_strains * 0.5) + 1.5))

y_pos = 0
y_labels = []
y_ticks = []

for asm in plot_asms:
    cluster = best_clusters[asm].sort_values(by="cds_start")
    min_x = cluster["cds_start"].min()
    max_x = cluster["cds_end"].max()
    span = max(max_x - min_x, 3000)
    ax.plot([0, span], [y_pos, y_pos], color="#CCCCCC", lw=1.5, zorder=1)
    
    for _, row in cluster.iterrows():
        g_name = row["target_gene"]
        start_rel = row["cds_start"] - min_x
        end_rel = row["cds_end"] - min_x
        width = end_rel - start_rel
        color = GENE_COLORS.get(g_name, "#999999")
        strand = row["strand"]
        
        head_width = width * 0.3 if width > 300 else width * 0.5
        y_bottom = y_pos - 0.25
        y_top = y_pos + 0.25
        y_mid = y_pos
        
        if strand == "-":
            pts = [
                [start_rel + head_width, y_top],
                [end_rel, y_top],
                [end_rel, y_bottom],
                [start_rel + head_width, y_bottom],
                [start_rel, y_mid]
            ]
        else:
            pts = [
                [start_rel, y_top],
                [end_rel - head_width, y_top],
                [end_rel, y_mid],
                [end_rel - head_width, y_bottom],
                [start_rel, y_bottom]
            ]
        poly = patches.Polygon(pts, closed=True, facecolor=color, edgecolor="black", lw=0.6, zorder=2)
        ax.add_patch(poly)
        
        pident_val = row["pident"]
        ax.text(start_rel + width / 2, y_pos + 0.32, f"{g_name}\n({pident_val:.1f}%)", 
                ha="center", va="bottom", fontsize=7, color="#222222")

    org_name = summary_df.loc[summary_df["assembly_accession"] == asm, "organism"].values[0]
    short_org = " ".join(str(org_name).split()[:2])
    g_cnt = summary_df.loc[summary_df["assembly_accession"] == asm, "gene_count"].values[0]
    
    y_labels.append(f"[{g_cnt} Genes] {asm}\n({short_org})")
    y_ticks.append(y_pos)
    y_pos -= 1.0

ax.set_yticks(y_ticks)
ax.set_yticklabels(y_labels, fontsize=8)
ax.set_xlabel("Relative Genomic Position (bp)", fontsize=11, fontweight="bold")
ax.set_title(f"Operon Synteny Organization (Intergenic Distance <= {MAX_OPERON_DISTANCE} bp)", fontsize=13, fontweight="bold", pad=15)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.grid(axis="x", linestyle="--", alpha=0.4)

legend_elements = [
    Line2D([0], [0], marker='s', color='w', label=f"{g}", markerfacecolor=c, markersize=10, markeredgecolor='black')
    for g, c in GENE_COLORS.items()
]
ax.legend(handles=legend_elements, loc="upper right", frameon=True, fontsize=9, title="Homologous Genes")

plt.tight_layout()
plt.savefig(OUTPUT_FIG_PNG, dpi=300)
plt.savefig(OUTPUT_FIG_PDF)
plt.close()

print(f"[✓] 基因排布示意图重新绘制完成！")
print(f" -> PNG 预览图: {OUTPUT_FIG_PNG}")
print(f" -> PDF 矢量图: {OUTPUT_FIG_PDF}")