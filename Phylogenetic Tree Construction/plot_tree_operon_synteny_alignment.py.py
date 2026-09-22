import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D

# ==================== 0. ETE3 环境初始化 ====================
try:
    import cgi
except ImportError:
    import html
    import types
    fake_cgi = types.ModuleType("cgi")
    fake_cgi.escape = html.escape
    sys.modules["cgi"] = fake_cgi

from ete3 import Tree

# ==================== 1. 文件路径配置 ====================
WORK_DIR = "/mnt/chuand/jwchen/P_tree"
NWK_FILE = os.path.join(WORK_DIR, "operon_strains_taxonomy.nwk")
MAPPING_CSV = os.path.join(WORK_DIR, "tree_metadata_mapping.csv")

OUT_FIG_PNG = os.path.join(WORK_DIR, "aligned_tree_operons_synteny.png")
OUT_FIG_PDF = os.path.join(WORK_DIR, "aligned_tree_operons_synteny.pdf")
ITOL_DIR = os.path.join(WORK_DIR, "itol_annotations")
os.makedirs(ITOL_DIR, exist_ok=True)

# 统一基因颜色
GENE_COLORS = {
    "kinase_CDS": "#E41A1C",  # 红色
    "PepSY_25.6": "#377EB8",  # 蓝色
    "PepSY":      "#4DAF4A",  # 绿色
    "RR_8.7":     "#984EA3",  # 紫色
}
GENE_NAMES = ["kinase_CDS", "PepSY_25.6", "PepSY", "RR_8.7"]

# ==================== 2. 读取数据与系统发育树 ====================
print("[*] 正在载入进化树与菌株共线性特征映射表...")
if not os.path.exists(NWK_FILE) or not os.path.exists(MAPPING_CSV):
    print("[x] 找不到 NWK 树文件或元数据表！")
    sys.exit(1)

meta_df = pd.read_csv(MAPPING_CSV, low_memory=False)
meta_dict = {row["tree_leaf_node_id"]: row for _, row in meta_df.iterrows()}

t = Tree(NWK_FILE, format=1)
leaf_names = [leaf.name for leaf in t.get_leaves()]
print(f"[+] 树包含叶节点数: {len(leaf_names)}")

# ==================== 3. 生成 iTOL 交互式全量注释配置文件 ====================
print(f"[*] 正在为全部 {len(leaf_names)} 株生成 iTOL 在线可视化注释配置集...")

tab_char = "\t"
gene_labels_str = tab_char.join(GENE_NAMES)
gene_colors_str = tab_char.join([GENE_COLORS[g] for g in GENE_NAMES])
gene_shapes_str = tab_char.join(["1" for _ in GENE_NAMES])

# 3.1 基因数量分类色条 (iTOL ColorStrip)
itol_colors = {4: "#E41A1C", 3: "#FF7F00", 2: "#377EB8", 1: "#999999"}
with open(os.path.join(ITOL_DIR, "itol_gene_count_strip.txt"), "w") as f:
    f.write("DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tGene_Count_Class\nCOLOR\t#E41A1C\nDATA\n")
    for leaf in leaf_names:
        if leaf in meta_dict:
            cnt = int(meta_dict[leaf]["gene_count"])
            f.write(f"{leaf}\t{itol_colors.get(cnt, '#CCCCCC')}\t{cnt}_Genes\n")

# 3.2 四基因存在/缺失二元矩阵 (iTOL Binary Matrix)
with open(os.path.join(ITOL_DIR, "itol_operon_binary_matrix.txt"), "w") as f:
    f.write("DATASET_BINARY\nSEPARATOR TAB\nDATASET_LABEL\tOperon_Genes_Presence\n")
    f.write(f"FIELD_LABELS\t{gene_labels_str}\n")
    f.write(f"FIELD_COLORS\t{gene_colors_str}\n")
    f.write(f"FIELD_SHAPES\t{gene_shapes_str}\nDATA\n")
    for leaf in leaf_names:
        if leaf in meta_dict:
            row = meta_dict[leaf]
            p_genes = str(row.get("present_genes", "")).split(";")
            vals = ["1" if g in p_genes else "0" for g in GENE_NAMES]
            vals_str = tab_char.join(vals)
            f.write(f"{leaf}\t{vals_str}\n")

# 3.3 各基因相似度热图 (iTOL Heatmap)
with open(os.path.join(ITOL_DIR, "itol_pident_heatmap.txt"), "w") as f:
    f.write("DATASET_HEATMAP\nSEPARATOR TAB\nDATASET_LABEL\tGene_Identity_Pident\n")
    f.write(f"FIELD_LABELS\t{gene_labels_str}\n")
    f.write("COLOR_MIN\t#FFFFFF\nCOLOR_MAX\t#08519C\nDATA\n")
    for leaf in leaf_names:
        if leaf in meta_dict:
            row = meta_dict[leaf]
            pident_vals = []
            for g in GENE_NAMES:
                val = row.get(f"{g}_pident", np.nan)
                pident_vals.append(f"{val:.1f}" if pd.notna(val) else "-1")
            pident_str = tab_char.join(pident_vals)
            f.write(f"{leaf}\t{pident_str}\n")

print(f"[✓] iTOL 注释文件已输出至: {ITOL_DIR}")

# ==================== 4. 提取代表性株构建精细对齐主图 ====================
print("[*] 正在构建代表性株系统发育树与操纵子结构对齐图...")

SUB_SAMPLE_PER_GROUP = 12
sampled_leaves = []
for cnt in [4, 3, 2, 1]:
    group_df = meta_df[meta_df["gene_count"] == cnt]
    if not group_df.empty:
        sampled = group_df.groupby("phylum", group_keys=False).apply(lambda x: x.head(3))
        if len(sampled) > SUB_SAMPLE_PER_GROUP:
            sampled = sampled.head(SUB_SAMPLE_PER_GROUP)
        sampled_leaves.extend(sampled["tree_leaf_node_id"].tolist())

sampled_leaves = [l for l in sampled_leaves if l in set(leaf_names)]
sub_tree = t.copy()
sub_tree.prune(sampled_leaves, preserve_branch_length=False)

ordered_leaves = [leaf.name for leaf in sub_tree.get_leaves()]
n_strains = len(ordered_leaves)

# ==================== 5. 联合图绘制 (Matplotlib 树 + 共线性) ====================
fig = plt.figure(figsize=(18, max(8, n_strains * 0.45)))
gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 2.5], wspace=0.15)

ax_tree = fig.add_subplot(gs[0])
ax_operon = fig.add_subplot(gs[1])

y_coords = {name: i for i, name in enumerate(reversed(ordered_leaves))}

node_coords = {}
def calc_coords(node, current_x=0):
    if node.is_leaf():
        node_coords[node] = (current_x, y_coords[node.name])
    else:
        children_y = []
        for child in node.children:
            calc_coords(child, current_x + 1)
            children_y.append(node_coords[child][1])
        node_coords[node] = (current_x, np.mean(children_y))

calc_coords(sub_tree)

def draw_branches(node):
    x_curr, y_curr = node_coords[node]
    for child in node.children:
        x_child, y_child = node_coords[child]
        ax_tree.plot([x_curr, x_curr], [y_curr, y_child], color="#444444", lw=1.2)
        ax_tree.plot([x_curr, x_child], [y_child, y_child], color="#444444", lw=1.2)
        draw_branches(child)

draw_branches(sub_tree)

ax_tree.set_ylim(-0.5, n_strains - 0.5)
ax_tree.set_xlim(-0.5, max([coord[0] for coord in node_coords.values()]) + 0.5)
ax_tree.axis("off")
ax_tree.set_title("Phylogenetic Taxonomy", fontsize=12, fontweight="bold", pad=10)

y_labels = []
y_ticks = []

for leaf_id in ordered_leaves:
    y_pos = y_coords[leaf_id]
    row = meta_dict[leaf_id]
    
    org = str(row.get("organism", "Unknown"))
    short_org = " ".join(org.split()[:2])
    g_cnt = int(row.get("gene_count", 0))
    y_labels.append(f"[{g_cnt}G] {short_org}")
    y_ticks.append(y_pos)
    
    ax_operon.plot([0, 5000], [y_pos, y_pos], color="#E0E0E0", lw=2, zorder=1)
    
    gene_records = []
    for gname in GENE_NAMES:
        start_val = row.get(f"{gname}_start")
        end_val = row.get(f"{gname}_end")
        if pd.notna(start_val) and pd.notna(end_val):
            gene_records.append({
                "gene": gname,
                "start": int(start_val),
                "end": int(end_val),
                "strand": str(row.get(f"{gname}_strand", "+")),
                "pident": float(row.get(f"{gname}_pident", 0))
            })
            
    if gene_records:
        gene_records.sort(key=lambda x: x["start"])
        min_base = gene_records[0]["start"]
        
        for g in gene_records:
            s_rel = g["start"] - min_base
            e_rel = g["end"] - min_base
            width = max(e_rel - s_rel, 200)
            color = GENE_COLORS.get(g["gene"], "#999999")
            strand = g["strand"]
            
            head_width = width * 0.35 if width > 300 else width * 0.5
            y_bot = y_pos - 0.22
            y_tp = y_pos + 0.22
            y_md = y_pos
            
            if strand == "-":
                pts = [
                    [s_rel + head_width, y_tp],
                    [e_rel, y_tp],
                    [e_rel, y_bot],
                    [s_rel + head_width, y_bot],
                    [s_rel, y_md]
                ]
            else:
                pts = [
                    [s_rel, y_tp],
                    [e_rel - head_width, y_tp],
                    [e_rel, y_md],
                    [e_rel - head_width, y_bot],
                    [s_rel, y_bot]
                ]
            poly = patches.Polygon(pts, closed=True, facecolor=color, edgecolor="black", lw=0.6, zorder=2)
            ax_operon.add_patch(poly)
            
            ax_operon.text(s_rel + width / 2, y_pos + 0.26, f"{g['pident']:.0f}%", 
                           ha="center", va="bottom", fontsize=6.5, color="#111111")

ax_operon.set_yticks(y_ticks)
ax_operon.set_yticklabels(y_labels, fontsize=8)
ax_operon.set_ylim(-0.5, n_strains - 0.5)
ax_operon.set_xlim(-200, 5200)
ax_operon.set_xlabel("Relative Physical Genomic Distance (bp)", fontsize=10, fontweight="bold")
ax_operon.set_title("Operon Synteny & Gene Identity Alignment", fontsize=12, fontweight="bold", pad=10)
ax_operon.spines["top"].set_visible(False)
ax_operon.spines["right"].set_visible(False)
ax_operon.grid(axis="x", linestyle="--", alpha=0.3)

legend_elements = [
    Line2D([0], [0], marker='s', color='w', label=g, markerfacecolor=c, markersize=9, markeredgecolor='black')
    for g, c in GENE_COLORS.items()
]
ax_operon.legend(handles=legend_elements, loc="upper right", frameon=True, fontsize=8, title="Genes")

plt.tight_layout()
plt.savefig(OUT_FIG_PNG, dpi=300)
plt.savefig(OUT_FIG_PDF)
plt.close()

print(f"[✓] 进化树与操纵子共线性联合对齐图绘制完成！")
print(f" -> PNG 图: {OUT_FIG_PNG}")
print(f" -> PDF 矢量图: {OUT_FIG_PDF}")