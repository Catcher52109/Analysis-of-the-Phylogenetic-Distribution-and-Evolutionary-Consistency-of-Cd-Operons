import os
import sys
import pandas as pd
import numpy as np

# ==================== 0. ETE3 环境兼容与初始化 ====================
try:
    import cgi
except ImportError:
    import html
    import types
    fake_cgi = types.ModuleType("cgi")
    fake_cgi.escape = html.escape
    sys.modules["cgi"] = fake_cgi

try:
    from ete3 import NCBITaxa, Tree
    LOCAL_TAXDUMP_DIR = "/mnt/chuand/Acr/nr/taxdump"
    tar_path = os.path.join(LOCAL_TAXDUMP_DIR, "taxdump.tar.gz")
    if os.path.exists(tar_path):
        ncbi = NCBITaxa(taxdump_file=tar_path)
    else:
        ncbi = NCBITaxa()
except Exception as e:
    print(f"[x] ETE3 加载失败: {e}")
    sys.exit(1)

# ==================== 1. 文件与参数配置 ====================
WORK_DIR = "/mnt/chuand/jwchen/P_tree"
SUMMARY_CSV = os.path.join(WORK_DIR, "operon_synteny_summary.csv")
OUTPUT_NWK = os.path.join(WORK_DIR, "operon_strains_taxonomy.nwk")
OUTPUT_TREE_MAPPING = os.path.join(WORK_DIR, "tree_metadata_mapping.csv")

# 抽样策略：
# 4 基因（1521 株）、3 基因（1171 株）、2 基因（1503 株）全部纳入！
# 1 基因仅抽取代表株作为系统发育进化背景对照
MAX_SINGLE_GENE_SAMPLES = 500

# ==================== 2. 读取第一阶段汇总表 ====================
print("[*] 正在载入 100 bp 筛选后的操纵子汇总表...")
if not os.path.exists(SUMMARY_CSV):
    print(f"[x] 找不到文件: {SUMMARY_CSV}")
    sys.exit(1)

df = pd.read_csv(SUMMARY_CSV, low_memory=False)
df = df[df["taxid"].notna() & (df["taxid"] != "N/A") & (df["taxid"] != "None")].copy()
df["taxid_clean"] = df["taxid"].astype(str).str.split(".").str[0].astype(int)

# 保留全部多基因操纵子
df_multi = df[df["gene_count"] >= 2].copy()
df_single = df[df["gene_count"] == 1].copy()

if len(df_single) > MAX_SINGLE_GENE_SAMPLES:
    df_single_sampled = df_single.sample(n=MAX_SINGLE_GENE_SAMPLES, random_state=42)
else:
    df_single_sampled = df_single

tree_df = pd.concat([df_multi, df_single_sampled], ignore_index=True)
print(f"[+] 纳入进化树构建的菌株总数: {len(tree_df)}")
print(f" -> 包含 4 基因高保真操纵子: {len(tree_df[tree_df['gene_count'] == 4])} 株")
print(f" -> 包含 3 基因操纵子: {len(tree_df[tree_df['gene_count'] == 3])} 株")
print(f" -> 包含 2 基因操纵子: {len(tree_df[tree_df['gene_count'] == 2])} 株")
print(f" -> 孤立单基因背景对照: {len(df_single_sampled)} 株")

# ==================== 3. 构造系统发育拓扑 ====================
print("[*] 正在基于本地 NCBI Taxonomy 生成系统发育树骨架...")
unique_taxids = tree_df["taxid_clean"].unique().tolist()

try:
    topology = ncbi.get_topology(unique_taxids)
except Exception as e:
    print(f"[!] 批量解析拓扑异常，正在逐一剔除不可识别节点: {e}")
    valid_taxids = []
    for t in unique_taxids:
        try:
            ncbi.get_lineage(t)
            valid_taxids.append(t)
        except Exception:
            continue
    topology = ncbi.get_topology(valid_taxids)
    tree_df = tree_df[tree_df["taxid_clean"].isin(valid_taxids)].copy()

# ==================== 4. 扩展叶节点为 taxID_assembly_accession ====================
print("[*] 正在将树节点扩展为真实菌株级别 (taxID_assembly_accession)...")

taxid_to_strains = {}
for _, row in tree_df.iterrows():
    tid = row["taxid_clean"]
    asm = row["assembly_accession"]
    leaf_id = f"{tid}_{asm}"
    taxid_to_strains.setdefault(tid, []).append(leaf_id)

tree_copy = topology.copy()

for leaf in tree_copy.get_leaves():
    tid = int(leaf.name)
    matched_leaves = taxid_to_strains.get(tid, [])
    
    if not matched_leaves:
        leaf.detach()
    elif len(matched_leaves) == 1:
        leaf.name = matched_leaves[0]
    else:
        parent = leaf.up
        leaf.detach()
        sub_node = parent.add_child(name=f"taxid_{tid}")
        for m_leaf in matched_leaves:
            sub_node.add_child(name=m_leaf)

# 清理冗余单子代内部节点
for node in tree_copy.traverse():
    if not node.is_leaf() and not node.is_root() and len(node.children) == 1:
        child = node.children[0]
        parent = node.up
        child.detach()
        node.detach()
        parent.add_child(child)

# 保存标准 Newick 树
tree_copy.write(outfile=OUTPUT_NWK, format=1)
print(f"[✓] 系统发育进化树构建完成！标准 NWK 保存至: {OUTPUT_NWK}")
print(f" -> 最终进化树叶节点总数: {len(tree_copy.get_leaves())}")

# ==================== 5. 导出与树严格对应的元数据映射表 ====================
print("[*] 正在导出与进化树严格对齐的元数据表...")
tree_leaves = set(tree_copy.get_leaf_names())

tree_df["tree_leaf_node_id"] = tree_df["taxid_clean"].astype(str) + "_" + tree_df["assembly_accession"]
tree_mapping_df = tree_df[tree_df["tree_leaf_node_id"].isin(tree_leaves)].copy()

front_cols = [
    "tree_leaf_node_id", "assembly_accession", "taxid", "organism", "phylum", 
    "gene_count", "present_genes", "operon_strand_consistent", "gene_order"
]
other_cols = [c for c in tree_mapping_df.columns if c not in front_cols]
tree_mapping_df = tree_mapping_df[front_cols + other_cols]

tree_mapping_df.to_csv(OUTPUT_TREE_MAPPING, index=False, encoding="utf-8-sig")
print(f"[✓] 进化树对接元数据表已保存至: {OUTPUT_TREE_MAPPING}")