import os
import re
import sys
import time
import subprocess
import pandas as pd
from Bio import Entrez
from xml.etree import ElementTree as ET

# ==================== 0. ETE3 本地分类学数据库兼容与初始化 ====================
try:
    import cgi
except ImportError:
    import html
    import types
    fake_cgi = types.ModuleType("cgi")
    fake_cgi.escape = html.escape
    sys.modules["cgi"] = fake_cgi

try:
    from ete3 import NCBITaxa
    LOCAL_TAXDUMP_DIR = "/mnt/chuand/Acr/nr/taxdump"
    tar_path = os.path.join(LOCAL_TAXDUMP_DIR, "taxdump.tar.gz")
    if os.path.exists(tar_path):
        ncbi_taxa = NCBITaxa(taxdump_file=tar_path)
    else:
        ncbi_taxa = NCBITaxa()
    ETE3_AVAILABLE = True
except Exception as e:
    print(f"[!] ETE3 模块初始化提示: {e}")
    ETE3_AVAILABLE = False


# ==============================================================================
#                      【当前蛋白配置: PepSY_25.6】
# ==============================================================================
PROTEIN_NAME = "PepSY_25.6"
TARGET_DIR_NAME = "PepSY_25.6"

AMINO_ACID_SEQ = """
MKTLTALTLASIIGLTAGTVHARDLGPDEALRLRDAGTIVSFEKLNATALAKHPGSTITDTELEEQYGKYIYQIELRDPQGLEWDLELDAVSGQVLKDHQDT
"""

# 本地数据库与比对参数配置
NR_DB_PATH = "/mnt/chuand/Acr/nr/nr"
THREADS = 16                                 # CPU 核心数
ITERATIONS = 3                               # PSI-BLAST 迭代轮数
EVALUE_THRESH = 0.005                        # 纳入迭代的阈值

# 目录与输出路径自动绑定
BASE_WORK_DIR = "/mnt/chuand/jwchen"
SUB_DIR = os.path.join(BASE_WORK_DIR, TARGET_DIR_NAME)
os.makedirs(SUB_DIR, exist_ok=True)

QUERY_FASTA = os.path.join(SUB_DIR, f"{PROTEIN_NAME}.fasta")
BLAST_OUTPUT_TSV = os.path.join(SUB_DIR, f"{PROTEIN_NAME}_psiblast_hits.tsv")
OUTPUT_CSV = os.path.join(SUB_DIR, f"{PROTEIN_NAME}_genomic_synteny_ete3.csv")

# NCBI 邮箱配置（用于 IPG 物理坐标解析）
Entrez.email = "your_email@example.com"
Entrez.api_key = os.getenv("NCBI_API_KEY", None)
# ==============================================================================


def clean_accession(acc_str):
    """提取标准蛋白质 Accession (如 ref|WP_010982341.1| -> WP_010982341.1)"""
    acc_str = str(acc_str).strip()
    if "|" in acc_str:
        parts = [p for p in acc_str.split("|") if p]
        for p in reversed(parts):
            if not p.isdigit() and len(p) > 3:
                return p
    return acc_str

def get_base_acc(acc_str):
    return re.sub(r"\.\d+$", "", clean_accession(acc_str))


def run_local_psiblast(query_fa, out_tsv, db_path, num_iterations=3, threads=16, evalue=0.005):
    """调用本地多线程 PSI-BLAST，执行 3 轮矩阵迭代"""
    print(f"[*] 启动本地 PSI-BLAST（蛋白: {PROTEIN_NAME}, 数据库: {db_path}, 迭代: {num_iterations} 次, 线程: {threads}）...")
    
    db_dir = os.path.dirname(db_path)
    env = os.environ.copy()
    env["BLASTDB"] = f"{db_dir}:{env.get('BLASTDB', '')}"
    
    outfmt = "6 sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle"
    
    cmd = [
        "psiblast",
        "-query", query_fa,
        "-db", db_path,
        "-num_iterations", str(num_iterations),
        "-inclusion_ethresh", str(evalue),
        "-evalue", str(evalue),
        "-num_threads", str(threads),
        "-outfmt", outfmt,
        "-out", out_tsv
    ]
    
    res = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"[x] 本地 PSI-BLAST 运行失败:\n{res.stderr}")
        sys.exit(1)
        
    print(f"[+] 本地 PSI-BLAST 运行完成，比对结果已保存至: {out_tsv}")
    
    cols = ["sacc", "pident", "length", "mismatch", "gapopen", 
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxids", "stitle"]
    df = pd.read_csv(out_tsv, sep="\t", header=None, names=cols)
    df = df.drop_duplicates(subset=["sacc"])
    print(f"[+] 提取到 {len(df)} 个同源蛋白唯一命中。")
    return df


def fetch_ipg_metadata(blast_df, batch_size=50, max_retries=3):
    """根据命中同源蛋白批量检索 NCBI IPG 基因组物理定位（Contig, 坐标, 正负链）"""
    blast_dict = {}
    for _, row in blast_df.iterrows():
        c_acc = clean_accession(row["sacc"])
        b_acc = get_base_acc(c_acc)
        score_info = {
            "query_protein_acc": c_acc,
            "pident": row["pident"],
            "evalue": row["evalue"],
            "bitscore": row["bitscore"],
            "staxid": str(row["staxids"]).split(";")[0] if pd.notna(row["staxids"]) else None
        }
        blast_dict[c_acc] = score_info
        blast_dict[b_acc] = score_info

    unique_accs = list(set([clean_accession(a) for a in blast_df["sacc"] if clean_accession(a)]))
    total = len(unique_accs)
    records = []
    print(f"[*] 开始通过 IPG 检索 {total} 个同源蛋白在全基因组中的物理定位...")

    for i in range(0, total, batch_size):
        batch = unique_accs[i:i + batch_size]
        for attempt in range(1, max_retries + 1):
            try:
                handle = Entrez.efetch(db="protein", id=",".join(batch), rettype="ipg", retmode="xml")
                xml_data = handle.read()
                handle.close()
                
                root = ET.fromstring(xml_data)
                reports = root.findall(".//IPGReport")
                if not reports and root.tag.endswith("IPGReport"):
                    reports = [root]

                for report in reports:
                    parsed_prot_acc = report.get("accver") or report.get("product_accver")
                    if not parsed_prot_acc:
                        for tag in [".//Product", ".//Protein"]:
                            node = report.find(tag)
                            if node is not None:
                                parsed_prot_acc = node.get("accver") or node.get("accession")
                                if parsed_prot_acc:
                                    break
                    
                    matched_scores = {}
                    final_protein_acc = None
                    if parsed_prot_acc:
                        final_protein_acc = parsed_prot_acc
                        matched_scores = blast_dict.get(parsed_prot_acc) or blast_dict.get(get_base_acc(parsed_prot_acc), {})
                    
                    if not matched_scores:
                        for b_id in batch:
                            if b_id in xml_data.decode("utf-8", errors="ignore"):
                                final_protein_acc = b_id
                                matched_scores = blast_dict.get(b_id, {})
                                break
                    
                    if not final_protein_acc:
                        final_protein_acc = batch[0] if len(batch) == 1 else "Unknown"

                    cds_nodes = report.findall(".//CDS")
                    for cds in cds_nodes:
                        assembly_id = cds.get("assembly")
                        taxid = cds.get("taxid") or matched_scores.get("staxid")
                        records.append({
                            "protein_accession": matched_scores.get("query_protein_acc", final_protein_acc),
                            "pident": matched_scores.get("pident", None),
                            "evalue": matched_scores.get("evalue", None),
                            "bitscore": matched_scores.get("bitscore", None),
                            "organism": cds.get("org"),
                            "taxid": str(taxid) if taxid else None,
                            "assembly_accession": assembly_id if assembly_id else "Unassembled/WGS",
                            "assembly_url": f"https://www.ncbi.nlm.nih.gov/datasets/genome/{assembly_id}/" if assembly_id else "N/A",
                            "replicon_contig": cds.get("accver"),
                            "cds_start": int(cds.get("start")) if cds.get("start") else None,
                            "cds_end": int(cds.get("stop")) if cds.get("stop") else None,
                            "strand": cds.get("strand")
                        })
                break
            except Exception as e:
                time.sleep(attempt * 2)
                
        time.sleep(0.3)
        print(f" -> 进度: 已解析 {min(i + batch_size, total)}/{total} 个蛋白，已提取 {len(records)} 条定位数据", end="\r")

    print("\n")
    return pd.DataFrame(records)


def get_ete3_lineage(taxid_series):
    """解析系统发育分类层级与进化树标注"""
    if not ETE3_AVAILABLE:
        return pd.DataFrame()

    print("[*] 正在解析 ETE3 系统发育层级与进化树标注...")
    lineage_records = []
    unique_taxids = taxid_series.dropna().unique()

    for tid in unique_taxids:
        try:
            tid_int = int(float(tid))
            lineage = ncbi_taxa.get_lineage(tid_int)
            ranks = ncbi_taxa.get_rank(lineage)
            names = ncbi_taxa.get_taxid_translator(lineage)

            rank_dict = {rank: names.get(t) for t, rank in ranks.items()}
            genus = rank_dict.get("genus", "Unclassified")
            species = rank_dict.get("species", "sp")
            
            lineage_records.append({
                "taxid": str(tid_int),
                "phylum": rank_dict.get("phylum", "Unclassified"),
                "class": rank_dict.get("class", "Unclassified"),
                "order": rank_dict.get("order", "Unclassified"),
                "family": rank_dict.get("family", "Unclassified"),
                "genus": genus,
                "species": species,
                "tree_leaf_label": f"{tid_int}_{genus}_{species}".replace(" ", "_")
            })
        except Exception:
            lineage_records.append({"taxid": str(tid)})

    return pd.DataFrame(lineage_records)


if __name__ == "__main__":
    clean_seq = "".join(AMINO_ACID_SEQ.split())
    if not clean_seq:
        print("[x] 氨基酸序列为空，请检查！")
        sys.exit(1)

    # 1. 写入子目录下的 FASTA
    with open(QUERY_FASTA, "w", encoding="utf-8") as f:
        f.write(f">{PROTEIN_NAME}\n{clean_seq}\n")

    # 2. 本地多线程 PSI-BLAST 迭代 3 轮
    blast_hits_df = run_local_psiblast(
        query_fa=QUERY_FASTA,
        out_tsv=BLAST_OUTPUT_TSV,
        db_path=NR_DB_PATH,
        num_iterations=ITERATIONS,
        threads=THREADS,
        evalue=EVALUE_THRESH
    )

    # 3. 提取基因组物理定位与共线性参数 (IPG)
    ipg_df = fetch_ipg_metadata(blast_hits_df, batch_size=50)

    if ipg_df.empty:
        print("[x] 未能匹配到基因组定位记录。")
        sys.exit(1)

    # 4. 关联本地 ETE3 谱系信息
    if ETE3_AVAILABLE and "taxid" in ipg_df.columns:
        ete3_df = get_ete3_lineage(ipg_df["taxid"])
        if not ete3_df.empty:
            ipg_df = pd.merge(ipg_df, ete3_df, on="taxid", how="left")

    # 5. 排序并导出最终 CSV
    ipg_df.sort_values(by=["assembly_accession", "replicon_contig", "cds_start"], inplace=True)
    ipg_df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")
    
    print(f"\n[✓] {PROTEIN_NAME} 挖掘分析全部完成！")
    print(f"[*] 结果表格保存至: {OUTPUT_CSV}")
    print("\n前 3 行预览：")
    preview_cols = [c for c in ["protein_accession", "pident", "evalue", "assembly_accession", 
                               "replicon_contig", "cds_start", "cds_end", "strand", "phylum", "tree_leaf_label"] if c in ipg_df.columns]
    print(ipg_df[preview_cols].head(3))