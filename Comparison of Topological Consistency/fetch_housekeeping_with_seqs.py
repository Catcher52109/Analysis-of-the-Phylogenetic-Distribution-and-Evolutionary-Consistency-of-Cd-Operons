#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本名称: fetch_housekeeping_with_seqs.py
功能: 
  1. 读取目标菌株清单 (TaxID, assembly_accession)。
  2. 提交 NCBI 获取 16S rRNA 核酸、RpoB、GyrB、AtpD 蛋白的 Accession 与全长序列。
  3. 16S rRNA 字段紧跟在 TaxID, assembly_accession 之后输出。
  4. 失败自动指数退避重试 10 次 (适合过夜稳定运行)。
  5. 自动导出 16S rRNA 及三大看家蛋白的标准 FASTA 和串联 FASTA，供下游建树与拓扑分析。
"""

import os
import sys
import time
import random
import json
import urllib.request
import urllib.parse
import urllib.error
import pandas as pd

# ==================== 1. 文件与参数配置 ====================
if os.path.exists("target_operon_strains_for_rpoB_2.csv"):
    INPUT_CSV = "target_operon_strains_for_rpoB_2.csv"
elif os.path.exists("target_operon_strains_for_rpoB.csv"):
    INPUT_CSV = "target_operon_strains_for_rpoB.csv"
else:
    INPUT_CSV = "target_operon_strains_for_rpoB.csv"

OUTPUT_CSV = "strains_with_housekeeping_and_16s.csv"
ACCESSION_LIST_FILE = "all_target_accessions_for_blastdbcmd.txt"

# NCBI 接口配置
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "researcher@example.com")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")

# 动态延时与过夜重试配置
MIN_SLEEP = 0.5
MAX_SLEEP = 1.0
MAX_RETRIES = 10

# 蛋白基因检索词定义 (保证全长序列质量)
PROTEIN_SPECS = {
    "RpoB": {
        "query": '(rpoB[Gene Name] OR "RNA polymerase subunit beta"[Protein Name] OR "DNA-directed RNA polymerase subunit beta"[Protein Name]) AND 800:1600[Sequence Length]',
        "backup": '(rpoB[Gene Name] OR "RNA polymerase subunit beta"[Protein Name])'
    },
    "GyrB": {
        "query": '(gyrB[Gene Name] OR "DNA gyrase subunit B"[Protein Name] OR "DNA topoisomerase (ATP-hydrolyzing) subunit B"[Protein Name]) AND 500:1100[Sequence Length]',
        "backup": '(gyrB[Gene Name] OR "DNA gyrase subunit B"[Protein Name])'
    },
    "AtpD": {
        "query": '(atpD[Gene Name] OR "ATP synthase subunit beta"[Protein Name] OR "F0F1 ATP synthase subunit beta"[Protein Name]) AND 350:700[Sequence Length]',
        "backup": '(atpD[Gene Name] OR "ATP synthase subunit beta"[Protein Name])'
    }
}

# 16S rRNA 核酸检索配置 (位于 nucleotide 库)
R16S_QUERY_TEMPLATE = 'txid{taxid}[Organism:noexp] AND (16S[Gene Name] OR "16S ribosomal RNA"[Title] OR "16S rRNA"[Title]) AND 1000:2000[Sequence Length]'
R16S_BACKUP_TEMPLATE = '"{organism}"[Organism] AND ("16S ribosomal RNA" OR "16S rRNA") AND 1000:2000[Sequence Length]'

# ==================== 2. NCBI 访问核心工具 ====================
def ncbi_request(url, params, max_retries=MAX_RETRIES):
    """带 10 次指数退避重试机制的通用 NCBI API 接口"""
    full_url = f"{url}?{urllib.parse.urlencode(params)}"
    headers = {"User-Agent": f"HousekeepingMiner/3.0 ({NCBI_EMAIL})"}
    
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(full_url, headers=headers)
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read().decode("utf-8")
                if params.get("retmode") == "json":
                    return json.loads(data)
                return data
        except urllib.error.HTTPError as e:
            wait_time = min(3 + attempt * 3, 30)
            print(f"  [!] NCBI HTTP {e.code} 异常 (重试 {attempt+1}/{max_retries})，休眠 {wait_time}s...", flush=True)
            time.sleep(wait_time)
        except Exception as e:
            wait_time = min(2 + attempt * 2, 25)
            print(f"  [!] 网络波动: {e} (重试 {attempt+1}/{max_retries})，休眠 {wait_time}s...", flush=True)
            time.sleep(wait_time)
            
    print(f"  [x] 重试 {max_retries} 次后依然失败，跳过当前项。", flush=True)
    return None

def search_ncbi(db, term):
    """在指定数据库执行 esearch 获取 UIDs"""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": db,
        "term": term,
        "retmode": "json",
        "retmax": 5,
        "sort": "relevance",
        "email": NCBI_EMAIL
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
        
    data = ncbi_request(url, params)
    if data and "esearchresult" in data:
        return data["esearchresult"].get("idlist", [])
    return []

def get_summaries(db, uids):
    """使用 esummary 获取条目基本信息"""
    if not uids:
        return {}
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": db,
        "id": ",".join(uids),
        "retmode": "json",
        "email": NCBI_EMAIL
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
        
    data = ncbi_request(url, params)
    if data and "result" in data:
        return data["result"]
    return {}

def fetch_fasta_sequence(db, accession_or_uid):
    """使用 efetch 下载完整 FASTA 序列文本"""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": db,
        "id": str(accession_or_uid),
        "rettype": "fasta",
        "retmode": "text",
        "email": NCBI_EMAIL
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
        
    fasta_text = ncbi_request(url, params)
    if fasta_text and fasta_text.startswith(">"):
        lines = fasta_text.strip().split("\n")
        seq = "".join(l.strip() for l in lines[1:])
        return seq
    return ""

def pick_best_hit(summary_result, uids, is_nucleotide=False):
    """筛选最具代表性、全长的标准登录号"""
    candidates = []
    for uid in uids:
        info = summary_result.get(uid, {})
        acc = info.get("accessionversion") or info.get("caption") or ""
        title = info.get("title", "").lower()
        if not acc:
            continue
        is_partial = "partial" in title or "fragment" in title
        if is_nucleotide:
            score = (0 if is_partial else 10) + (5 if "complete" in title else 0)
        else:
            is_refseq = acc.startswith(("WP_", "NP_", "YP_"))
            score = (0 if is_partial else 10) + (5 if is_refseq else 0)
        candidates.append((score, acc, uid))
        
    if candidates:
        candidates.sort(key=lambda x: x[0], reverse=True)
        return candidates[0][1]
    return ""

# ==================== 3. 单菌株目标检索 ====================
def query_16s_rRNA(taxid, assembly, organism):
    """在 nucleotide 数据库中检索该菌株的 16S rRNA"""
    queries = [
        R16S_QUERY_TEMPLATE.format(taxid=taxid),
        f'"{assembly}" AND ("16S ribosomal RNA" OR "16S rRNA") AND 1000:2000[Sequence Length]' if assembly and assembly != "nan" else None,
        R16S_BACKUP_TEMPLATE.format(organism=organism) if organism and organism != "nan" else None,
        f'txid{taxid}[Organism:exp] AND (16S[Gene Name] OR "16S ribosomal RNA"[Title]) AND 1000:2000[Sequence Length]'
    ]
    uids = []
    for q in queries:
        if not q:
            continue
        uids = search_ncbi("nucleotide", q)
        time.sleep(random.uniform(MIN_SLEEP, MAX_SLEEP))
        if uids:
            break
            
    if uids:
        sums = get_summaries("nucleotide", uids)
        time.sleep(random.uniform(MIN_SLEEP, MAX_SLEEP))
        best_acc = pick_best_hit(sums, uids, is_nucleotide=True)
        if best_acc:
            seq = fetch_fasta_sequence("nucleotide", best_acc)
            time.sleep(random.uniform(MIN_SLEEP, MAX_SLEEP))
            return best_acc, seq
            
    return "", ""

def query_protein(taxid, assembly, organism, gene_name):
    """在 protein 数据库中检索核心看家蛋白"""
    g_spec = PROTEIN_SPECS[gene_name]
    gene_q = g_spec["query"]
    queries = [
        f"txid{taxid}[Organism:noexp] AND {gene_q}",
        f'"{assembly}" AND {gene_q}' if assembly and assembly != "nan" else None,
        f'"{organism}"[Organism] AND {gene_q}' if organism and organism != "nan" else None,
        f"txid{taxid}[Organism:exp] AND {gene_q}",
        f"txid{taxid}[Organism] AND {g_spec['backup']}"
    ]
    uids = []
    for q in queries:
        if not q:
            continue
        uids = search_ncbi("protein", q)
        time.sleep(random.uniform(MIN_SLEEP, MAX_SLEEP))
        if uids:
            break
            
    if uids:
        sums = get_summaries("protein", uids)
        time.sleep(random.uniform(MIN_SLEEP, MAX_SLEEP))
        best_acc = pick_best_hit(sums, uids, is_nucleotide=False)
        if best_acc:
            seq = fetch_fasta_sequence("protein", best_acc)
            time.sleep(random.uniform(MIN_SLEEP, MAX_SLEEP))
            return best_acc, seq
            
    return "", ""

# ==================== 4. 自动生成建树 FASTA 与本地 NR 清单 ====================
def export_downstream_files(df_out):
    print("\n[*] 正在自动生成下游进化树构建及本地 NR 所需文件...", flush=True)
    all_accs = set()
    
    # 1. 导出 16S rRNA FASTA (供构建 16S 物种参考树)
    r16s_file = "16S_rRNA.fasta"
    c_16s = 0
    with open(r16s_file, "w", encoding="utf-8") as f_16s:
        for _, row in df_out.iterrows():
            leaf_id = f"{row['TaxID']}_{row['assembly_accession']}"
            acc = str(row.get("16S_rRNA_accession", "")).strip()
            seq = str(row.get("16S_rRNA_sequence", "")).strip()
            if acc and seq and seq != "nan":
                f_16s.write(f">{leaf_id}\n{seq}\n")
                all_accs.add(acc)
                c_16s += 1
    print(f"  -> 生成 16S 物种建树序列文件: {r16s_file} (共 {c_16s} 条序列)", flush=True)

    # 2. 导出三大看家蛋白独立 FASTA
    for gene in ["RpoB", "GyrB", "AtpD"]:
        fa_file = f"{gene}_proteins.fasta"
        c_gene = 0
        with open(fa_file, "w", encoding="utf-8") as f_fa:
            for _, row in df_out.iterrows():
                leaf_id = f"{row['TaxID']}_{row['assembly_accession']}"
                acc = str(row.get(f"{gene}_protein", "")).strip()
                seq = str(row.get(f"{gene}_sequence", "")).strip()
                if acc and seq and seq != "nan":
                    f_fa.write(f">{leaf_id}\n{seq}\n")
                    all_accs.add(acc)
                    c_gene += 1
        print(f"  -> 生成蛋白建树序列文件: {fa_file} (共 {c_gene} 条序列)", flush=True)

    # 3. 导出三蛋白串联多位点超矩阵 FASTA
    concat_file = "concatenated_RpoB_GyrB_AtpD.fasta"
    c_concat = 0
    with open(concat_file, "w", encoding="utf-8") as f_c:
        for _, row in df_out.iterrows():
            leaf_id = f"{row['TaxID']}_{row['assembly_accession']}"
            r_seq = str(row.get("RpoB_sequence", "")).strip()
            g_seq = str(row.get("GyrB_sequence", "")).strip()
            a_seq = str(row.get("AtpD_sequence", "")).strip()
            if r_seq and g_seq and a_seq and "nan" not in [r_seq, g_seq, a_seq]:
                f_c.write(f">{leaf_id}\n{r_seq + g_seq + a_seq}\n")
                c_concat += 1
    print(f"  -> 生成多基因串联超矩阵文件: {concat_file} (共 {c_concat} 株拥有完整三蛋白的菌株)", flush=True)

    # 4. 导出本地 NR 反查 Accession 清单
    with open(ACCESSION_LIST_FILE, "w", encoding="utf-8") as f_list:
        for acc in sorted(all_accs):
            f_list.write(f"{acc}\n")
    print(f"  -> 生成本地 NR 反查清单: {ACCESSION_LIST_FILE} (包含 {len(all_accs)} 个唯一编号)", flush=True)

# ==================== 5. 主执行逻辑 ====================
def main():
    if not os.path.exists(INPUT_CSV):
        print(f"[错误] 未找到输入菌株文件: {INPUT_CSV}", flush=True)
        return

    print(f"[*] 载入目标菌株清单: {INPUT_CSV}", flush=True)
    df_target = pd.read_csv(INPUT_CSV)
    total_strains = len(df_target)
    print(f"[+] 成功载入 {total_strains} 株目标菌株。", flush=True)

    # 严格按照需求将 16S rRNA 字段排在 TaxID, assembly_accession 之后
    final_columns = [
        "TaxID",
        "assembly_accession",
        "16S_rRNA_accession",
        "16S_rRNA_sequence",
        "operon_type",
        "organism",
        "genus",
        "species",
        "RpoB_protein",
        "RpoB_sequence",
        "GyrB_protein",
        "GyrB_sequence",
        "AtpD_protein",
        "AtpD_sequence",
        "tree_leaf_id"
    ]

    # 支持断点续传
    if os.path.exists(OUTPUT_CSV):
        print(f"[*] 检测到历史文件 {OUTPUT_CSV}，正在无缝加载以继续断点续传...", flush=True)
        try:
            df_out = pd.read_csv(OUTPUT_CSV, dtype=str).fillna("")
            for col in final_columns:
                if col not in df_out.columns:
                    df_out[col] = ""
        except Exception:
            df_out = df_target.copy()
            for col in final_columns:
                if col not in df_out.columns:
                    df_out[col] = ""
    else:
        df_out = df_target.copy()
        for col in final_columns:
            if col not in df_out.columns:
                df_out[col] = ""

    # 规范列顺序
    existing_cols = [c for c in final_columns if c in df_out.columns] + [c for c in df_out.columns if c not in final_columns]
    df_out = df_out[existing_cols]

    print(f"[*] 启动多基因与 16S rRNA 全量抓取流程 (重试上限: {MAX_RETRIES} 次，防封延时: {MIN_SLEEP}~{MAX_SLEEP}s)...", flush=True)

    for idx, row in df_out.iterrows():
        taxid = str(row.get("TaxID", "")).strip()
        assembly = str(row.get("assembly_accession", "")).strip()
        organism = str(row.get("organism", "")).strip()

        # 检查是否全部已抓取
        r16_acc = str(row.get("16S_rRNA_accession", "")).strip()
        r16_seq = str(row.get("16S_rRNA_sequence", "")).strip()
        r_acc = str(row.get("RpoB_protein", "")).strip()
        r_seq = str(row.get("RpoB_sequence", "")).strip()
        g_acc = str(row.get("GyrB_protein", "")).strip()
        g_seq = str(row.get("GyrB_sequence", "")).strip()
        a_acc = str(row.get("AtpD_protein", "")).strip()
        a_seq = str(row.get("AtpD_sequence", "")).strip()

        if r16_acc and r16_seq and r_acc and r_seq and g_acc and g_seq and a_acc and a_seq:
            continue

        print(f"[{idx+1}/{total_strains}] ({(idx+1)/total_strains*100:.1f}%) 菌株: TaxID {taxid} | Assembly {assembly} | {organism}", flush=True)

        # 1. 抓取 16S rRNA (核酸)
        if not (r16_acc and r16_seq):
            r16_acc, r16_seq = query_16s_rRNA(taxid, assembly, organism)
            df_out.at[idx, "16S_rRNA_accession"] = r16_acc
            df_out.at[idx, "16S_rRNA_sequence"] = r16_seq

        # 2. 抓取 RpoB (蛋白)
        if not (r_acc and r_seq):
            r_acc, r_seq = query_protein(taxid, assembly, organism, "RpoB")
            df_out.at[idx, "RpoB_protein"] = r_acc
            df_out.at[idx, "RpoB_sequence"] = r_seq

        # 3. 抓取 GyrB (蛋白)
        if not (g_acc and g_seq):
            g_acc, g_seq = query_protein(taxid, assembly, organism, "GyrB")
            df_out.at[idx, "GyrB_protein"] = g_acc
            df_out.at[idx, "GyrB_sequence"] = g_seq

        # 4. 抓取 AtpD (蛋白)
        if not (a_acc and a_seq):
            a_acc, a_seq = query_protein(taxid, assembly, organism, "AtpD")
            df_out.at[idx, "AtpD_protein"] = a_acc
            df_out.at[idx, "AtpD_sequence"] = a_seq

        print(f"  -> 16S: {r16_acc or '未命中'} ({len(r16_seq)} bp) | RpoB: {r_acc or '未命中'} | GyrB: {g_acc or '未命中'} | AtpD: {a_acc or '未命中'}", flush=True)

        # 每 5 株自动保存一次
        if (idx + 1) % 5 == 0 or (idx + 1) == total_strains:
            df_out.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")

    df_out.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")
    export_downstream_files(df_out)
    print(f"\n[+] 全部抓取任务圆满完成！最终整合汇总表格保存在: {OUTPUT_CSV}", flush=True)

if __name__ == "__main__":
    main()