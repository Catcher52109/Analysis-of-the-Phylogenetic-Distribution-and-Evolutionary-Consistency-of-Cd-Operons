#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
步骤 1: 批量提取并格式化 FASTA 序列
"""

import os
import sys
import subprocess
import pandas as pd

NR_DB_PATH = "/mnt/chuand/Acr/nr/nr"
OPERON_CSV = "strains_1594_phylogeny_and_cad_operon.csv"
HOUSE_CSV = "strains_with_housekeeping_and_16s.csv"

def main():
    if not os.path.exists(OPERON_CSV):
        print(f"[错误] 未找到文件: {OPERON_CSV}")
        return

    print("[1/3] 读取 1594 株目标菌株信息...")
    df_operon = pd.read_csv(OPERON_CSV, dtype=str).fillna("")
    print(f" -> 目标菌株数量: {len(df_operon)}")

    # 1. 导出看家基因 FASTA (RpoB, GyrB, AtpD)
    if os.path.exists(HOUSE_CSV):
        print("[2/3] 从看家基因表中提取 RpoB, GyrB, AtpD 氨基酸序列...")
        df_house = pd.read_csv(HOUSE_CSV, dtype=str).fillna("")
        df_house_map = df_house.set_index('tree_leaf_id')
        
        for gene in ['RpoB', 'GyrB', 'AtpD']:
            fa_name = f"{gene}_proteins.fasta"
            seq_col = f"{gene}_sequence"
            written = 0
            with open(fa_name, 'w', encoding='utf-8') as f_out:
                for _, row in df_operon.iterrows():
                    leaf_id = row['tree_leaf_id']
                    if leaf_id in df_house_map.index:
                        seq = df_house_map.loc[leaf_id, seq_col]
                        if isinstance(seq, pd.Series):
                            seq = seq.iloc[0]
                        seq = str(seq).strip()
                        if seq and seq != "nan":
                            f_out.write(f">{leaf_id}\n{seq}\n")
                            written += 1
            print(f"  -> 生成 {fa_name}: 写入 {written} 条序列")

    # 2. 从本地 NR 提取 CadS, PepSY, CadG, CadJ 序列
    print("[3/3] 正在从本地 NR 数据库批量提取 Cad 操纵子基因序列...")
    cad_genes = ['CadS', 'PepSY', 'CadG', 'CadJ']
    all_cad_accs = set()
    for g in cad_genes:
        col = f"{g}_protein"
        accs = df_operon[col].replace("", None).dropna().tolist()
        all_cad_accs.update(accs)

    acc_list_file = "cad_operon_accessions.txt"
    with open(acc_list_file, 'w', encoding='utf-8') as f:
        for acc in sorted(all_cad_accs):
            f.write(f"{acc}\n")
    print(f"  -> 收集到 {len(all_cad_accs)} 个唯一 Cad 操纵子蛋白登录号，已保存至 {acc_list_file}")

    dump_fasta = "cad_operon_dump.fasta"
    cmd = [
        "blastdbcmd",
        "-db", NR_DB_PATH,
        "-entry_batch", acc_list_file,
        "-out", dump_fasta
    ]
    print(f"  -> 执行本地提取命令: {' '.join(cmd)}")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"  [!] blastdbcmd 提示:\n{res.stderr}")
    else:
        print(f"  [+] 成功导出序列文件: {dump_fasta}")

    # 解析提取出的 FASTA 文件并建立字典
    seq_dict = {}
    if os.path.exists(dump_fasta):
        cur_acc = None
        cur_seq = []
        with open(dump_fasta, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if cur_acc:
                        seq_dict[cur_acc] = "".join(cur_seq)
                    header = line[1:].split()[0]
                    cur_acc = header.split(".")[0]
                    seq_dict[header] = ""
                    cur_seq = []
                else:
                    cur_seq.append(line)
            if cur_acc:
                seq_dict[cur_acc] = "".join(cur_seq)

    # 分别导出 CadS, PepSY, CadG, CadJ 独立 FASTA
    for g in cad_genes:
        col = f"{g}_protein"
        fa_name = f"{g}_proteins.fasta"
        written = 0
        with open(fa_name, 'w', encoding='utf-8') as f_out:
            for _, row in df_operon.iterrows():
                leaf_id = row['tree_leaf_id']
                acc = row[col].strip()
                if not acc:
                    continue
                seq = seq_dict.get(acc) or seq_dict.get(acc.split(".")[0], "")
                if seq:
                    f_out.write(f">{leaf_id}\n{seq}\n")
                    written += 1
        print(f"  -> 生成 {fa_name}: 写入 {written} 条序列")

    print("\n[+] 步骤 1 执行完成！所有建树所需的标准 FASTA 文件均已就绪。")

if __name__ == "__main__":
    main()