#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本名称: mega_protein.py
功能: 提取 1594 株拥有完整三看家蛋白的菌株，并合并操纵子基因 (CadJ, PepSY, CadG, CadS)。
"""

import os
import pandas as pd

def find_existing_file(candidates):
    for fn in candidates:
        if os.path.exists(fn):
            return fn
    return None

def extract_target_protein(df, default_col_name):
    candidates = [
        col for col in df.columns 
        if any(k in col.lower() for k in ['prot', 'acc', 'target', 'anchor', 'query']) 
        and not any(k in col.lower() for k in ['assembly', 'tax', 'leaf', 'tree'])
    ]
    if candidates:
        return candidates[0]
    return default_col_name if default_col_name in df.columns else df.columns[-1]

def build_combined_phylogeny_table():
    housekeeping_file = "strains_with_housekeeping_and_16s.csv"
    fasta_1594_file = "concatenated_RpoB_GyrB_AtpD.fasta"
    output_file = "strains_1594_phylogeny_and_cad_operon.csv"

    # 自动同时兼容原名与带 _2 后缀的文件名
    operon_file_candidates = {
        'CadJ_protein': ['CadJ_genomic_synteny_ete3.csv', 'CadJ_genomic_synteny_ete3_2.csv'],
        'PepSY_protein': ['PepSY_genomic_synteny_ete3.csv', 'PepSY_genomic_synteny_ete3_2.csv'],
        'CadG_protein': ['CadG_genomic_synteny_ete3.csv', 'CadG_genomic_synteny_ete3_2.csv'],
        'CadS_protein': ['CadS_genomic_synteny_ete3.csv', 'CadS_genomic_synteny_ete3_2.csv']
    }

    if not os.path.exists(housekeeping_file):
        print(f"[错误] 未找到基础文件: {housekeeping_file}")
        return

    print(f"[*] 读取基础管家基因表: {housekeeping_file}")
    df_main = pd.read_csv(housekeeping_file, dtype=str).fillna("")

    # 1. 锁定三看家蛋白完整的 1594 株菌
    if os.path.exists(fasta_1594_file):
        print(f"[*] 检测到串联超矩阵文件 {fasta_1594_file}，以此文件提取确切的 1594 个菌株 ID...")
        valid_ids = set()
        with open(fasta_1594_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    valid_ids.add(line.strip().lstrip('>').split()[0])
        
        id_col = 'tree_leaf_id' if 'tree_leaf_id' in df_main.columns else 'assembly_accession'
        df_1594 = df_main[df_main[id_col].astype(str).isin(valid_ids)].copy()
    else:
        print("[*] 未找到串联 FASTA，通过过滤 RpoB、GyrB、AtpD 非空值锁定菌株...")
        cond = (
            (df_main['RpoB_protein'] != "") &
            (df_main['GyrB_protein'] != "") &
            (df_main['AtpD_protein'] != "")
        )
        df_1594 = df_main[cond].copy()

    print(f"[*] 已确定目标菌株数: {len(df_1594)} (预期 1594)")

    join_key = 'assembly_accession'
    df_1594[join_key] = df_1594[join_key].astype(str).str.strip()

    # 2. 逐个合并 Cad 操纵子/同线性文件
    for target_col, candidates in operon_file_candidates.items():
        filepath = find_existing_file(candidates)
        if not filepath:
            print(f"[!] 警告: 未找到文件 {candidates[0]}，该列将填充为空")
            df_1594[target_col] = ""
            continue

        print(f"[*] 正在关联同线性文件: {filepath} -> {target_col}")
        df_sub = pd.read_csv(filepath, dtype=str).fillna("")

        sub_join_key = 'assembly_accession' if 'assembly_accession' in df_sub.columns else df_sub.columns[0]
        df_sub[sub_join_key] = df_sub[sub_join_key].astype(str).str.strip()

        prot_col = extract_target_protein(df_sub, target_col)
        print(f"    - 使用键列: [{sub_join_key}], 提取蛋白列: [{prot_col}]")

        # 按 bitscore 降序排序去重
        if 'bitscore' in df_sub.columns:
            df_sub['bitscore_num'] = pd.to_numeric(df_sub['bitscore'], errors='coerce').fillna(0)
            df_sub = df_sub.sort_values(by='bitscore_num', ascending=False)

        df_sub_clean = (
            df_sub[[sub_join_key, prot_col]]
            [df_sub[prot_col] != ""]
            .drop_duplicates(subset=[sub_join_key], keep='first')
        )
        df_sub_clean = df_sub_clean.rename(columns={
            sub_join_key: join_key,
            prot_col: target_col
        })

        df_1594 = pd.merge(df_1594, df_sub_clean[[join_key, target_col]], on=join_key, how='left')

    final_columns = [
        'tree_leaf_id',
        'TaxID',
        'assembly_accession',
        'operon_type',
        'organism',
        'RpoB_protein',
        'GyrB_protein',
        'AtpD_protein',
        'CadJ_protein',
        'PepSY_protein',
        'CadG_protein',
        'CadS_protein'
    ]

    for col in final_columns:
        if col not in df_1594.columns:
            df_1594[col] = ""

    df_result = df_1594[final_columns].fillna("")
    df_result.to_csv(output_file, index=False, encoding='utf-8-sig')

    print(f"\n[✓] 处理完成！输出文件: {output_file}")
    print(f"[✓] 总行数: {len(df_result)}")
    print("[*] 各功能基因在 1594 株菌中的检出统计:")
    for gene in ['CadJ_protein', 'PepSY_protein', 'CadG_protein', 'CadS_protein']:
        detected = (df_result[gene] != "").sum()
        print(f"    - {gene}: {detected} 株携带 ({detected / len(df_result) * 100:.2f}%)")

if __name__ == "__main__":
    build_combined_phylogeny_table()