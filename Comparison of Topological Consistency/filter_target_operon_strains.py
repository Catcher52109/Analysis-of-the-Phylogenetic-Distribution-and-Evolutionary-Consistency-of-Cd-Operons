import os
import re
import pandas as pd

def main():
    infile = 'itol_operon_binary_matrix_reordered.txt'
    outfile = 'target_operon_strains_for_rpoB.csv'

    if not os.path.exists(infile):
        print(f"[错误] 未在当前目录找到二进制矩阵文件: {infile}")
        return

    print("[1/3] 正在解析 iTOL 二进制矩阵文件...")
    strains_data = []
    with open(infile, 'r', encoding='utf-8') as f:
        start_data = False
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line == 'DATA':
                start_data = True
                continue
            if start_data:
                parts = line.split('\t')
                if len(parts) == 5:
                    strains_data.append(parts)

    df_matrix = pd.DataFrame(strains_data, columns=['full_id', 'CadJ', 'PepSY', 'CadG', 'CadS'])
    for col in ['CadJ', 'PepSY', 'CadG', 'CadS']:
        df_matrix[col] = df_matrix[col].astype(int)

    # 筛选规则：
    # 1. 完整4基因：1-1-1-1
    # 2. CadJ缺失型（CadS+PepSY+CadG）：0-1-1-1
    cond_4genes = (df_matrix['CadJ'] == 1) & (df_matrix['PepSY'] == 1) & (df_matrix['CadG'] == 1) & (df_matrix['CadS'] == 1)
    cond_del_cadj = (df_matrix['CadJ'] == 0) & (df_matrix['PepSY'] == 1) & (df_matrix['CadG'] == 1) & (df_matrix['CadS'] == 1)

    df_matrix['operon_type'] = 'Other'
    df_matrix.loc[cond_4genes, 'operon_type'] = 'Complete_4_Genes'
    df_matrix.loc[cond_del_cadj, 'operon_type'] = 'CadS_PepSY_CadG (CadJ-lacking)'

    selected_df = df_matrix[df_matrix['operon_type'] != 'Other'].copy()

    print(f"[2/3] 正在提取 TaxID 与 Assembly 编号 (筛选出 {len(selected_df)} 条记录)...")
    def parse_id(full_id):
        m = re.match(r'^(\d+)_(GC[AF]_\d+\.\d+)$', full_id)
        if m:
            return m.group(1), m.group(2)
        parts = full_id.split('_')
        return parts[0], '_'.join(parts[1:])

    parsed = selected_df['full_id'].apply(parse_id)
    selected_df['TaxID'] = [p[0] for p in parsed]
    selected_df['assembly_accession'] = [p[1] for p in parsed]

    print("[3/3] 正在从同目录下的 synteny CSV 中映射物种学名信息...")
    org_map = {}
    synteny_files = ['CadG_genomic_synteny_ete3.csv', 'CadS_genomic_synteny_ete3.csv', 'CadJ_genomic_synteny_ete3.csv', 'PepSY_genomic_synteny_ete3.csv']
    for sf in synteny_files:
        if os.path.exists(sf):
            try:
                sdf = pd.read_csv(sf, usecols=['assembly_accession', 'organism', 'genus', 'species']).drop_duplicates(subset=['assembly_accession'])
                for _, r in sdf.iterrows():
                    acc = r['assembly_accession']
                    if acc not in org_map and pd.notna(r['organism']):
                        org_map[acc] = (r['organism'], r.get('genus', ''), r.get('species', ''))
            except Exception as e:
                print(f"  [提示] 读取 {sf} 时跳过部分列: {e}")

    selected_df['organism'] = selected_df['assembly_accession'].apply(lambda x: org_map.get(x, ('', '', ''))[0])
    selected_df['genus'] = selected_df['assembly_accession'].apply(lambda x: org_map.get(x, ('', '', ''))[1])
    selected_df['species'] = selected_df['assembly_accession'].apply(lambda x: org_map.get(x, ('', '', ''))[2])
    selected_df['tree_leaf_id'] = selected_df['full_id']

    out_cols = ['TaxID', 'assembly_accession', 'operon_type', 'CadJ', 'PepSY', 'CadG', 'CadS', 'organism', 'genus', 'species', 'tree_leaf_id']
    final_df = selected_df[out_cols]
    final_df.to_csv(outfile, index=False, encoding='utf-8-sig')
    
    print(f"\n[成功] 汇总文件已生成: {outfile}")
    print(f" -> 包含完整4基因菌株: {(final_df['operon_type']=='Complete_4_Genes').sum()} 个")
    print(f" -> 包含CadJ缺失型菌株: {(final_df['operon_type'].str.contains('CadJ-lacking')).sum()} 个")
    print(f" -> 总计记录数: {len(final_df)}")

if __name__ == '__main__':
    main()