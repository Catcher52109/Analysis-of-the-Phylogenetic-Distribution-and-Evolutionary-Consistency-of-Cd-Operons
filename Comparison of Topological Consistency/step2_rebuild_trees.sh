#!/usr/bin/env bash
set -e

echo "=== [1/2] 正在串联 RpoB 与 AtpD 对齐序列 ==="
python3 - << 'EOF'
def read_fasta(file_path):
    seqs = {}
    cur_id = None
    cur_seq = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if cur_id:
                    seqs[cur_id] = "".join(cur_seq)
                cur_id = line
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id:
            seqs[cur_id] = "".join(cur_seq)
    return seqs

rpo_dict = read_fasta("RpoB_aligned.fasta")
atp_dict = read_fasta("AtpD_aligned.fasta")

out_file = "species_RpoB_AtpD_aligned.fasta"
count = 0
with open(out_file, "w") as f:
    for header in rpo_dict:
        if header in atp_dict:
            f.write(f"{header}\n{rpo_dict[header]}{atp_dict[header]}\n")
            count += 1
print(f" -> 成功串联 RpoB + AtpD 超矩阵: {out_file} (包含 {count} 株菌)")
EOF

echo "=== [2/2] 使用 FastTree 构建物种树与 GyrB 对照树 ==="
echo " -> 正在构建物种参考树 (RpoB + AtpD)..."
FastTree species_RpoB_AtpD_aligned.fasta > species_tree_RpoB_AtpD.nwk

echo " -> 正在构建 GyrB 对照组单基因树..."
FastTree GyrB_aligned.fasta > GyrB_tree.nwk

echo "=== 建树完成！==="
ls -lh species_tree_RpoB_AtpD.nwk GyrB_tree.nwk