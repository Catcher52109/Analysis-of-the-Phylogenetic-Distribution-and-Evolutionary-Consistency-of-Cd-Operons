#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本名称: step3_eval_topology_rpoB_atpD.py
功能:
  - 物种树: RpoB + AtpD 串联树
  - 对照组: GyrB 单基因树 vs 物种树
  - 实验组 1: CadS, PepSY, CadG, CadJ vs 物种树 (HGT 检验)
  - 实验组 2: 操纵子内部基因两两对比 (模块协同演化检验)
  - 输出英文标签图，避免 Linux 环境方块乱码
"""

import os
import sys
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 兼容 Python 3.11+ 的 cgi 模块
try:
    import cgi
except ImportError:
    import html, types
    fake_cgi = types.ModuleType("cgi")
    fake_cgi.escape = html.escape
    sys.modules["cgi"] = fake_cgi

from ete3 import Tree

TREE_FILES = {
    "Species": "species_tree_RpoB_AtpD.nwk",
    "GyrB":    "GyrB_tree.nwk",
    "CadS":    "CadS_tree.nwk",
    "PepSY":   "PepSY_tree.nwk",
    "CadG":    "CadG_tree.nwk",
    "CadJ":    "CadJ_tree.nwk"
}

COMPARISONS = [
    # Control group
    {"group": "Control (Baseline)", "pair_name": "GyrB vs Species Tree", "tree1": "GyrB", "tree2": "Species"},
    
    # Test Group 1: Operon genes vs Species Tree (HGT testing)
    {"group": "Test 1 (HGT Test)", "pair_name": "CadS vs Species Tree", "tree1": "CadS", "tree2": "Species"},
    {"group": "Test 1 (HGT Test)", "pair_name": "PepSY vs Species Tree", "tree1": "PepSY", "tree2": "Species"},
    {"group": "Test 1 (HGT Test)", "pair_name": "CadG vs Species Tree", "tree1": "CadG", "tree2": "Species"},
    {"group": "Test 1 (HGT Test)", "pair_name": "CadJ vs Species Tree", "tree1": "CadJ", "tree2": "Species"},
    
    # Test Group 2: Within-operon co-evolution
    {"group": "Test 2 (Module Co-evolution)", "pair_name": "CadS vs PepSY", "tree1": "CadS", "tree2": "PepSY"},
    {"group": "Test 2 (Module Co-evolution)", "pair_name": "CadS vs CadG", "tree1": "CadS", "tree2": "CadG"},
    {"group": "Test 2 (Module Co-evolution)", "pair_name": "CadS vs CadJ", "tree1": "CadS", "tree2": "CadJ"},
    {"group": "Test 2 (Module Co-evolution)", "pair_name": "PepSY vs CadG", "tree1": "PepSY", "tree2": "CadG"},
    {"group": "Test 2 (Module Co-evolution)", "pair_name": "PepSY vs CadJ", "tree1": "PepSY", "tree2": "CadJ"},
    {"group": "Test 2 (Module Co-evolution)", "pair_name": "CadG vs CadJ", "tree1": "CadG", "tree2": "CadJ"},
]

def sample_quartet_similarity(t1, t2, common_leaves, n_samples=3000):
    """基于四点拓扑距离计算 Quartet 相似度"""
    if len(common_leaves) < 4:
        return np.nan
        
    leaves_list = list(common_leaves)
    matches = 0
    resolved = 0
    
    for _ in range(n_samples):
        a, b, c, d = random.sample(leaves_list, 4)
        try:
            d1_ab = t1.get_distance(a, b, topology_only=True)
            d1_cd = t1.get_distance(c, d, topology_only=True)
            d1_ac = t1.get_distance(a, c, topology_only=True)
            d1_bd = t1.get_distance(b, d, topology_only=True)
            d1_ad = t1.get_distance(a, d, topology_only=True)
            d1_bc = t1.get_distance(b, c, topology_only=True)
            
            s1_1, s1_2, s1_3 = d1_ab + d1_cd, d1_ac + d1_bd, d1_ad + d1_bc
            top1 = "ab|cd" if (s1_1 < s1_2 and s1_1 < s1_3) else ("ac|bd" if (s1_2 < s1_1 and s1_2 < s1_3) else ("ad|bc" if (s1_3 < s1_1 and s1_3 < s1_2) else None))
            if not top1:
                continue

            d2_ab = t2.get_distance(a, b, topology_only=True)
            d2_cd = t2.get_distance(c, d, topology_only=True)
            d2_ac = t2.get_distance(a, c, topology_only=True)
            d2_bd = t2.get_distance(b, d, topology_only=True)
            d2_ad = t2.get_distance(a, d, topology_only=True)
            d2_bc = t2.get_distance(b, c, topology_only=True)
            
            s2_1, s2_2, s2_3 = d2_ab + d2_cd, d2_ac + d2_bd, d2_ad + d2_bc
            top2 = "ab|cd" if (s2_1 < s2_2 and s2_1 < s2_3) else ("ac|bd" if (s2_2 < s2_1 and s2_2 < s2_3) else ("ad|bc" if (s2_3 < s2_1 and s2_3 < s2_2) else None))
            if not top2:
                continue
                
            resolved += 1
            if top1 == top2:
                matches += 1
        except Exception:
            continue
            
    return (matches / resolved) if resolved > 0 else np.nan

def plot_publication_figure(df_res, out_png="topology_comparison_rpoB_atpD.png"):
    """绘制纯英文发表级对比图"""
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 5.5), dpi=300)
    
    color_map = {
        "Control (Baseline)": "#2ca02c",
        "Test 1 (HGT Test)": "#d62728",
        "Test 2 (Module Co-evolution)": "#1f77b4"
    }
    bar_colors = [color_map[g] for g in df_res["Group"]]
    x_pos = np.arange(len(df_res))
    
    # Subplot A: Normalized RF Distance
    ax1.bar(x_pos, df_res["Normalized_RF"], color=bar_colors, width=0.65, edgecolor='black', linewidth=0.8)
    ax1.set_title("A. Topological Distance (Normalized RF Distance)", fontsize=12.5, fontweight='bold', pad=12)
    ax1.set_ylabel("Normalized RF Distance (Lower = More Congruent)", fontsize=10.5)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(df_res["Pair"], rotation=45, ha="right", fontsize=9.5)
    ax1.set_ylim(0, 1.05)
    
    # Subplot B: Quartet Congruence Score
    q_vals = pd.to_numeric(df_res["Quartet_Score"], errors='coerce').fillna(0)
    ax2.bar(x_pos, q_vals, color=bar_colors, width=0.65, edgecolor='black', linewidth=0.8)
    ax2.set_title("B. Clade Similarity (Quartet Congruence Score)", fontsize=12.5, fontweight='bold', pad=12)
    ax2.set_ylabel("Quartet Similarity (Higher = More Congruent)", fontsize=10.5)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(df_res["Pair"], rotation=45, ha="right", fontsize=9.5)
    ax2.set_ylim(0, 1.05)
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, edgecolor='black', label=g) for g, c in color_map.items()]
    fig.legend(handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, 1.03), ncol=3, fontsize=10.5, frameon=True)
    
    plt.tight_layout()
    plt.savefig(out_png, bbox_inches='tight')
    print(f"[+] 发表级插图已保存: {out_png}")

def main():
    print("[*] 正在执行系统发育拓扑评估 (Species: RpoB+AtpD | Control: GyrB)...")
    results = []
    
    for item in COMPARISONS:
        g = item["group"]
        pair_name = item["pair_name"]
        f1 = TREE_FILES[item["tree1"]]
        f2 = TREE_FILES[item["tree2"]]
        
        if not (os.path.exists(f1) and os.path.exists(f2)):
            print(f"  [!] 跳过 {pair_name}: 缺少文件 ({f1} 或 {f2})")
            continue
            
        print(f"  [*] 评估中: [{g}] {pair_name}...")
        t1 = Tree(f1)
        t2 = Tree(f2)
        
        leaves1 = set(t1.get_leaf_names())
        leaves2 = set(t2.get_leaf_names())
        common = leaves1.intersection(leaves2)
        
        if len(common) < 4:
            continue
            
        t1_p = t1.copy()
        t2_p = t2.copy()
        t1_p.prune(common, preserve_branch_length=True)
        t2_p.prune(common, preserve_branch_length=True)
        
        rf_data = t1_p.robinson_foulds(t2_p, unrooted_trees=True)
        rf = rf_data[0]
        max_rf = rf_data[1]
        norm_rf = rf / max_rf if max_rf > 0 else 0.0
        rf_sim = 1.0 - norm_rf
        q_sim = sample_quartet_similarity(t1_p, t2_p, common, n_samples=3000)
        
        results.append({
            "Group": g,
            "Pair": pair_name,
            "Common_Strains": len(common),
            "Absolute_RF": rf,
            "Max_RF": max_rf,
            "Normalized_RF": round(norm_rf, 4),
            "RF_Similarity": round(rf_sim, 4),
            "Quartet_Score": round(q_sim, 4) if not np.isnan(q_sim) else "N/A"
        })
        
    df_res = pd.DataFrame(results)
    out_csv = "topology_summary_rpoB_atpD.csv"
    df_res.to_csv(out_csv, index=False, encoding="utf-8-sig")
    print(f"\n[+] 评估完成！统计表格保存在: {out_csv}")
    print(df_res[["Group", "Pair", "Common_Strains", "Normalized_RF", "Quartet_Score"]].to_string(index=False))
    
    plot_publication_figure(df_res)

if __name__ == "__main__":
    main()