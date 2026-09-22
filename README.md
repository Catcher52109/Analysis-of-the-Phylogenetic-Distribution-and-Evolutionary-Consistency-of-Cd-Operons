# Code Description

This repository contains scripts for identifying homologous proteins, reconstructing the phylogenetic distribution of the Cad operon, analyzing operon synteny, and comparing the topological consistency between operon gene trees and the species phylogeny.

## 1. PSI-BLAST

The `PSI-Blast` directory contains scripts for homologous protein mining using PSI-BLAST against the NCBI non-redundant protein database and NCBI taxonomy resources.

### Scripts

- `CadJ.py`
- `PepSY.py`
- `CadG.py`
- `CadS.py`

These scripts are used to search for homologous proteins of CadJ, PepSY, CadG, and CadS, respectively, using PSI-BLAST.

For each homologous protein, the scripts retrieve and organize information including:

- Taxonomy ID (`TaxID`)
- Protein accession
- Genomic location
- Other associated sequence and taxonomy information

These outputs are subsequently used for operon identification, synteny analysis, and phylogenetic reconstruction.

---

## 2. Phylogenetic Tree Construction

The `Phylogenetic Tree Construction` directory contains scripts for identifying homologous operons, analyzing their genomic organization, reconstructing strain-level phylogenetic relationships, and generating visualization files.

### `operon_synteny_cluster_and_profiler.py`

This script performs physical clustering and synteny analysis of homologous Cad operons.

Its main functions include:

- Identifying physically clustered Cad homologs
- Detecting strains carrying similar operon architectures
- Characterizing the gene composition of candidate operons
- Assessing local gene order and microsynteny

### `build_strain_taxonomy_tree.py`

This script constructs a strain-level phylogenetic tree based on the NCBI Taxonomy database.

It is used to organize candidate strains according to their taxonomic relationships and to generate a phylogenetic framework for downstream operon distribution analysis.

### `plot_tree_operon_synteny_alignment.py`

This script performs the integrated visualization of:

- The phylogenetic tree
- Operon microsynteny
- Gene presence/absence patterns
- Sequence identity information
- iTOL-compatible annotation tracks

The generated files can be imported into [iTOL](https://itol.embl.de/) for visualization of the phylogenetic distribution and structural variation of the Cad operon.

---

## 3. Comparison of Topological Consistency

The `Comparison of Topological Consistency` directory contains scripts for selecting target strains, retrieving housekeeping genes, reconstructing species and gene trees, and quantitatively comparing phylogenetic tree topologies.

### `filter_target_operon_strains.py`

This script is used to identify strains carrying the target operon configurations and to prepare their associated taxonomy information.

The selected strains are used for subsequent species-tree and gene-tree reconstruction.

### [`fetch_housekeeping_with_seqs.py`](https://github.com/Catcher52109/Analysis-of-the-Phylogenetic-Distribution-and-Evolutionary-Consistency-of-Cd-Operons/blob/main/Comparison%20of%20Topological%20Consistency/fetch_housekeeping_with_seqs.py)

This script retrieves housekeeping genes from NCBI in batch mode.

The main outputs include:

- Protein accession numbers
- Full-length protein sequences
- Housekeeping-gene sequence datasets
- Sequence information required for species-tree supermatrix construction

### `Script for Generating a Summary File of Topological Consistency Comparisons.py`

This script summarizes the protein accession information for:

- Cad operon genes
- Housekeeping genes

across all selected candidate strains.

The resulting summary table is used to organize the sequence datasets required for phylogenetic tree reconstruction and topological consistency comparisons.

### [`step1_extract_fastas.py`](https://github.com/Catcher52109/Analysis-of-the-Phylogenetic-Distribution-and-Evolutionary-Consistency-of-Cd-Operons/blob/main/Comparison%20of%20Topological%20Consistency/step1_extract_fastas.py)

This script extracts the protein sequences for each target gene and generates individual FASTA files for downstream sequence alignment and phylogenetic analysis.

### [`step2_rebuild_trees.sh`](https://github.com/Catcher52109/Analysis-of-the-Phylogenetic-Distribution-and-Evolutionary-Consistency-of-Cd-Operons/blob/main/Comparison%20of%20Topological%20Consistency/step2_rebuild_trees.sh)

This shell script reconstructs the phylogenetic trees used in the topological consistency analysis.

It includes:

- Construction of the concatenated RpoB–AtpD species reference tree
- Construction of the GyrB gene tree
- Construction of individual CadS, PepSY, CadG, and CadJ gene trees

These phylogenies are subsequently used for quantitative topology comparisons.

### [`step3_eval_topology_rpoB_atpD.py`](https://github.com/Catcher52109/Analysis-of-the-Phylogenetic-Distribution-and-Evolutionary-Consistency-of-Cd-Operons/blob/main/Comparison%20of%20Topological%20Consistency/step3_eval_topology_rpoB_atpD.py)

This script evaluates the topological consistency among the reconstructed phylogenetic trees.

Two types of comparisons are performed:

1. **Gene tree vs. species tree**
   - GyrB vs. RpoB–AtpD species tree
   - CadS vs. RpoB–AtpD species tree
   - PepSY vs. RpoB–AtpD species tree
   - CadG vs. RpoB–AtpD species tree
   - CadJ vs. RpoB–AtpD species tree

2. **Pairwise comparisons among Cad operon gene trees**
   - CadS vs. PepSY
   - CadS vs. CadG
   - CadS vs. CadJ
   - PepSY vs. CadG
   - PepSY vs. CadJ
   - CadG vs. CadJ

The script quantifies phylogenetic topological similarity using metrics such as:

- Normalized Robinson–Foulds distance
- Quartet similarity score

These analyses are used to evaluate the evolutionary consistency among Cad operon genes and to compare their evolutionary patterns with the strain-level species phylogeny.

---

## Workflow Overview

```text
PSI-BLAST homolog search
        |
        v
Homologous protein and taxonomy retrieval
        |
        v
Operon physical clustering and synteny analysis
        |
        v
Identification of strains carrying Cad operons
        |
        v
Strain-level taxonomy-based phylogenetic reconstruction
        |
        v
Housekeeping gene retrieval
        |
        v
FASTA extraction and sequence alignment
        |
        v
RpoB–AtpD species-tree reconstruction
        |
        +-----------------------------+
        |                             |
        v                             v
Cad gene-tree reconstruction      GyrB tree reconstruction
        |                             |
        +-------------+---------------+
                      |
                      v
          Topological consistency analysis
                      |
                      v
   Normalized RF distance and Quartet similarity
> **Note:** The code in this repository was developed with assistance from Gemini.
