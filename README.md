# incidence_catalog

Data files and R scripts for "A catalog of new incidence estimates of monogenic neurodevelopmental disorders"

## Introduction
A large fraction of rare and severe neurodevelopmental disorders (NDDs) are caused by sporadic de novo variants. Epidemiological disease estimates are not available for the vast majority of these de novo, monogenic NDDs due to phenotypic heterogeneity and the absence of large-scale genomic screens. Yet, knowledge of disease incidence is important for clinicians and researchers to guide health policy planning. Here, we adjusted a statistical method based on genetic data to predict, for the first time, the incidences of 101 known de novo variant-associated neurodevelopmental disorders as well as 3,106 putative monogenic disorders. 

## Contents

- **incidence_catalog.R**: R script used to conduct analysis.
- **constraint_gnomAD.csv** : Data extracted from gnomAD.
- **DNV_NDD_genes_from_papers.csv** : list of 101 de novo variant associated neurodevelopmental disorder genes collected from Homsy et al., 2015; Deciphering Developmental Disorders Study, 2017; Furey et al., 2018; and Heyne et al., 2018.
- **ExAC_pmuts_fs_canonical_transcripts_div_Feb07.txt** : de novo mutational model from Samocha et al., 2014.
- **Homo sapiens_consolidated.csv** : data regarding human gene essentiality from the Online GEne Essentiality database (extracted October, 2018).
- **missense_collapsed.tsv** : total per gene number of unique pathogenic or likely pathogenic missense variants reported in ClinVar and HGMD.
- **ptv_collapsed.tsv** : total per gene number of unique pathogenic or likely pathogenic protein-truncating variants reported in ClinVar and HGMD.
- **gene_panel_data_lindy.csv** : diagnostic gene panel results from Lindy et al., 2018.
- **gene_panel_data_truty.csv** : diagnostic gene panel results from Truty er al., N.D.
- **list_mouse_het_lethal_genes.tsv** : list of genes with heterozygous embryonically lethal orthologs in mice.
- **CEGv2_subset_universe.tsv** : List of genes essential for human cells based on CRISPR screens.

## Usage

1) Place all files in the same folder.
2) On R studio, run incidence_catalog.R

incidence_catalog.R will perform the following:

**Section 1 : Incidence Estimation**
Combines gnomAD constraint metric data from constraint_gnomAD.csv and de novo mutation probabilities from ExAC_pmuts_fs_canonical_transcripts_div_Feb07.txt (Samoch et al. 2014) to estimate per gene monogenic disorder incidences for the 101 neurodevelopmental disorder genes in DNV_NDD_genes_from_papers.csv as well as 3,106 variant constrained genes from gnomAD. 

**Section 2 : Assessment of gene essentiality**
Uses human gene essentiality data from `Homo sapiens_consolidated.csv` (from Online GEne Essentiality database; extracted October, 2018) to perform a fisher's test to asses enrichment for essential or conditionally essential genes in the 3,106 variant intolerant genes extracted from gnomAD for which incidence was predicted.

**Section 3 : Disease association and correlation analysis**
Uses data on the per-gene number of missense (missense_collapsed.tsv) and protein truncating (ptv_collapsed.tsv) patient pathogenic variants reported on ClinVar and HGMD to perform a fisher's test to asses enrichment for disease-associated variants in the 3,106 variant intolerant genes for which incidence was predicted. Also uses a Kendall test to evaluate any correlation between the number of reported patient variants and the predicted incidence estimates for the 101 NDD genes, 3,106 variant intolerant genes, and for all 3,207 genes combined.

**Section 4 : Assessing the role of structural variants**
Uses diagnostc gene panel testing results (gene_panel_data_lindy.csv and gene_panel_data_truty.csv) to evaluate the proportion of pathogenic patient variants composed by copy number variants for 28 epilepsy genes.

**Section 5 : Figure 1**
Generates figure 1 from "A catalog of new incidence estimates of monogenic neurodevelopmental disorders". 

**Section 6 : Mouse lethality, human cell essentiality, and make supplementary files**
Incorporates data regarding heterozygous embryonic lethality of orthologous genes in mice (list_mouse_het_lethal_genes.tsv) and gene essentiality for human cell viability from CRISPR screens (CEGv2_subset_universe.tsv) and generates the supplementary data files with the full catalog of incidence estimates for all 3,207 genes for which incidence was predicted.
