###################################################################################################################

### Packages required ###

###################################################################################################################

require(readr)
require(dplyr)
require(ggplot2)
require(gdata)
require(ggpubr)
require(gridExtra)
require(annotables)
require(stringr)

###################################################################################################################

### FUNCTIONS ###

###################################################################################################################

#Assign a mutation intolerance type based on o/e and mis_Z. Truncating is only
#counted for genes with lof o/e <= 0.35, missense only counted for genes with mis_z >= 3.09
#requires vectors with pLI, mis_z
#returns a character vector of same length as inputs with gene types
assign_type <- function(lof, mis_z) {
  output <- vector("character", length(lof))
  for(n in seq_along(lof)) {
    if(lof[n] <= 0.35 & mis_z[n] >= 3.09) {
      output[n] <- "truncating + missense"
    }
    else if(lof[n] <= 0.35) {
      output[n] <- "truncating"
    }
    else if(mis_z[n] >= 3.09) {
      output[n] <- "missense"
    }
    else {
      output[n] <- "NA"
    }
  }
  factor(output)
}

#Assign a mutation intolerance type for known DNV-associated genes based on reported variants
#genes with only PTVs are assigned truncating, those with only missense variants are assigned missense
#requires vectors with counts for missense and PTVs
#returns a character vector of same length as inputs with gene types
known_assign_type <- function(missense, PTV) {
  output <- vector("character", length(PTV))
  for(n in seq_along(PTV)) {
    if(PTV[n] == 0) {
      output[n] <- "missense"
    }
    else if(missense[n] == 0) {
      output[n] <- "truncating"
    }
    else if(missense[n] > 0 & PTV[n] > 0) {
      output[n] <- "truncating + missense"
    }
    else {
      output[n] <- "NA"
    }
  }
  factor(output)
}

#fixes format issue from known DNV genes input file where variant counts are formatted as "number (number)",
#adds the two numerical values together
convert_counts <- function(variants) {
  output <- vector("double", length(variants))
  for(n in seq_along(variants)) {
    num1 <- as.double(str_split_fixed(variants[n], " ", 2)[1])
    num2 <- as.double(str_split_fixed(str_split_fixed(variants[n], " \\(", 2)[2], "\\)", 2)[1])
    if(is.na(num2)) {
      output[n] <- num1
    }
    else {
      output[n] <- num1 + num2
    }
  }
  output
}



#Calculates total number of de novo variants (PTV or missense) per 100k individuals, adjusted for gnomAD o/e's
#requires input vectors with de novo probabilities from Samocha et al., 2014; Chromosome, 
#and gnomAD o/e (can be LoF or missense)
#same function is used to calculate upper and lower bounds
#returns a double vector with number of variants/100k for that variant type
calc_variants <- function(prob, chr, oe) {
  output <- vector("double", length(prob))
  for(n in seq_along(prob)) {
    if(oe[n] > 1) {
      output[n] <- 0    
    }
    else {
      if(chr[n] == "Y") {
        output[n] <- prob[n]*50000*(1-oe[n])
      }
      else if(chr[n] == "X") {
        output[n] <- prob[n]*50000*(1-oe[n]) + prob[n]*2*50000*(1-oe[n])
      }
      else if(chr[n] == "NA") {
        output[n] <- 0
      }
      else {
        output[n] <- prob[n]*2*100000*(1-oe[n])
      }
    }
  }
  output
}

#calculate total incidence per 100k based on gene type
#requires vectors with gene_type, trunc variants per 100k and miss variants per 100k all of same length
#returns a numerical vector of same length as inputs with predicted incidence values
calc_incidence <- function(gene_type, trunc_var, miss_var) { 
  output <- vector("double", length(gene_type))
  for(n in seq_along(gene_type)) {
    if(gene_type[n] == "truncating + missense") {
      output[n] <- trunc_var[n] + miss_var[n]
    }
    else if(gene_type[n] == "truncating") {
      output[n] <- trunc_var[n]
    }
    else if(gene_type[n] == "missense") {
      output[n] <- miss_var[n]
    }
    else {
      output[n] <- 0
    }
  }
  output
}

###################################################################################################################

### Data Input files ###

###################################################################################################################

#Make sure all the required files are in the working directory!

#Extracts the constraint data from gnomAD
gnomad_data <- read_tsv("constraint_gnomAD.csv")

#Extracts probability of de novo mutation from Samocha et al., 2014
#extracts only the necessary variables: transcript ID, gene name, chromosome (as characters)
#transcript length, probability of missense and probability of loss-of-function (as doubles)
dnv_prob_data <- read_tsv("ExAC_pmuts_fs_canonical_transcripts_div_Feb07.txt", col_types = "c_c__d___d___d") %>%
  mutate(transcript = str_replace(transcript,
                                  pattern = ".[0-9]+$",
                                  replacement = ""))

#Extracts list of known de novo genes compiled from literature as well as the variant counts from the papers
#and assigns gene type based on reported variants
known_de_novo <- read_csv("DNV_NDD_genes_from_papers.csv") %>%
  filter(!is.na(Gene)) %>%
  mutate(mis_c = convert_counts(Missense)) %>%
  mutate(PTV_c = convert_counts(PTV)) %>%
  mutate(gene_type = known_assign_type(mis_c, PTV_c))

#Load human gene essentiality data from OGEE
hgenes <- read_csv("Homo Sapiens_consolidated.csv") %>%
  select(1, 5:6) #keep only gene ID and essentiality
#Rename columns to something easier to work with
colnames(hgenes)[2] <- "ess_stat"
colnames(hgenes)[3] <- "essentiality"

#Clinvar/HGMD data files
clinvar_mis <- read_tsv("missense_collapsed.tsv", col_names = FALSE, col_types = "cccccccccc")
clinvar_ptv <- read_tsv("ptv_collapsed.tsv", col_names = FALSE, col_types = "cccccccccc")

#Diagnostic gene panel results data files
gpanel_data <- read_csv("gene_panel_data_lindy.csv")
gpanel_data_truty <- read_csv("gene_panel_data_truty.csv")

#Data on mouse lethal genes
mouse_het_lethal <- read_tsv("list_mouse_het_lethal_genes.tsv", col_names = FALSE)%>%
  mutate(het_lethal = str_to_upper(X1))

#Data for essential genes based on CRISPR screens
cell_essential <- read_tsv("CEGv2_subset_universe.tsv", col_names = FALSE)


###################################################################################################################

### Incidence Estimation ###

###################################################################################################################

#Combine mutation probability data with gnomAD data
gnomad_probs <- gnomad_data %>%
  filter(canonical == TRUE) %>%
  left_join(dnv_prob_data, by = "transcript") 

#Make df with gnomAD data for the known de novo genes
gnomad_known <- gnomad_probs %>%
  semi_join(known_de_novo, by = c("gene" = "Gene")) %>% 
  left_join(known_de_novo[,-(2:4)], by = c("gene" = "Gene")) %>%
  mutate(ptv_100k = calc_variants(p_lof,chr, oe_lof), ptv_100k_lo = calc_variants(p_lof,chr, oe_lof_upper), ptv_100k_hi = calc_variants(p_lof,chr, oe_lof_lower)) %>% 
  mutate(mv_100k = calc_variants(p_mis, chr, oe_mis), mv_100k_lo = calc_variants(p_mis, chr, oe_mis_upper), mv_100k_hi = calc_variants(p_mis, chr, oe_mis_lower)) %>%
  mutate(incid = calc_incidence(gene_type, ptv_100k, mv_100k), incid_lo = calc_incidence(gene_type, ptv_100k_lo, mv_100k_lo), incid_hi = calc_incidence(gene_type, ptv_100k_hi, mv_100k_hi))

#Filter gnomAD data to identify variant intolerant genes based on oe_lof and missense Zscore
gnomad_filt <- gnomad_probs %>%
  filter(oe_lof_upper <= 0.35 | mis_z >= 3.09) %>%
  anti_join(known_de_novo, by = c("gene" = "Gene")) #gets rid of any known de novo genes in the list

#replaces any NAs (generated when combining the probabilities to the gnomAD data) present in VI genes for which there
#is no mutation probabilities available
gnomad_filt[is.na(gnomad_filt)] <- 0

#Preliminary dataset which still includes VI genes with gnomAD issues
gnomad_filt_w_issues <- gnomad_filt %>%
  #assign type of mutation intolerance to each gene based on oe_lof/mis_z
  mutate(gene_type = assign_type(oe_lof_upper, mis_z)) %>%
  #expected variant per 100k for each type and prevalence based on gene_type, with 90% CIs:
  mutate(ptv_100k = calc_variants(p_lof,chr, oe_lof), ptv_100k_lo = calc_variants(p_lof,chr, oe_lof_upper), ptv_100k_hi = calc_variants(p_lof,chr, oe_lof_lower)) %>% 
  mutate(mv_100k = calc_variants(p_mis, chr, oe_mis), mv_100k_lo = calc_variants(p_mis, chr, oe_mis_upper), mv_100k_hi = calc_variants(p_mis, chr, oe_mis_lower)) %>%
  mutate(incid = calc_incidence(gene_type, ptv_100k, mv_100k), incid_lo = calc_incidence(gene_type, ptv_100k_lo, mv_100k_lo), incid_hi = calc_incidence(gene_type, ptv_100k_hi, mv_100k_hi))

#Eliminates genes with gnomAD issues (except synonym outliers)
gnomad_filt2 <- gnomad_filt_w_issues %>% filter(gene_issues == "[]" | gene_issues == '["syn_outlier"]')

#For incidence statistics, gets rid of VI genes with no probabilities (i.e., unable to calculate incidence for these)
gnomad_filt3 <- gnomad_filt2 %>% filter(incid > 0)

#makes a data frame with the combined total estimated incidence from all known/VI genes
tot_inc <- data_frame( type = c("Missense", "Truncation", "Truncation and Missense"),
                       VI_incid = c(sum(gnomad_filt3$incid[gnomad_filt3$gene_type == "missense"]),
                                    sum(gnomad_filt3$incid[gnomad_filt3$gene_type == "truncating"]),
                                    sum(gnomad_filt3$incid[gnomad_filt3$gene_type == "truncating + missense"])),
                       VI_incid_lo = c(sum(gnomad_filt3$incid_lo[gnomad_filt3$gene_type == "missense"]),
                                       sum(gnomad_filt3$incid_lo[gnomad_filt3$gene_type == "truncating"]),
                                       sum(gnomad_filt3$incid_lo[gnomad_filt3$gene_type == "truncating + missense"])),
                       VI_incid_hi = c(sum(gnomad_filt3$incid_hi[gnomad_filt3$gene_type == "missense"]),
                                       sum(gnomad_filt3$incid_hi[gnomad_filt3$gene_type == "truncating"]),
                                       sum(gnomad_filt3$incid_hi[gnomad_filt3$gene_type == "truncating + missense"])),
                       known_incid = c(sum(gnomad_known$incid[gnomad_known$gene_type == "missense"]),
                                       sum(gnomad_known$incid[gnomad_known$gene_type == "truncating"]),
                                       sum(gnomad_known$incid[gnomad_known$gene_type == "truncating + missense"])),
                       known_incid_lo = c(sum(gnomad_known$incid_lo[gnomad_known$gene_type == "missense"]),
                                          sum(gnomad_known$incid_lo[gnomad_known$gene_type == "truncating"]),
                                          sum(gnomad_known$incid_lo[gnomad_known$gene_type == "truncating + missense"])),
                       known_incid_hi = c(sum(gnomad_known$incid_hi[gnomad_known$gene_type == "missense"]),
                                          sum(gnomad_known$incid_hi[gnomad_known$gene_type == "truncating"]),
                                          sum(gnomad_known$incid_hi[gnomad_known$gene_type == "truncating + missense"]))
) 

###################################################################################################################

### Assessment of gene essentiality ###

###################################################################################################################

g_essent <- gnomad_filt2 %>%
  left_join(grch37_tx2gene, by = c("transcript" = "enstxp")) %>%
  left_join(hgenes, by = c("ensgene" = "locus")) %>%
  #replaces NA (i.e., missing) essentiality with "Unknown"
  mutate(ess_stat = replace(ess_stat, is.na(ess_stat), "Unknown")) %>%    
  mutate(essentiality = replace(essentiality, is.na(essentiality), "Unknown")) 

known_g_essent <- gnomad_known %>%
  left_join(grch37_tx2gene, by = c("transcript" = "enstxp")) %>%
  left_join(hgenes, by = c("ensgene" = "locus")) %>%
  #replaces NA (i.e., missing) essentiality with "Unknown"
  mutate(ess_stat = replace(ess_stat, is.na(ess_stat), "Unknown")) %>%    
  mutate(essentiality = replace(essentiality, is.na(essentiality), "Unknown"))

#Make a data frame with all the genes in human genome which are not VI genes
nonVI_genes <- hgenes %>%
  anti_join(g_essent, by = c("locus" = "ensgene"))

#data frame with counts for enrichment analysis
ess_enrichment <- data_frame(
  "essentiality" = c("Essential", "Conditional", "Nonessential"),
  "n_in_VI" = c(sum(g_essent$essentiality == "Essential"), sum(g_essent$essentiality == "Conditional"), sum(g_essent$essentiality == "Nonessential")),
  "r_in_VI" = rep(3104, 3) - c(sum(g_essent$essentiality == "Essential"), sum(g_essent$essentiality == "Conditional"), sum(g_essent$essentiality == "Nonessential")),
  "n_in_nonVI" = c(sum(nonVI_genes$essentiality == "Essential"), sum(nonVI_genes$essentiality == "Conditional"), sum(nonVI_genes$essentiality == "Nonessential")),
  "r_in_nonVI" = rep (length(nonVI_genes$essentiality), 3) - c(sum(nonVI_genes$essentiality == "Essential"), sum(nonVI_genes$essentiality == "Conditional"), sum(nonVI_genes$essentiality == "Nonessential")),
  "p_value" = vector("double", 3),
  "OR" = vector("double", 3),
  "OR_lo" = vector("double", 3),
  "OR_hi" = vector("double", 3)
)

#Enrichment analysis comparing VI genes to the rest of the genome

for(n in seq_along(ess_enrichment$essentiality)) {
  ess_ttest <- matrix(c(as.numeric(ess_enrichment[n,2]), as.numeric(ess_enrichment[n,4]), as.numeric(ess_enrichment[n,3]), as.numeric(ess_enrichment[n,5])), byrow = TRUE, nrow = 2)
  res <- fisher.test(ess_ttest)
  ess_enrichment$p_value[n] <- res$p.value
  ess_enrichment$OR[n] <- res$estimate[[1]]
  ess_enrichment$OR_lo[n] <- res$conf.int[1]
  ess_enrichment$OR_hi[n] <- res$conf.int[2]
}

###################################################################################################################

### Disease association and correlation analysis ###

###################################################################################################################


##get number of pathogenic variants for VI genes
clinvar_data_VI <- data_frame(
  gene = vector("character", length(gnomad_filt2$gene)),
  patho_mis_vars = vector("double", length(gnomad_filt2$gene)),
  patho_ptv_vars = vector("double", length(gnomad_filt2$gene)),
  tot_patho_vars = vector("double", length(gnomad_filt2$gene))
)

for(n in seq_along(gnomad_filt2$gene)) {
  clinvar_data_VI$gene[n] <- gnomad_filt2$gene[n]
  clinvar_data_VI$patho_mis_vars[n] <- sum(clinvar_mis$X8 == gnomad_filt2$gene[n])
  clinvar_data_VI$patho_ptv_vars[n] <- sum(clinvar_ptv$X8 == gnomad_filt2$gene[n])
  clinvar_data_VI$tot_patho_vars[n] <- clinvar_data_VI$patho_mis_vars[n] + clinvar_data_VI$patho_ptv_vars[n]
}

clinvar_data_VI2 <- bind_cols(gnomad_filt2, select(clinvar_data_VI, 2:4))

##get number of pathogenic variants for known DNV intolerant
clinvar_data_known <- data_frame(
  gene = vector("character", length(gnomad_known$gene)),
  patho_mis_vars = vector("double", length(gnomad_known$gene)),
  patho_ptv_vars = vector("double", length(gnomad_known$gene)),
  tot_patho_vars = vector("double", length(gnomad_known$gene))
)

for(n in seq_along(gnomad_known$gene)) {
  clinvar_data_known$gene[n] <- gnomad_known$gene[n]
  clinvar_data_known$patho_mis_vars[n] <- sum(clinvar_mis$X8 == gnomad_known$gene[n])
  clinvar_data_known$patho_ptv_vars[n] <- sum(clinvar_ptv$X8 == gnomad_known$gene[n])
  clinvar_data_known$tot_patho_vars[n] <- clinvar_data_known$patho_mis_vars[n] + clinvar_data_known$patho_ptv_vars[n]
}

clinvar_data_known2 <- bind_cols(gnomad_known, select(clinvar_data_known, 2:4))

#all unique genes with variants in clinvar/HGMD
all_clin_genes <- unique(c(clinvar_mis$X8, clinvar_ptv$X8))

#contingency matrix with number of VI genes/non-VI genes with/without pathogenic variants in ClinVar/HGMD
clinmatrix <- matrix(c(1014, 3276, 2092, 13220), 2, 2)

#Combine ClinVar/HGMD data for known and VI genes
clinvar_data_all <- bind_rows("DD_known" = select(clinvar_data_known2, "gene", "gene_type", "incid", "patho_mis_vars", "patho_ptv_vars", "tot_patho_vars"),
                              "VI" = select(clinvar_data_VI2, "gene", "gene_type", "incid", "patho_mis_vars", "patho_ptv_vars", "tot_patho_vars"),
                              .id = "gene_source")
#Only known/VI genes with pathogenic variants
clinvar_only <- clinvar_data_all %>%
  filter(tot_patho_vars > 0)

#Kendall correlation test for all NDD/VI genes with annotated pathogenic variants vs. incidence
correlation_all <- cor.test(clinvar_only$tot_patho_vars, clinvar_only$incid, method = "kendall")

#Kendall correlation test for only known NDD genes with annotated pathogenic variants vs. incidence
clinvar_only_known <- filter(clinvar_only, gene_source == "DD_known")
correlation_known <- cor.test(clinvar_only_known$tot_patho_vars, clinvar_only_known$incid, method = "kendall")

#Kendall correlation test for only VI genes with annotated pathogenic variants vs. incidence
clinvar_only_VI <- filter(clinvar_only, gene_source == "VI")
correlation_VI <- cor.test(clinvar_only_VI$tot_patho_vars, clinvar_only_VI$incid, method = "kendall")


###################################################################################################################

### Assessing the role of Structural Variants ###

###################################################################################################################

gpanel_data_2 <- gpanel_data %>%
  mutate(included = ifelse(Gene %in% gnomad_filt2$gene, "VI", 
                           ifelse(Gene %in% gnomad_known$gene, "known", NA)),
         tot_V = SNV + CNV,
         perc_SNV = SNV/tot_V
  )

gpanel_data_truty_2 <- gpanel_data_truty %>%
  filter(interpretation == "Pathogenic" | interpretation == "Likely Pathogenic") %>%
  mutate(included = ifelse(gene %in% gnomad_filt2$gene, "VI", 
                           ifelse(gene %in% gnomad_known$gene, "known", NA))) %>%
  filter(!is.na(included))

gpanel_data_truty_3 <- gpanel_data_truty_2 %>%  
  group_by(gene, included, inheritance) %>%
  summarise(
    SNVs = sum(heterozygous[!str_detect(`variant type`, "CNV")], homozygous[!str_detect(`variant type`, "CNV")], hemizygous[!str_detect(`variant type`, "CNV")]),
    CNVs = sum(`copy number = 2`[str_detect(`variant type`, "CNV")], `copy number = 3`[str_detect(`variant type`, "CNV")],
               `copy number = 4`[str_detect(`variant type`, "CNV")],
               `copy number = 5`[str_detect(`variant type`, "CNV")], `copy number = 6`[str_detect(`variant type`, "CNV")],
               heterozygous[str_detect(`variant type`, "CNV")], homozygous[str_detect(`variant type`, "CNV")],
               hemizygous[str_detect(`variant type`, "CNV")])
  ) %>%
  mutate(tot_V = SNVs + CNVs,
         perc_SNV = SNVs/tot_V) 

gpanel_data_sum <- gpanel_data_2 %>%
  filter(!is.na(included)) %>%
  left_join(gpanel_data_truty_3, by = c("Gene" = "gene", "included" = "included")) %>%
  mutate(tot_SNV = SNVs + SNV,
         tot_CNV = CNVs + CNV,
         tot_vars = tot_SNV + tot_CNV,
         tot_propSNV = tot_SNV/tot_vars,
         tot_prop_lCI = NA,
         tot_prop_uCI = NA) %>%
  select(1, 4, 2, 3, 5, 6, 8:17)

for(n in seq_along(gpanel_data_sum$Gene)) {
  bino <- binom.test(gpanel_data_sum$tot_SNV[n], gpanel_data_sum$tot_vars[n], conf.level = 0.95)
  gpanel_data_sum$tot_prop_lCI[n] <- bino$conf.int[[1]]
  gpanel_data_sum$tot_prop_uCI[n] <- bino$conf.int[[2]]
}

gpanel_data_sum2 <- gpanel_data_sum %>%
  mutate_at(c(14:16), round, digits = 2) %>%
  mutate_at(c(11:12), as.integer) %>%
  mutate(`Proportion of SNVs (95% CI)`= paste(tot_propSNV, " (", tot_prop_lCI, "-", tot_prop_uCI, ")", sep = "")) %>%
  select(1:2, 11:12, 17) %>%
  arrange(included)

###################################################################################################################

### Figure 1 ###

###################################################################################################################

panel_a <- ggplot(filter(bar_data, group == "DD"), aes(x = gene_type2, y = perc, label = perc_labels)) + 
  geom_bar(stat = "identity", fill = "#450053") + 
  geom_text(size = 3, position = position_dodge(width = 0.9), vjust = -0.40) +
  ylim(c(0,100)) +
  labs(tag = "A", x = "Variant Intolerance", y = "Percentage of Genes") + 
  theme_pubr(base_size = 8) +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_text(face = "bold", margin= margin(0,0,3,0)),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.position = c(0.5,0.9),
        plot.tag = element_text(size = 14, face = "bold"),
        plot.margin = unit(c(0,5,3,0), "mm")) +
  scale_x_discrete(label = c("MV", "PTV", "MV & PTV"))

panel_b <- ggplot(tot_inc, aes(x= type, y = known_incid)) + 
  geom_point(size = 2.5, color = "#450053", fill = "#450053") +
  geom_errorbar(aes(ymax = known_incid_hi, ymin = known_incid_lo), color = "#450053", width = 0.2, size = 0.6) +
  labs(tag = "B", x = "Variant Intolerance", y = "Incidence per\n100,000 births") +# title = "Estimated disease incidence due to known de novo genes") +
  guides(fill = "none") +
  theme_pubr(base_size = 8) +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_text(face = "bold", margin = margin(2,0,0,0)),
        axis.title.y = element_text(face = "bold"),
        plot.tag = element_text(size = 14, face = "bold"),
        plot.margin = unit(c(0,0,0,0), "mm")) +
  scale_x_discrete(label = c("MV", "PTV", "MV & PTV")) +
  scale_y_continuous(breaks = c(0,50,100,150,200,250))


fig1_layout <- rbind(c(1,2), c(1,2))
fig1 <- grid.arrange(panel_a, panel_b, layout_matrix = brain_layout)
ggsave(filename = "fig1_brain.jpeg", plot = fig1_brain, dpi = 500, width = 178, height = 60, units = "mm")



###################################################################################################################

### Mouse lethality, human cell essentiality, and make supplementary files ###

###################################################################################################################


NDD_supplement <- gnomad_known %>%
  select(-(3:5), -(9:19), -(21:25)) %>%
  left_join(clinvar_data_known, by = "gene") %>%
  unique() %>%
  #add column to supplement regarding mouse lethality and CRISPR essentiality
  mutate(het_lethal_mice = ifelse(gene %in% mouse_het_lethal$het_lethal, "heterozygous lethal", "nonlethal"),
         CRISPR_KO_essential = ifelse(gene %in% cell_essential$X1, "essential", "nonessential"),
         essentiality = known_g_essent$essentiality) 

write_excel_csv(NDD_supplement, "NDD_Supplement.csv")

VI_supplement <- g_essent %>%
  select(-(3:5), -(9:19), -(21:25), -(41:42), -(31)) %>%
  left_join(clinvar_data_VI, by = "gene") %>%
  unique() %>%
  #add column to supplement regarding mouse lethality and CRISPR essentiality
  mutate(het_lethal_mice = ifelse(gene %in% mouse_het_lethal$het_lethal, "heterozygous lethal", "nonlethal"),
         CRISPR_KO_essential = ifelse(gene %in% cell_essential$X1, "essential", "nonessential"))

write_excel_csv(VI_supplement, "VI_Supplement.csv")
