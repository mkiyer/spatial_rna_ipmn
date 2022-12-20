##################################################
# Nanostring GeoMx Spatial RNA Profiling Analysis
#
# Author: Matthew K. Iyer
# Copyright: 2022
##################################################
library(tidyverse)
library(readxl)
library(writexl)
library(edgeR)
library(Rtsne)
library(ggridges)
library(ggrepel)
library(eulerr)
library(pheatmap)
library(fgsea)
library(randomForest)
library(qvalue)
library(igraph)
library(msigdbr)
library(clusterProfiler)

##################################################
# Configuration
##################################################

# base directories for analysis
working_dir = "~/OneDrive - Duke University/research/projects/st_ipmn/current"
data_dir = "~/OneDrive - Duke University/research/projects/st_ipmn/data"

# spatial rna input file
input_xlsx_file = file.path(data_dir, "ipmn_counts_manuscript_2022-07-17.xlsx")
# spatial rna sheet names
metadata_sheet <- "SegmentProperties"
metadata_id_key <- "SegmentDisplayName"
count_data_sheet <- "BioProbeCountMatrix"

# external resource files
hpa_path_file <- file.path(data_dir, "hpa_pathology.tsv")
cptac_panc_file <- file.path(data_dir, "cptac_panc_supp_mmc3.xlsx")
oncotarget_panc_file <- file.path(data_dir, "oncotarget-08-42537-s002.xlsx")

# segments included
analysis_segments = c("panck")

# qc filtering parameters
qc_min_counts = 100000
qc_min_auc = 0.65
qc_negprobe_lod_quantile = 0.9
qc_min_negprobe_geomean = 1
qc_min_frac_expr = 0.2

# de params
de_padj_cutoff = 0.05
de_log2fc_cutoff = 1.0

# gsea params
gsea_padj_cutoff <- 0.05
gsea_log2fc_cutoff <- 1

# output files
processed_xlsx_file <- "processed_data.xlsx"
tsne_results_xlsx_file <- "tsne_allgenes.xlsx"
de_results_xlsx <- "de_results.xlsx"
gsea_results_xlsx <- "gsea_results.xlsx"
cta_bg_gene_symbols_txt <- "cta_gene_symbols.txt"
community_enrichment_xlsx <- "community_enrichment.xlsx"

# output directories
if (!dir.exists(working_dir)) {
  dir.create(working_dir)
}
setwd(working_dir)

plot_dir <- file.path(working_dir, "plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}
gene_plot_dir <- file.path(working_dir, "gene_plots")
if (!dir.exists(gene_plot_dir)) {
  dir.create(gene_plot_dir)
}
gsea_plot_dir <- file.path(working_dir, "gsea_plots")
if (!dir.exists(gsea_plot_dir)) {
  dir.create(gsea_plot_dir)
}
de_result_dir <- file.path(working_dir, "de_results")
if (!dir.exists(de_result_dir)) {
  dir.create(de_result_dir)
}

get_color_scales <- function(slide_ids) {
  color_scales = list(
    path = c(lgd="#AAAAAA", hgd="#DD0000"),
    pathbw = c(lgd="#cccccc", hgd="#333333"),
    histology = c(pb = "#32cd32", intestinal = "#670092", gf="#1677c4"),
    histpath = c(panck_intlgd="#f0ceff", panck_inthgd="#b042ff", 
                 panck_pblgd="#1677c4", panck_pbhgd="#00ab08"),
    carcinoma = c(no = "#666666", yes = "#ff0000"),
    prognostic = c(no="#aaaaaa", fav="#00FFFF", unfav="#FFFF00"),
    de_cptac_rna = c(no="#aaaaaa", dn="#00FFFF", up="#FFFF00", na="#FFFFFF"), 
    de_cptac_prot = c(no="#aaaaaa", dn="#00FFFF", up="#FFFF00", na="#FFFFFF"),
    grade_subtype = c(none = "#888888", hgd = "#FF0000", lgd = "#0000FF", pbhgd = "#32cd32",
                      inthgd = "#670092", pblgd = "#1677c4", intlgd = "#ff6700"),
    hist_subtype = c(none = "#888888", pb = "#32cd32",
                      int = "#670092", gf = "#1677c4",
                      pbgf = "#ff6700", pbint = "#ff0000")
  )
  x <- pals::glasbey(length(unique(slide_ids)))
  color_scales$slide_id = setNames(x, unique(slide_ids))
  return(color_scales)
}


##################################################
# Data processing
##################################################

process_input <- function(xlsx_file,
                          metadata_sheet = "SegmentProperties",
                          count_data_sheet = "BioProbeCountMatrix",
                          metadata_id_key = "SegmentDisplayName",
                          negprobe_lod_quantile = 0.9) {
  # read input data
  metadata <- read_excel(xlsx_file, sheet=metadata_sheet)
  counts <- read_excel(xlsx_file, sheet=count_data_sheet)
  
  # check input data
  metadata_ids <- select(metadata, !!metadata_id_key) %>% pull()
  count_aoi_cols <- select(counts, all_of(metadata_ids)) %>% colnames()
  count_meta_cols <- select(counts, !all_of(count_aoi_cols)) %>% colnames()
  stopifnot(metadata_ids == count_aoi_cols)
  
  # aoi_id (study, slide, roi, segment)
  metadata <- metadata %>%
    mutate(aoi_id = sprintf("%s_s%02d_r%s_%s", study_id, slide_id, ROILabel, segment))
  colnames(counts) <- c(count_meta_cols, metadata$aoi_id)
  
  # annotate negative probes (background)
  counts <- mutate(counts, probe_type = ifelse(CodeClass == "Negative", 0, 1))
  
  # signal vs background AUC
  compute_auc <- function(responses, predictors) {
    # calculate AUC of signal vs background (negative probe)
    myroc <- pROC::roc(responses, predictors, levels=c(0, 1), direction="<")
    myauc <- as.numeric(myroc$auc)
    tibble(auc = myauc)
    #myci <- pROC::ci.auc(myroc)
    #tibble(auc=myauc, auc_ci_lo=myci[1], auc_ci_hi=myci[3])
  }
  aoi_aucs <- counts %>%
    select(probe_type, metadata$aoi_id) %>%
    pivot_longer(metadata$aoi_id, names_to="aoi_id", values_to="count") %>%
    group_by(aoi_id) %>%
    summarize(compute_auc(probe_type, count))
  metadata <- metadata %>%
    inner_join(aoi_aucs, by=c("aoi_id"="aoi_id"))
  
  # background signal statistics
  compute_negprobe_stats <- function(m) {
    m %>%
      select(probe_type, metadata$aoi_id) %>%
      filter(probe_type == 0) %>%
      pivot_longer(metadata$aoi_id, names_to="aoi_id", values_to="count") %>%
      group_by(aoi_id) %>%
      summarise(negprobe_counts = sum(count),
                negprobe_geomean = exp(mean(log(count))),
                negprobe_lod = quantile(count, negprobe_lod_quantile, names=FALSE))
  }
  # append negative probe stats to metadata
  metadata <- metadata %>% 
    inner_join(compute_negprobe_stats(counts), by=c("aoi_id"="aoi_id"))
  
  # number of probes expressed over background
  num_gene_probes <- sum(counts$probe_type == 1)
  aoi_expressed <- counts %>%
    select(probe_type, metadata$aoi_id) %>%
    filter(probe_type == 1) %>%
    pivot_longer(metadata$aoi_id, names_to="aoi_id", values_to="count") %>%
    inner_join(select(metadata, aoi_id, negprobe_lod), by=c("aoi_id"="aoi_id")) %>%
    group_by(aoi_id) %>%
    summarize(pct_expressed=sum(count > negprobe_lod) / num_gene_probes)
  metadata <- metadata %>% 
    inner_join(aoi_expressed, by=c("aoi_id"="aoi_id"))
  
  # compute upper quartile across libraries
  metadata$q3 <- counts %>%
    filter(probe_type == 1) %>%
    summarize(across(metadata$aoi_id, quantile, 0.75, names = FALSE)) %>%
    unlist()
  
  # metadata$q2 <- counts %>%
  #   filter(probe_type == 1) %>%
  #   summarize(across(metadata$aoi_id, quantile, 0.50, names = FALSE)) %>%
  #   unlist()
  
  # calculate signal to noise ratios
  metadata$num_counts <- colSums(counts[, metadata$aoi_id])
  metadata$snr <- metadata$num_counts / metadata$negprobe_counts
  metadata$q3snr <- metadata$q3 / metadata$negprobe_geomean
  metadata$q3diff <- pmax(0, metadata$q3 - metadata$negprobe_geomean)
  
  list(metadata = metadata,
       counts = counts)
}

geomean <- function(x) {
  return(exp(mean(log(x))))
}

merge_probes_to_genes <- function(x, sample_cols) {
  x <- x %>%
    group_by(entrez_id) %>%
    summarise(
      gene_symbol = unique(gene_symbol),
      nprobes = n(),
      probe_ids = str_flatten(probe_id, collapse=","),
      probe_type = min(probe_type),
      frac_expr = mean(frac_expr),
      frac_expr_min = min(frac_expr),
      frac_expr_max = max(frac_expr),
      across(all_of(sample_cols), geomean)
    )
}

calc_frac_expr <- function(m, bg_lod) {
  # subtract limit of detection and clip at zero
  x <- sweep(m, 2, bg_lod, FUN="-")
  x <- apply(x, 2, pmax, 0)
  # true/false matrix of expressed/not expressed
  x <- x > 0
  # calc fraction of samples that express each probe
  frac_expr <- (rowSums(x > 0) / ncol(x))
  frac_expr
}

background_subtract <- function(m, offset=0, min.count=1) {
  x <- sweep(m, 2, offset, FUN="-")
  x <- apply(x, 2, pmax, min.count)
  x
}

normalize_quantile <- function(x) {
  # convert to cpm
  num_counts <- colSums(x)
  xcpm <- sweep(x * 1e6, 2, num_counts, "/")
  
  # use row sums as a tiebreaker for genes with equal counts
  row_rank <- rank(rowSums(xcpm), ties.method="random")
  row_ord <- apply(xcpm, 2, order, row_rank)
  # rank matrix with ties broken by row sums
  xrank <- apply(row_ord, 2, order)
  
  # now standard quantile norm (geomean)
  xsort <- apply(xcpm, 2, sort)
  xgeomean <- apply(xsort, 1, geomean)

  index_to_value <- function(my_index, my_value){
    return(my_value[my_index])
  }
  
  xnorm <- apply(xrank, 2, index_to_value, xgeomean)
  return(xnorm)
  # convert back to raw counts
  #xraw <- sweep(xnorm / 1e6, 2, num_counts, "*")
  #xraw <- apply(xraw, 2, pmax, 1)
  #return(xraw)
}

get_bgsub_counts <- function(x, sample_meta) {
  count_meta <- select(x, -all_of(sample_meta$aoi_id))
  m <- select(x, all_of(sample_meta$aoi_id))
  m <- background_subtract(m, sample_meta$negprobe_geomean)
  m <- bind_cols(count_meta, m)
  return(m)
}

get_bgsub_qnorm_cpm <- function(x, sample_meta) {
  count_meta <- select(x, -all_of(sample_meta$aoi_id))
  m <- select(x, all_of(sample_meta$aoi_id))
  m <- background_subtract(m, sample_meta$negprobe_geomean)
  m <- normalize_quantile(m)
  m <- bind_cols(count_meta, m)
  return(m)
}


# read input data
ipmn <- process_input(input_xlsx_file, 
                      negprobe_lod_quantile = qc_negprobe_lod_quantile)

colnames(ipmn$metadata)
table(ipmn$metadata$slide_id, ipmn$metadata$histpath)

#
# samples 
#
# filter aois based on qc parameters
metadata <- ipmn$metadata %>% 
  filter(segment %in% analysis_segments) %>%
  mutate(qcfail = (num_counts <= qc_min_counts | auc <= qc_min_auc | 
                     negprobe_geomean <= qc_min_negprobe_geomean))
fmetadata <- filter(metadata, !qcfail)

# filter aois in probe matrix
counts_probe <- ipmn$counts %>%
  select(-all_of(ipmn$metadata$aoi_id), all_of(fmetadata$aoi_id)) %>%
  rename(probe_id = ProbeName,
         entrez_id = GeneID) %>%
  select(probe_id, probe_type, entrez_id, all_of(fmetadata$aoi_id))

# setup samples
samples <- fmetadata
samples$slide_id <- factor(samples$slide_id)
samples$aoi_id <- factor(samples$aoi_id)
samples$path <- factor(samples$path, levels=c("lgd", "hgd"))
samples$histology <- factor(samples$histology)
samples$carcinoma <- factor(samples$carcinoma)
samples$histpath <- factor(samples$histpath)
samples$roi_group <- factor(samples$roi_group)

table(samples$segment, samples$histology)
table(samples$histpath)
summary(samples$pct_expressed)

# color scales for plots
color_scales <- get_color_scales(unique(samples$slide_id))


#
# count data
# 
# separate negative probes
counts_negprobe <- counts_probe %>%
  filter(probe_type == 0) %>%
  mutate(entrez_id = 99990000 + as.integer(probe_id),
         gene_symbol = paste0("GEOMX_CTA_NEGPROBE_", probe_id), 
         .after = "entrez_id")
counts_probe <- counts_probe %>%
  filter(probe_type == 1)

# annotate gene probes
gene_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                      keys=as.character(counts_probe$entrez_id),
                                      column="SYMBOL", 
                                      keytype="ENTREZID", 
                                      multiVals="first")
counts_probe <- counts_probe %>% 
  mutate(gene_symbol = gene_symbols,
         .after = "entrez_id")

# combine gene probes with negative probes
counts_probe <- bind_rows(counts_probe, counts_negprobe)

# fraction of aois expressing each genes
frac_expr <- calc_frac_expr(select(counts_probe, samples$aoi_id), samples$negprobe_lod)
counts_probe <- counts_probe %>%
  mutate(frac_expr = frac_expr, .after = "gene_symbol") %>%
  group_by(entrez_id) %>%
  mutate(frac_expr_avg = mean(frac_expr)) %>%
  ungroup() %>%
  mutate(expressed = case_when(probe_type == 0 ~ "bg",
                               frac_expr_avg <= qc_min_frac_expr ~ "no",
                               TRUE ~ "yes"))
summary(counts_probe$frac_expr)
table(counts_probe$expressed)


# filter probes with low expression
# fcounts_probe <- filter(counts_probe, (probe_type == 0) | (frac_expr > qc_min_frac_expr))
fcounts_probe <- filter(counts_probe, probe_type == 1, frac_expr_avg > qc_min_frac_expr)

#
# normalized count matrices
#
#
# bgsub qnorm normalization
#
probe_bgsub_qnorm_cpm <- get_bgsub_qnorm_cpm(fcounts_probe, samples)
gene_bgsub_qnorm_cpm <- merge_probes_to_genes(probe_bgsub_qnorm_cpm, samples$aoi_id)

# gene metadata
gene_meta <- select(gene_bgsub_qnorm_cpm, -samples$aoi_id)
nrow(gene_meta)

# normalized count data
gene_bgsub_qnorm_cpm <- select(gene_bgsub_qnorm_cpm, gene_symbol, all_of(samples$aoi_id)) %>%
  column_to_rownames("gene_symbol")
gene_bgsub_qnorm_dgelist <- DGEList(counts = gene_bgsub_qnorm_cpm,
                                    samples = samples,
                                    genes = gene_meta)
gene_bgsub_qnorm_lcpm <- log2(gene_bgsub_qnorm_cpm)

#
# write processed data
#
write_xlsx(
  list(unfiltered_metadata = metadata,
       unfiltered_counts_probe = counts_probe,
       unfiltered_counts_negprobe = counts_negprobe,
       metadata = fmetadata,
       gene_meta = gene_meta,
       gene_expr_norm_log2 = as_tibble(gene_bgsub_qnorm_lcpm, rownames="gene_symbol")),
  file.path(working_dir, processed_xlsx_file)
)

##################################################
#
# Matrix for downstream analysis
#
##################################################

norm_count_dgelist <- gene_bgsub_qnorm_dgelist
norm_count_lcpm <- gene_bgsub_qnorm_lcpm

##################################################
#
# Setup gene annotation
#
##################################################

read_cptac_pdac_data <- function(filename, log2fc_cutoff, padj_cutoff) {
  rna <- read_excel(filename, sheet="Diff_exp_RNA")
  rna <- rna %>%
    dplyr::rename(gene_symbol = GeneSymbol,
                  log2fc = MedianLog2FoldChange,
                  pval = `p value`,
                  fdr = FDR) %>%
    dplyr::mutate(analysis = "cptac_rna", 
                  de = case_when(fdr > padj_cutoff ~ "no",
                                 abs(log2fc) <= log2fc_cutoff ~ "no",
                                 log2fc < -log2fc_cutoff ~ "dn",
                                 log2fc > log2fc_cutoff ~ "up",
                                 TRUE ~ "no"))
  
  # GeneSymbol	MedianLog2FoldChange	p value	FDR	
  # MedianLog2FoldChange.duct	p value.duct	FDR.duct	
  # Expression Coefficient	p value of adjusted expression	FDR of adjusted expression
  # add additional data?
  prot <- read_excel(cptac_panc_file, sheet="Diff_exp_prot")
  prot <- prot %>%
    dplyr::rename(
      gene_symbol = GeneSymbol,
      log2fc = MedianLog2FoldChange,
      pval = `p value`,
      fdr = FDR
    ) %>%
    select(gene_symbol, log2fc, pval, fdr) %>%
    mutate(analysis = "cptac_prot",
           de = case_when(fdr > padj_cutoff ~ "no",
                          abs(log2fc) <= log2fc_cutoff ~ "no",
                          log2fc < -log2fc_cutoff ~ "dn",
                          log2fc > log2fc_cutoff ~ "up",
                          TRUE ~ "no"))
  return(bind_rows(rna, prot))
}

get_cptac_gene_sets <- function(cptac_de_data) {
  cptac_rna_up <- cptac_de_data %>% filter(analysis == "cptac_rna", de == "up") %>% select(gene_symbol) %>% pull()
  cptac_rna_dn <- cptac_de_data %>% filter(analysis == "cptac_rna", de == "dn") %>% select(gene_symbol) %>% pull()
  cptac_prot_up <- cptac_de_data %>% filter(analysis == "cptac_prot", de == "up") %>% select(gene_symbol) %>% pull()
  cptac_prot_dn <- cptac_de_data %>% filter(analysis == "cptac_prot", de == "dn") %>% select(gene_symbol) %>% pull()
  cptac_de_gs <- list("CPTAC_RNA_UP" = cptac_rna_up,
                      "CPTAC_RNA_DN" = cptac_rna_dn,
                      "CPTAC_PROT_UP" = cptac_prot_up,
                      "CPTAC_PROT_DN" = cptac_prot_dn)
}


read_hpa_data <- function(hpa_path_file) {
  x <- read.table(hpa_path_file, header=TRUE, sep="\t")
  colnames(x) <- c("ens_gene_id", "gene_symbol", "cancer_type",
                   "ihc_high", "ihc_medium", "ihc_low", "ihc_not_detected",
                   "prog_fav", "unprog_fav", "prog_unfav", "unprog_unfav")
  x <- as_tibble(x) %>%
    distinct(gene_symbol, cancer_type, .keep_all=TRUE) %>%
    mutate(cancer_type = str_replace(cancer_type, " ", "_")) %>%
    mutate(prognostic = case_when(!is.na(prog_fav) ~ "fav",
                                  !is.na(prog_unfav) ~ "unfav",
                                  TRUE ~ "no")) %>%
    select(ens_gene_id, gene_symbol, cancer_type, prognostic, prog_fav, prog_unfav)
  return(x)
}

get_hpa_gene_sets <- function(hpa) {
  gs = list()
  for (catype in unique(hpa$cancer_type)) {
    unfav <- filter(hpa, cancer_type == catype, prognostic == "unfav")  %>% select(gene_symbol) %>% pull()
    fav <- filter(hpa, cancer_type == catype, prognostic == "fav")  %>% select(gene_symbol) %>% pull()
    unfav_name <- paste0("HPA_", toupper(catype), "_UNFAV_PROGNOSIS")
    fav_name <- paste0("HPA_",  toupper(catype), "_FAV_PROGNOSIS")
    x <- setNames(list(unfav, fav), c(unfav_name, fav_name))
    gs <- append(gs, x)
  }
  gs
}

# define background gene set
cta_genes <- unique(gene_meta$gene_symbol)
length(cta_genes)
# write cta genes to file (background gene set)
write.table(cta_genes, file=file.path(working_dir, cta_bg_gene_symbols_txt),
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# read cptac data
cptac_de_data <- read_cptac_pdac_data(cptac_panc_file, gsea_log2fc_cutoff, gsea_padj_cutoff)
cptac_de_data <- filter(cptac_de_data, gene_symbol %in% cta_genes)
# merge gene annotation data
cptac_de_merged <- cptac_de_data %>%
  pivot_wider(names_from = analysis, values_from = c(log2fc, pval, fdr, de))

# read hpa data
hpa <- read_hpa_data(hpa_path_file)
hpa <- filter(hpa, gene_symbol %in% cta_genes)
hpa_pdac <- filter(hpa, cancer_type == "pancreatic_cancer")

# read oncotarget data
oncotarget_data <- read_excel(oncotarget_panc_file)
oncotarget_data <- oncotarget_data %>% filter(gene_symbol %in% cta_genes)
oncotarget_data <- oncotarget_data %>% filter(abs(log_ratio) > 1, fdr < 0.01)
oncotarget_up <- oncotarget_data %>% filter(log_ratio > 0) %>% pull(gene_symbol)
oncotarget_dn <- oncotarget_data %>% filter(log_ratio < 0) %>% pull(gene_symbol)
oncotarget_gs <- list("MAO_ONCOTARGET_PDAC_RNA_UP" = oncotarget_up,
                      "MAO_ONCOTARGET_PDAC_RNA_DN" = oncotarget_dn)

gene_annot <- gene_meta %>%
  left_join(hpa_pdac, by=c("gene_symbol"="gene_symbol")) %>%
  left_join(cptac_de_merged, by=c("gene_symbol" = "gene_symbol")) %>%
  mutate(de_cptac_rna = replace_na(de_cptac_rna, "na"),
         de_cptac_prot = replace_na(de_cptac_prot, "na"))

# subtype specific scaled expression
inthgd_aois <- samples %>% filter(histpath == "panck_inthgd") %>% pull(aoi_id)
intlgd_aois <- samples %>% filter(histpath == "panck_intlgd") %>% pull(aoi_id)
pbhgd_aois <- samples %>% filter(histpath == "panck_pbhgd") %>% pull(aoi_id)
pblgd_aois <- samples %>% filter(histpath == "panck_pblgd") %>% pull(aoi_id)
x <- t(scale(t(norm_count_lcpm)))
gene_annot <- bind_cols(gene_annot,
                        inthgd_scaled_mean = rowMeans(x[, inthgd_aois]),
                        intlgd_scaled_mean = rowMeans(x[, intlgd_aois]),
                        pbhgd_scaled_mean = rowMeans(x[, pbhgd_aois]),
                        pblgd_scaled_mean = rowMeans(x[, pblgd_aois]))

##################################################
#
# Quality Control Plots
#
##################################################

table(metadata$segment)
table(metadata$segment, metadata$slide_id)
summary(metadata$AOINucleiCount)
summary(metadata$num_counts)
summary(metadata$pct_expressed)
table(metadata$qcfail, metadata$segment)
table(counts_probe$frac_expr > qc_min_frac_expr)
nrow(gene_meta)
length(unique(ipmn$counts$GeneID))

# figure: AOINucleiCount versus num_counts
x <- cor.test(metadata$AOINucleiCount, metadata$num_counts)
corlab <- paste0("r = ", round(x$estimate, 2))
p <- ggplot(metadata, aes(x=AOINucleiCount, y=num_counts, color=qcfail)) +
  geom_point() +
  geom_smooth(method = 'lm', 
              formula = y ~ x, 
              se = FALSE, 
              color = "black", 
              linetype = "dashed") +
  scale_x_continuous(trans="log10") +
  scale_y_continuous(trans="log10") +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "grey")) +
  xlab("Nuclei Count") +
  ylab("Total Counts") +
  labs(color = "QC Filtering") +
  annotate("text", x=1e3, y=1e4, label = corlab) +
  theme_minimal()
p
f <- file.path(plot_dir, "scatter_nuclei_vs_counts.pdf")
ggsave(f, plot=p, device="pdf", width = 4, height = 3)

# figure: AOISurfaceArea versus num_counts
x <- cor.test(metadata$AOISurfaceArea, metadata$num_counts)
corlab <- paste0("r = ", round(x$estimate, 2))
p <- ggplot(metadata, aes(x=AOISurfaceArea, y=num_counts, color=qcfail)) +
  geom_point() +
  geom_smooth(method = 'lm', 
              formula = y ~ x, 
              se = FALSE, 
              color = "black", 
              linetype = "dashed") +
  scale_x_continuous(trans="log10") +
  scale_y_continuous(trans="log10") +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "grey")) +
  xlab("Surface Area (um2)") +
  ylab("Total Counts") +
  labs(color = "QC Filtering") +
  annotate("text", x=1e5, y=1e4, label = corlab) +
  theme_minimal()
p
f <- file.path(plot_dir, "scatter_area_vs_counts.pdf")
ggsave(f, plot=p, device="pdf", width = 4, height = 3)

# figure: total counts vs snAUC
p <- ggplot(metadata, aes(x=num_counts, y=auc, color=qcfail)) +
  geom_point() + 
  geom_hline(yintercept = qc_min_auc, linetype = "dashed", color = "red") +
  geom_vline(xintercept = qc_min_counts, linetype = "dashed", color = "red") +
  scale_x_continuous(trans="log10") +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "grey")) +
  xlab("Total Counts") +
  ylab("snAUC") +
  labs(color = "Segment") +
  theme_minimal()
#  facet_wrap(~factor(segment, levels=c("cd45", "sma", "panck")))
p
f <- file.path(plot_dir, "scatter_counts_vs_auc_bysegment.pdf")
ggsave(f, plot=p, device="pdf", width = 4, height = 3)


# figure: boxplot path versus total counts
p <- ggplot(filter(fmetadata, segment == "panck"), aes(x=path, y=num_counts)) +
  geom_boxplot(aes(fill=path), width=0.5, outlier.shape=NA) +
  geom_point(aes(fill=path), color="black", size=2, position=position_jitterdodge(jitter.width=0.15)) +
  scale_y_continuous(trans="log10") +
  scale_fill_manual(values=color_scales$path) +
  scale_color_manual(values=color_scales$path) +
  xlab("Pathology") +
  ylab("Total Counts") +
  labs(fill = "Pathology") +
  theme_minimal() +
  theme(legend.position="bottom")
  #facet_grid(cols = vars(slide_id)) 
p
f <- file.path(plot_dir, "boxplot_path_vs_counts.pdf")
ggsave(f, plot=p, device="pdf", width = 2, height = 3)

# figure: boxplot path versus auc
p <- ggplot(filter(fmetadata, segment == "panck"), aes(x=path, y=auc)) +
  geom_boxplot(aes(fill=path), width=0.5, outlier.shape=NA) +
  geom_point(aes(fill=path), color="black", size=2, position=position_jitterdodge(jitter.width=0.15)) +
  scale_y_continuous(trans="log10") +
  scale_fill_manual(values=color_scales$path) +
  scale_color_manual(values=color_scales$path) +
  xlab("Pathology") +
  ylab("snAUC") +
  labs(fill = "Pathology") +
  theme_minimal() +
  theme(legend.position="bottom")
  #facet_grid(cols = vars(slide_id)) 
p
f <- file.path(plot_dir, "boxplot_path_vs_auc.pdf")
ggsave(f, plot=p, device="pdf", width = 2, height = 3)

# figure: boxplot path versus aoi density
p <- ggplot(filter(fmetadata, segment == "panck"), aes(x=path, y=100 * AOINucleiCount / AOISurfaceArea)) +
  geom_boxplot(aes(fill=path), width=0.5, outlier.shape=NA) +
  geom_point(aes(fill=path), color="black", size=2, position=position_jitterdodge(jitter.width=0.15)) +
  scale_y_continuous(trans="log10") +
  scale_fill_manual(values=color_scales$path) +
  scale_color_manual(values=color_scales$path) +
  xlab("Pathology") +
  ylab("AOI Density (Nuclei / 100um2)") +
  labs(fill = "Pathology") +
  theme_minimal() +
  theme(legend.position="bottom")
#facet_grid(cols = vars(slide_id)) 
p
f <- file.path(plot_dir, "boxplot_path_vs_aoi_density.pdf")
ggsave(f, plot=p, device="pdf", width = 2, height = 3)


# figure: boxplot histology versus aoi density
p <- ggplot(filter(fmetadata, segment == "panck"), aes(x=histology, y=100 * AOINucleiCount / AOISurfaceArea)) +
  geom_boxplot(aes(fill=histology), width=0.5, outlier.shape=NA) +
  geom_point(aes(fill=histology), color="black", size=2, position=position_jitterdodge(jitter.width=0.15)) +
  scale_y_continuous(trans="log10") +
  scale_fill_manual(values=color_scales$histology) +
  scale_color_manual(values=color_scales$histology) +
  xlab("Pathology") +
  ylab("AOI Density (Nuclei / 100um2)") +
  labs(fill = "Pathology") +
  theme_minimal() +
  theme(legend.position="bottom")
p
f <- file.path(plot_dir, "boxplot_histology_vs_aoidensity.pdf")
ggsave(f, plot=p, device="pdf", width = 2, height = 3)


# figure: boxplot histology versus total counts
p <- ggplot(filter(fmetadata, segment == "panck"), aes(x=histology, y=num_counts)) +
  geom_boxplot(aes(fill=histology), width=0.5, outlier.shape=NA) +
  geom_point(aes(fill=histology), color="black", size=2, position=position_jitterdodge(jitter.width=0.15)) +
  scale_y_continuous(trans="log10") +
  scale_fill_manual(values=color_scales$histology) +
  scale_color_manual(values=color_scales$histology) +
  xlab("Histology") +
  ylab("Total Counts") +
  labs(fill = "Histology") +
  theme_minimal() +
  theme(legend.position="bottom")
p
f <- file.path(plot_dir, "boxplot_histology_vs_counts.pdf")
ggsave(f, plot=p, device="pdf", width = 2, height = 3)

# figure: boxplot histology versus auc
p <- ggplot(filter(fmetadata, segment == "panck"), aes(x=histology, y=auc)) +
  geom_boxplot(aes(fill=histology), width=0.5, outlier.shape=NA) +
  geom_point(aes(fill=histology), color="black", size=2, position=position_jitterdodge(jitter.width=0.15)) +
  scale_y_continuous(trans="log10") +
  scale_fill_manual(values=color_scales$histology) +
  scale_color_manual(values=color_scales$histology) +
  xlab("Histology") +
  ylab("snAUC") +
  labs(fill = "Histology") +
  theme_minimal() +
  theme(legend.position="bottom")
#facet_grid(cols = vars(slide_id)) 
p
f <- file.path(plot_dir, "boxplot_histology_vs_auc.pdf")
ggsave(f, plot=p, device="pdf", width = 2, height = 3)


# figure: gene filtering fraction expressed cutoff
p <- ggplot(filter(counts_probe, probe_type == 1), aes(x=frac_expr)) +
  geom_histogram(binwidth=0.05, color="black", fill="grey") +
  geom_vline(xintercept = qc_min_frac_expr, linetype="dashed", color="red") +
  theme_minimal() +
  xlab("Fraction Expressed Above Background") +
  ylab("# of Probes") +
  facet_wrap(~ probe_type)
p
f <- file.path(plot_dir, "histogram_frac_expr.pdf")
ggsave(f, plot=p, device="pdf", width = 4, height = 4)

#
# count ridge plots
#
x <- counts_probe
y <- select(fmetadata, aoi_id, slide_id, segment, histology, histpath, path)
x <- x %>%
  pivot_longer(fmetadata$aoi_id, names_to="aoi_id", values_to="count") %>%
  inner_join(y, by=c("aoi_id"="aoi_id"))

# figure: ridge plot unnormalized counts (expressed bg/no/yes)
p <- x %>%
  ggplot(aes(x=count, y=factor(aoi_id), fill=factor(expressed))) +
  stat_density_ridges(alpha=0.7, scale=10, rel_min_height=0.001, quantile_lines=TRUE) + 
  scale_x_continuous(trans="log10") +
  xlab("Raw Counts") +
  ylab("AOIs") +
  labs(fill = "Expressed") +
  theme_ridges() + 
  theme(axis.text = element_text(size=5))
p
f <- file.path(plot_dir, "ridge_plot_rawcounts_by_aoi.pdf")
ggsave(f, plot=p, device="pdf", width = 8, height = 8)


# cpm (without background subtract)
x <- counts_probe %>% select(fmetadata$aoi_id)
x <- sweep(x, 2, colSums(x) / 1e6, "/")
x <- bind_cols(select(counts_probe, -fmetadata$aoi_id), x)
y <- select(fmetadata, aoi_id, slide_id, segment, histology, histpath, path)
x <- x %>%
  pivot_longer(fmetadata$aoi_id, names_to="aoi_id", values_to="count") %>%
  inner_join(y, by=c("aoi_id"="aoi_id"))

# ridge plot showing variability in background noise
p <- x %>% filter(expressed == "bg") %>%
  ggplot(aes(x=count, y=reorder(factor(aoi_id), count), fill=factor(slide_id))) +
  stat_density_ridges(alpha=0.7, scale=4, rel_min_height=0.001, quantile_lines=TRUE) + 
  scale_fill_manual(values = pals::cols25()) +
  xlab("CPM") +
  ylab("AOIs") +
  labs(fill = "Neg. Probes") +
  theme_ridges() + 
  theme(axis.text = element_text(size=5)) + 
  facet_grid(slide_id ~ ., scales="free_y")
p
f <- file.path(plot_dir, "ridge_plot_bg_cpm_by_aoi.pdf")
ggsave(f, plot=p, device="pdf", width = 8, height = 8)


# plot: compared median counts to bg counts
x <- x %>%
  group_by(aoi_id, expressed) %>%
  summarise(slide_id = unique(slide_id),
            mediancpm = median(count),
            q3cpm = quantile(count, 0.75)) %>%
  ungroup()
p <- ggplot(x, aes(x=mediancpm, y=reorder(factor(aoi_id), mediancpm), color=factor(expressed))) +
  geom_point(size=2) +
  scale_color_manual(values = pals::cols25()) +
  labs(color = "Expressed", x="Median CPM", y="ROI") + 
  theme_ridges() + 
  theme(axis.text.y = element_blank())
p
f <- file.path(plot_dir, "qc_point_plot_cpm_by_aoi.pdf")
ggsave(f, plot=p, device="pdf")


x <- x %>% 
  pivot_wider(names_from = "expressed", values_from=c("mediancpm", "q3cpm"))
corlab <- paste0("r = ", round(cor.test(x$mediancpm_bg, x$mediancpm_yes)$estimate, 2))
p <- ggplot(x, aes(x=mediancpm_bg, y=mediancpm_yes)) +
  geom_point() + 
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color="red", linetype="dashed") +
  annotate("text", x=60, y=60, label = corlab) +
  theme_minimal() + 
  labs(x="median CPM background", y="median CPM expressed")
p
f <- file.path(plot_dir, "scatter_qc_bg_vs_expressed.pdf")
ggsave(f, plot=p, device="pdf", width = 4, height = 3)


# cpm (with background subtract)
x <- counts_probe %>% select(fmetadata$aoi_id)
x <- background_subtract(x, fmetadata$negprobe_geomean)
x <- sweep(x, 2, colSums(x) / 1e6, "/")
x <- bind_cols(select(counts_probe, -fmetadata$aoi_id), x)
y <- select(fmetadata, aoi_id, slide_id, segment, histology, histpath, path)
x <- x %>%
  pivot_longer(fmetadata$aoi_id, names_to="aoi_id", values_to="count") %>%
  inner_join(y, by=c("aoi_id"="aoi_id"))
x <- x %>%
  group_by(aoi_id, expressed) %>%
  summarise(slide_id = unique(slide_id),
            mediancpm = median(count),
            q3cpm = quantile(count, 0.75)) %>%
  ungroup()
p <- ggplot(x, aes(x=mediancpm, y=reorder(factor(aoi_id), mediancpm), color=factor(expressed))) +
  geom_point(size=2) +
  scale_color_manual(values = pals::cols25()) +
  labs(color = "Expressed", x="Median CPM", y="ROI") + 
  theme_ridges() + 
  theme(axis.text.y = element_blank())
p
f <- file.path(plot_dir, "qc_point_plot_bgsub_cpm_by_aoi.pdf")
ggsave(f, plot=p, device="pdf")


x <- x %>% 
  pivot_wider(names_from = "expressed", values_from=c("mediancpm", "q3cpm"))
corlab <- paste0("r = ", round(cor.test(x$mediancpm_bg, x$mediancpm_yes)$estimate, 2))
p <- ggplot(x, aes(x=mediancpm_bg, y=mediancpm_yes)) +
  geom_point() + 
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color="red", linetype="dashed") + 
  annotate("text", x=10, y=40, label = corlab) +
  theme_minimal() + 
  labs(x="median CPM background", y="median CPM expressed")
p
f <- file.path(plot_dir, "scatter_qc_bgsub_bg_vs_expressed.pdf")
ggsave(f, plot=p, device="pdf", width = 4, height = 3)


# cpm (with qnorm)
x <- get_bgsub_qnorm_cpm(counts_probe, samples)
x <- x %>%
  pivot_longer(fmetadata$aoi_id, names_to="aoi_id", values_to="count") %>%
  inner_join(y, by=c("aoi_id"="aoi_id"))
x <- x %>%
  group_by(aoi_id, expressed) %>%
  summarise(slide_id = unique(slide_id),
            mediancpm = median(count),
            q3cpm = quantile(count, 0.75)) %>%
  ungroup()
p <- ggplot(x, aes(x=mediancpm, y=reorder(factor(aoi_id), mediancpm), color=factor(expressed))) +
  geom_point(size=2) +
  scale_color_manual(values = pals::cols25()) +
  labs(color = "Expressed", x="Median CPM", y="ROI") + 
  theme_ridges() + 
  theme(axis.text.y = element_blank())
p
f <- file.path(plot_dir, "qc_point_plot_bgsub_qnorm_cpm_by_aoi.pdf")
ggsave(f, plot=p, device="pdf")


#
# normalized count ridge plots
#
# qnorm
x <- get_bgsub_qnorm_cpm(counts_probe, samples)
x <- DGEList(counts = log2(select(x, samples$aoi_id)), 
             samples = samples, 
             genes = select(x, -samples$aoi_id))
x <- bind_cols(x$genes, x$counts)
y <- select(samples, aoi_id, slide_id, segment, histology, path)
x <- x %>%
  pivot_longer(samples$aoi_id, names_to="aoi_id", values_to="count") %>%
  inner_join(y, by=c("aoi_id"="aoi_id"))

# by aoi
p <- x %>%
  ggplot(aes(x=count, y=factor(aoi_id), fill=factor(expressed))) +
  stat_density_ridges(alpha=0.7, scale=10, rel_min_height=0.001, quantile_lines=TRUE) + 
  xlab("log2(Counts per million reads)") +
  ylab("Slides") +
  labs(fill = "Expressed") +
  theme_ridges() + 
  theme(axis.text = element_text(size=5))
p
f <- file.path(plot_dir, "ridge_plot_unfiltered_probes_by_aoi.pdf")
ggsave(f, plot=p, device="pdf", width = 8, height = 8)

# by path
p <- x %>%
  ggplot(aes(x=count, y=factor(path), fill=factor(expressed))) +
  stat_density_ridges(alpha=0.7, scale=10, rel_min_height=0.001, quantile_lines=TRUE) + 
  xlab("log2(Counts per million reads)") +
  ylab("Slides") +
  labs(fill = "Expressed") +
  theme_ridges() + 
  theme(axis.text = element_text(size=5))
p
f <- file.path(plot_dir, "ridge_plot_unfiltered_probes_by_path.pdf")
ggsave(f, plot=p, device="pdf", width = 8, height = 8)

# by histology
p <- x %>%
  ggplot(aes(x=count, y=factor(histology), fill=factor(expressed))) +
  stat_density_ridges(alpha=0.7, scale=10, rel_min_height=0.001, quantile_lines=TRUE) + 
  xlab("log2(Counts per million reads)") +
  ylab("Slides") +
  labs(fill = "Expressed") +
  theme_ridges() + 
  theme(axis.text = element_text(size=5))
p
f <- file.path(plot_dir, "ridge_plot_unfiltered_probes_by_histology.pdf")
ggsave(f, plot=p, device="pdf", width = 8, height = 8)


#
# filtered normalized gene ridge plots
#
x <- norm_count_lcpm
x <- as_tibble(x, rownames("gene_symbol"))
y <- select(fmetadata, aoi_id, slide_id, segment, histology, path, histpath)
x <- x %>%
  pivot_longer(fmetadata$aoi_id, names_to="aoi_id", values_to="count") %>%
  inner_join(y, by=c("aoi_id"="aoi_id"))

p <- x %>%
  ggplot(aes(x=count, y=factor(aoi_id), fill=factor(slide_id))) +
  stat_density_ridges(alpha=0.7, scale=10, rel_min_height=0.001, quantile_lines=TRUE) + 
  scale_fill_manual(values = color_scales$slide_id) +
  xlab("log2(cpm)") +
  ylab("AOIs") +
  labs(fill = "Slide #") +
  theme_ridges() +
  theme(axis.text = element_text(size=5))
p
f <- file.path(plot_dir, "ridge_plot_log2cpm_by_aoislide.pdf")
ggsave(f, plot=p, device="pdf", width = 8, height = 8)

p <- filter(x, segment == "panck") %>%
  ggplot(aes(x=count, y=path, fill=path)) +
  stat_density_ridges(alpha=0.7, scale=3, rel_min_height=0.001, quantile_lines=TRUE) + 
  scale_fill_manual(values = color_scales$path) +
  xlab("log2(cpm)") +
  ylab("Pathology") +
  labs(fill = "Pathology") +
  theme_ridges() + 
  theme(legend.position="bottom")
p
f <- file.path(plot_dir, "ridge_plot_log2cpm_by_path.pdf")
ggsave(f, plot=p, device="pdf", width = 4, height = 4)

p <- filter(x, segment == "panck") %>%
  ggplot(aes(x=count, y=histology, fill=histology)) +
  stat_density_ridges(alpha=0.7, scale=3, rel_min_height=0.001, quantile_lines=TRUE) + 
  scale_fill_manual(values = color_scales$histology) +
  xlab("log2(cpm)") +
  ylab("Histology") +
  labs(fill = "Histology") +
  theme_ridges() + 
  theme(legend.position="bottom")
p
f <- file.path(plot_dir, "ridge_plot_log2cpm_by_hist.pdf")
ggsave(f, plot=p, device="pdf", width = 4, height = 4)


##################################################
#
# Dimensionality Reduction Analysis
#
##################################################

run_tsne <- function(x, meta, perp=9) {
  # tsne normalizes input by default
  x <- t(x)
  tsne_fit <- Rtsne(x, 
                    perplexity = perp,
                    theta = 0.0,
                    max_iter = 10000,
                    normalize = TRUE)
  # tsne_fit <- Rtsne(x, perplexity=perp)
  tsne_tbl <- tsne_fit$Y
  colnames(tsne_tbl) <- c("tsne1", "tsne2")
  tsne_tbl <- tsne_tbl %>%
    as_tibble() %>%
    mutate(aoi_id=meta$aoi_id) %>%
    inner_join(meta, by=c("aoi_id"="aoi_id"))
  return(list(tbl=tsne_tbl))
}

run_pca <- function(x, meta) {
  #x <- scale(t(x))
  #x <- normalize_input(t(x))
  x <- t(x)
  pr <- prcomp(x, center=TRUE, scale=FALSE, rank=10)
  pca_tbl <- as.data.frame(pr$x)
  pca_tbl <- pca_tbl %>%
    as_tibble(rownames="aoi_id") %>%
    inner_join(meta, by=c("aoi_id"="aoi_id"))
  return(list(pr = pr, tbl=pca_tbl))
}

# quantile normalization
y <- norm_count_lcpm
pca_res <- run_pca(y, samples)
pca_tbl <- pca_res$tbl
# tsne_res <- run_tsne(y, samples, perp=11)
# tsne_tbl <- tsne_res$tbl
#write_xlsx(tsne_tbl, file.path(working_dir, tsne_results_xlsx_file))
tsne_tbl <- read_xlsx(file.path(working_dir, tsne_results_xlsx_file), sheet = "Sheet1")

# pca
pr <- pca_res$pr
pct <- round(pr$sdev / sum(pr$sdev) * 100, 2)
pct <- paste(colnames(pr$x), " (", paste(as.character(pct), "%", ")", sep=""), sep="")
p <- ggplot(pca_tbl, aes(x=PC1, y=PC2, color=histology, shape=path)) +
  geom_point(size=3) +
  scale_color_manual(values=color_scales$histology) +
  theme_minimal() +
  labs(x=pct[1], y=pct[2])
p
f <- file.path(plot_dir, "pca_histology.pdf")
ggsave(f, plot=p, device="pdf", width=6, height=5)

p <- ggplot(pca_tbl, aes(x=PC1, y=PC2, color=slide_id, shape=path)) +
  geom_point(size=3) +
  scale_color_manual(values=color_scales$slide_id) +
  theme_minimal() +
  labs(x=pct[1], y=pct[2])
p
f <- file.path(plot_dir, "pca_slide.pdf")
ggsave(f, plot=p, device="pdf", width=6, height=5)

p <- ggplot(pca_tbl, aes(x=PC1, y=PC2, color=carcinoma, shape=path)) +
  geom_point(size=3) +
  scale_color_manual(values=color_scales$carcinoma) +
  theme_minimal() +
  labs(x=pct[1], y=pct[2])
p
f <- file.path(plot_dir, "pca_carcinoma.pdf")
ggsave(f, plot=p, device="pdf", width=6, height=5)


# tsne
p <- ggplot(tsne_tbl, aes(x=tsne1, y=tsne2, color=histology, shape=path)) +
  geom_point(size=3) +
  scale_color_manual(values=color_scales$histology) +
  theme_minimal()
p
f <- file.path(plot_dir, "tsne_panck_histology.pdf")
ggsave(f, plot=p, device="pdf", width=6, height=5)

# slide
p <- ggplot(tsne_tbl, aes(x=tsne1, y=tsne2, color=slide_id, shape=path)) +
  geom_point(size=3) +
  scale_color_manual(values=color_scales$slide_id) +
  theme_minimal()
p
f <- file.path(plot_dir, "tsne_panck_slide_id.pdf")
ggsave(f, plot=p, device="pdf", width=6, height=5)

# carcinoma
p <- ggplot(tsne_tbl, aes(x=tsne1, y=tsne2, color=carcinoma, shape=path)) +
  geom_point(size=3) +
  scale_color_manual(values=color_scales$carcinoma) +
  theme_minimal()
p
f <- file.path(plot_dir, "tsne_panck_carcinoma.pdf")
ggsave(f, plot=p, device="pdf", width=6, height=5)


##################################################
#
# Differential Expression Analysis
#
##################################################

formatDEResults <- function(res, analysis, method,
                            padj_cutoff, log2fc_cutoff) {
  mutate(res, 
         de = case_when(padj > padj_cutoff ~ "no",
                        log2fc < -log2fc_cutoff ~ "dn",
                        log2fc > log2fc_cutoff ~ "up",
                        TRUE ~ "no"),
         analysis = analysis,
         method = method)
}

saveDEResults <- function(file_prefix, res, file_path=NULL) {
  file_path <- ifelse(is.null(file_path), getwd(), file_path)
  # write results to tab-separated file
  f <- paste0(file_prefix, "_results.tsv")
  write.table(res, file=file.path(file_path, f), row.names=FALSE, sep="\t", quote=FALSE)
  # save de genes
  up <- filter(res, de == "up") %>% select(gene) %>% pull()
  dn <- filter(res, de == "dn") %>% select(gene) %>% pull()
  de <- filter(res, de != "no") %>% select(gene) %>% pull()
  write(up, file=file.path(file_path, paste0(file_prefix, "_de_up_genes.txt")), sep="\n")
  write(dn, file=file.path(file_path, paste0(file_prefix, "_de_dn_genes.txt")), sep="\n")
  write(de, file=file.path(file_path, paste0(file_prefix, "_de_genes.txt")), sep="\n")
}


limmaRun <- function(y, design, contrasts, trend=FALSE) {
  fit <- lmFit(y, design)
  fit <- contrasts.fit(fit, contrasts)
  fit <- eBayes(fit, trend=trend)
  return(list(fit=fit)) 
}

limmaRenameResults <- function(res) {
  as_tibble(res, rownames="gene") %>%
    select(gene,
           log2fc=logFC,
           avgexpr=AveExpr,
           pval=P.Value,
           padj=adj.P.Val)
}

processLimma <- function(fit, contrasts, method, padj_cutoff, log2fc_cutoff,
                         file_path=NULL) {
  x <- NULL
  for (coef in colnames(contrasts)) {
    prefix <- paste(method, coef, sep="_")
    res <- topTable(fit, coef=coef, number=Inf, sort.by="none")
    res <- limmaRenameResults(res) %>%
      formatDEResults(coef, method, padj_cutoff, log2fc_cutoff)
    saveDEResults(prefix, res, file_path)
    x <- bind_rows(x, res)
  }
  return(x)
}


#
# histology
#
design <- model.matrix(~0 + histology, data=samples)
colnames(design) <- levels(samples$histology)
colnames(design)
contrasts = makeContrasts(hist_pb_vs_gf = pb - gf,
                          hist_pb_vs_int = pb - intestinal,
                          hist_int_vs_gf = intestinal - gf,
                          hist_int = intestinal - (pb + gf)/2,
                          hist_pb = pb - (intestinal + gf)/2,
                          hist_gf = gf - (pb + intestinal)/2,
                          levels=design)

# qnorm
method <- "limma_trend_qnorm"
x <- limmaRun(gene_bgsub_qnorm_lcpm, design, contrasts, trend=TRUE)
res <- processLimma(x$fit, contrasts, method, de_padj_cutoff, de_log2fc_cutoff, de_result_dir)
hist_de <- res

# merged de analyses
hist_de_merged <- select(hist_de, gene, analysis, avgexpr, log2fc, padj, de) %>%
  pivot_wider(names_from = analysis, values_from = c(log2fc, padj, de))

hist_de_merged <- hist_de_merged %>%
  mutate(subtype = case_when(de_hist_pb_vs_gf == "up" & de_hist_pb_vs_int == "up" ~ "pb",
                             de_hist_pb_vs_gf == "dn" & de_hist_int_vs_gf == "dn" ~ "gf",
                             de_hist_int_vs_gf == "up" & de_hist_pb_vs_int == "dn" ~ "int",
                             de_hist_int_vs_gf == "dn" & de_hist_pb_vs_int == "up" ~ "pbgf",
                             de_hist_pb_vs_gf == "up" & de_hist_int_vs_gf == "up" ~ "pbint",
                             de_hist_pb_vs_gf == "dn" & de_hist_int_vs_gf == "up" ~ "intgf",
                             TRUE ~ "none"))
table(hist_de_merged$subtype)
nrow(filter(hist_de_merged, subtype != "none"))


#
# grade
#
table(samples$histpath)
design <- model.matrix(~0 + histpath, data=samples)
colnames(design)
colnames(design) <- levels(samples$histpath)
colnames(design)
contrasts = makeContrasts(path_hgd_vs_lgd = (panck_pbhgd + panck_inthgd)/2 - (panck_intlgd + panck_pblgd)/2,
                          path_pbhgd_vs_intlgd = panck_pbhgd - panck_intlgd,
                          path_pbhgd_vs_pblgd = panck_pbhgd - panck_pblgd,
                          path_inthgd_vs_intlgd = panck_inthgd - panck_intlgd,
                          path_inthgd_vs_pblgd = panck_inthgd - panck_pblgd,
                          path_intlgd_vs_pblgd = panck_intlgd - panck_pblgd,
                          path_pbhgd_vs_inthgd = panck_pbhgd - panck_inthgd,
                          levels=design)

# qnorm
method <- "limma_trend_qnorm"
x <- limmaRun(gene_bgsub_qnorm_lcpm, design, contrasts, trend=TRUE)
res <- processLimma(x$fit, contrasts, method, de_padj_cutoff, de_log2fc_cutoff, de_result_dir)
grade_de <- res

# merged de analyses
grade_de_merged <- select(grade_de, gene, analysis, avgexpr, log2fc, padj, de) %>%
  pivot_wider(names_from = analysis, values_from = c(log2fc, padj, de))
grade_de_merged <- grade_de_merged %>%
  mutate(subtype = case_when(
    de_path_hgd_vs_lgd == "up" ~ "hgd",
    de_path_hgd_vs_lgd == "dn" ~ "lgd",
    de_path_inthgd_vs_intlgd == "up" ~ "inthgd",
    de_path_inthgd_vs_intlgd == "dn" ~ "intlgd",
    de_path_pbhgd_vs_pblgd == "up" ~ "pbhgd",
    de_path_pbhgd_vs_pblgd == "dn" ~ "pblgd",
    TRUE ~ "none")
  )
table(grade_de_merged$subtype)

# add to gene annotation
gene_annot <- gene_annot %>%
  left_join(select(grade_de_merged, gene, grade_subtype = subtype), by=c("gene_symbol"="gene")) %>%
  left_join(select(hist_de_merged, gene, hist_subtype = subtype), by=c("gene_symbol"="gene"))
colnames(gene_annot)


# write gene expression data
write_xlsx(list(samples = samples,
                gene_meta = gene_meta,
                gene_annot = gene_annot,
                de = bind_rows(hist_de, grade_de),
                hist_de_merged = hist_de_merged,
                grade_de_merged = grade_de_merged),
           file.path(working_dir, de_results_xlsx))

##################################################
#
#
# Heatmaps
#
#
##################################################

# heatmap functions
plot_pheatmap <- function(m, 
                          annot_col, 
                          annot_row, 
                          annot_colors) {
  pheatmap(m,
           scale = "row",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           breaks = seq(-2, 2, length.out=51),
           annotation_col = annot_col,
           annotation_row = annot_row,
           annotation_colors = annot_colors,
           border_color = NA,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "ward.D2",
           fontsize_row = 6,
           fontsize_col = 6)
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

get_row_annot <- function(genes, gene_annot) {
  gene_annot %>%
    filter(gene_symbol %in% genes) %>%
    select(gene_symbol, prognostic, de_cptac_rna, de_cptac_prot, grade_subtype, hist_subtype) %>%
    column_to_rownames("gene_symbol") %>%
    as.data.frame()
}

# 
# histology analysis heatmap
#
# union of subtype specific genes
hist_genes <- filter(hist_de_merged, subtype != "none") %>% select(gene) %>% pull()
x <- norm_count_lcpm[hist_genes,]
nrow(x)

annot_row <- get_row_annot(hist_genes, gene_annot)
annot_col <- column_to_rownames(samples, "aoi_id") %>% 
  select(histology = histology,
         pathbw = path, 
         slide_id = slide_id) %>%
  as.data.frame()
f <- file.path(plot_dir,  paste0("heatmap_histology.pdf"))
p <- plot_pheatmap(x, annot_col, annot_row, color_scales)
h <- round(3 + (length(hist_genes) / 10))
save_pheatmap_pdf(p, filename=f, width=10, height=h)

# 
# grade analysis heatmap
#
# all hgd vs all lgd
hgd_lgd_genes <- filter(grade_de_merged, de_path_hgd_vs_lgd != "no") %>% select(gene) %>% pull()
nrow(filter(grade_de_merged, de_path_hgd_vs_lgd == "up"))
nrow(filter(grade_de_merged, de_path_hgd_vs_lgd == "dn"))

x <- norm_count_lcpm[hgd_lgd_genes,]
nrow(x)
annot_row <- get_row_annot(hgd_lgd_genes, gene_annot)
annot_col <- column_to_rownames(samples, "aoi_id") %>% 
  select(histology = histology,
         pathbw = path, 
         slide_id = slide_id) %>%
  as.data.frame()
f <- file.path(plot_dir, "heatmap_hgd_vs_lgd.pdf")
p <- plot_pheatmap(x, annot_col, annot_row, color_scales)
h <- round(3 + (length(hgd_lgd_genes) / 10))
save_pheatmap_pdf(p, filename=f, width=10, height=h)
dev.off()

# union of subtype hgd genes
table(grade_de_merged$subtype)
grade_genes <- filter(grade_de_merged, subtype != "none") %>% select(gene) %>% pull()
x <- norm_count_lcpm[grade_genes,]
nrow(x)

annot_row <- get_row_annot(grade_genes, gene_annot)
annot_col <- column_to_rownames(samples, "aoi_id") %>% 
  select(histology = histology,
         pathbw = path, 
         slide_id = slide_id) %>%
  as.data.frame()
f <- file.path(plot_dir, "heatmap_grade.pdf")
p <- plot_pheatmap(x, annot_col, annot_row, color_scales)
h <- round(3 + (length(grade_genes) / 10))
save_pheatmap_pdf(p, filename=f, width=10, height=h)
dev.off()

##################################################
#
#
# Euler Plots of DE Gene Sets
#
#
##################################################

# histology
x <- filter(hist_de, de != "no")
pb_int <- filter(x, analysis == "hist_pb_vs_int", de == "up") %>% select(gene) %>% pull()
gf_int <- filter(x, analysis == "hist_int_vs_gf", de == "dn") %>% select(gene) %>% pull()
int_pb <- filter(x, analysis == "hist_pb_vs_int", de == "dn") %>% select(gene) %>% pull()
int_gf <- filter(x, analysis == "hist_int_vs_gf", de == "up") %>% select(gene) %>% pull()
pb_gf <- filter(x, analysis == "hist_pb_vs_gf", de == "up") %>% select(gene) %>% pull()
gf_pb <- filter(x, analysis == "hist_pb_vs_gf", de == "dn") %>% select(gene) %>% pull()

f <- file.path(plot_dir, "euler_histology.pdf")
pdf(f)
p <- euler(list("PB vs INT" = pb_int, 
                "PB vs GF" = pb_gf, 
                "GF vs INT" = gf_int, 
                "GF vs PB" = gf_pb,
                "INT vs PB" = int_pb,
                "INT vs GF" = int_gf),
           shape = "ellipse")
print(plot(p, 
           fills = c("#90ee90", "#1BF118", 
                     "#89cff0", "#27BBE0", 
                     "#ffb200", "#ffff00"),
           quantities = TRUE), cex.lab=2)
dev.off()

# pathology / grade
x <- filter(grade_de, de != "no")
pb_hgd_up <- filter(x, analysis == "path_pbhgd_vs_pblgd", de == "up") %>% select(gene) %>% pull()
pb_hgd_dn <- filter(x, analysis == "path_pbhgd_vs_pblgd", de == "dn") %>% select(gene) %>% pull()
int_hgd_up <- filter(x, analysis == "path_inthgd_vs_intlgd", de == "up") %>% select(gene) %>% pull()
int_hgd_dn <- filter(x, analysis == "path_inthgd_vs_intlgd", de == "dn") %>% select(gene) %>% pull()
hgd_up <- filter(x, analysis == "path_hgd_vs_lgd", de == "up") %>% select(gene) %>% pull()
hgd_dn <- filter(x, analysis == "path_hgd_vs_lgd", de == "dn") %>% select(gene) %>% pull()

f <- file.path(plot_dir, "euler_grade.pdf")
pdf(f)
p <- euler(list("PB-HGD UP" = pb_hgd_up, "PB-HGD DN" = pb_hgd_dn, 
                "INT-HGD UP" = int_hgd_up, "INT-HGD DN" = int_hgd_dn,
                "HGD UP" = hgd_up, "HGD DN" = hgd_dn),
           shape = "ellipse")
print(plot(p, 
           fills = c("#90ee90", "#1BF118", 
                     "#ffff00", "#ffaa1d",
                     "#89cff0", "#27BBE0"),
           quantities = TRUE))
dev.off()


##################################################
#
#
# Volcano plot of DE genes
#
#
##################################################

plot_de_volcano <- function(x) {
  p <- ggplot(x, aes(x=log2fc, y=-log10(padj))) + 
    geom_point(data=subset(x, de == "no"), color="#888888", size=1, alpha=0.4) + 
    geom_point(data=subset(x, de == "up"), color="#ff0000", size=2, alpha=0.8) +
    geom_point(data=subset(x, de == "dn"), color="#0000ff", size=2, alpha=0.8) +
    geom_text_repel(data=subset(x, de != "no"), color="black", size=3, aes(label=gene), max.overlaps=Inf) +
    ylab("-log10(adjusted p-value)") +
    xlab("log2(Fold Change)") +
    theme_minimal() +
    theme(legend.position="bottom") + 
    theme(axis.line = element_line(color = "black"))
  return(p)
}

x <- filter(grade_de, analysis == "path_hgd_vs_lgd")
p <- plot_de_volcano(x) +
  ggtitle("DE genes (HGD vs LGD)")
p
f <- file.path(plot_dir, "volcano_de_hgd_vs_lgd.pdf")
ggsave(f, plot=p, device="pdf", width=5, height=5)

x <- filter(grade_de, analysis == "path_pbhgd_vs_pblgd")
p <- plot_de_volcano(x) +
  ggtitle("PB DE genes (PB-HGD vs PB-LGD)")
p
f <- file.path(plot_dir, "volcano_de_pbhgd_vs_lgd.pdf")
ggsave(f, plot=p, device="pdf", width=5, height=5)

x <- filter(grade_de, analysis == "path_inthgd_vs_intlgd")
p <- plot_de_volcano(x) +
  ggtitle("INT DE genes (INT-HGD vs INT-LGD)")
p
f <- file.path(plot_dir, "volcano_de_inthgd_vs_intlgd.pdf")
ggsave(f, plot=p, device="pdf", width=5, height=5)


##################################################
#
#
# 2D plot of DE genes
#
#
##################################################

#
# Pathology plot
#
x <- grade_de_merged
nrow(grade_de_merged %>% filter(subtype != "none"))
table(grade_de_merged$subtype)

color_scale = c(none = "#888888",
                hgd = "#FF0000",
                lgd = "#0000FF",
                pbhgd = "#32cd32",
                inthgd = "#670092",
                pblgd = "#1677c4",
                intlgd = "#ff6700")

size_scale = c(none = 1,
               hgd = 2,
               lgd = 2,
               pbhgd = 2,
               pblgd = 2,
               inthgd = 2,
               intlgd = 2)

annot_data <- x %>%
  filter(subtype != "none")
none_data <- x %>%
  filter(subtype == "none")

label_data <- filter(x, subtype != "none")
label_data <- bind_rows(filter(label_data, subtype == "hgd") %>% slice_max(log2fc_path_hgd_vs_lgd, n = 30),
                        filter(label_data, subtype == "lgd") %>% slice_min(log2fc_path_hgd_vs_lgd, n = 10),
                        filter(label_data, subtype == "pbhgd") %>% slice_max(log2fc_path_pbhgd_vs_pblgd, n = 10),
                        filter(label_data, subtype == "pblgd") %>% slice_min(log2fc_path_pbhgd_vs_pblgd, n = 10),
                        filter(label_data, subtype == "inthgd") %>% slice_max(log2fc_path_inthgd_vs_intlgd, n = 10),
                        filter(label_data, subtype == "intlgd") %>% slice_min(log2fc_path_inthgd_vs_intlgd, n = 10))
label_data <- distinct(label_data)
nrow(label_data)

p <- ggplot(x, aes(x=log2fc_path_pbhgd_vs_pblgd, 
                   y=log2fc_path_inthgd_vs_intlgd,
                   color=subtype, 
                   size=subtype)) +
  geom_point(data=none_data, alpha=0.4) + 
  geom_point(data=annot_data, alpha=0.8) +
  geom_text_repel(data=label_data, color="black", size=3, aes(label=gene), max.overlaps=Inf) +
  scale_color_manual(values = color_scale) +
  scale_size_manual(values = size_scale, guide = "none") +
  theme_minimal() +
  xlab("log2(Fold Change) PB-HGD vs PB-LGD") +
  ylab("log2(Fold Change) INT-HGD vs INT-LGD") +
  labs(color = "Pathology")
p
f <- file.path(plot_dir, "scatter_de_grade.pdf")
ggsave(f, plot=p, device="pdf", width=8, height=8)


#
# Histology Plot
#
x <- hist_de_merged

color_scale = c(none = "#888888",
                pb = "#32cd32",
                int = "#670092",
                gf = "#1677c4",
                pbgf = "#ff6700",
                pbint = "#ff0000")

size_scale = c(none = 1,
               pb = 2,
               int = 2,
               gf = 2,
               pbgf = 2,
               pbint = 2)

annot_data <- x %>%
  filter(subtype != "none")
none_data <- x %>%
  filter(subtype == "none")

label_data <- bind_rows(filter(x, subtype == "pb") %>% slice_max(log2fc_hist_pb_vs_gf, n = 10),
                        filter(x, subtype == "pb") %>% slice_max(log2fc_hist_pb_vs_int, n = 10),
                        filter(x, subtype == "gf") %>% slice_min(log2fc_hist_pb_vs_gf, n = 10),
                        filter(x, subtype == "gf") %>% slice_min(log2fc_hist_int_vs_gf, n = 10),
                        filter(x, subtype == "int") %>% slice_max(log2fc_hist_int_vs_gf, n = 10),
                        filter(x, subtype == "int") %>% slice_min(log2fc_hist_pb_vs_int, n = 10),
                        filter(x, subtype == "pbgf") %>% slice_max(log2fc_hist_pb_vs_int, n = 10),
                        filter(x, subtype == "pbgf") %>% slice_min(log2fc_hist_int_vs_gf, n = 10),
                        filter(x, subtype == "pbint") %>% slice_max(log2fc_hist_pb_vs_gf, n = 10),
                        filter(x, subtype == "pbint") %>% slice_max(log2fc_hist_int_vs_gf, n = 10))

# manually add MUC1 as this is featured in plot but not in top 10
label_data <- bind_rows(label_data, filter(x, gene == "MUC1"))
label_data <- bind_rows(label_data, filter(x, gene == "MUC4"))
label_data <- distinct(label_data)
nrow(label_data)

p <- ggplot(x, aes(x=log2fc_hist_pb_vs_gf, 
                   y=log2fc_hist_int,
                   color=subtype, 
                   size=subtype)) +
  geom_point(data=none_data, alpha=0.4) + 
  geom_point(data=annot_data, alpha=0.8) +
  geom_text_repel(data=label_data, color="black", size=3, aes(label=gene), max.overlaps=Inf) +
  scale_color_manual(values = color_scale) +
  scale_size_manual(values = size_scale, guide = "none") +
  theme_minimal() +
  xlab("log2(Fold Change) PB vs GF") +
  ylab("log2(Fold Change) INT vs PB-GF") +
  labs(color = "Subtype")
p
f <- file.path(plot_dir, "scatter_de_histology.pdf")
ggsave(f, plot=p, device="pdf", width=8, height=8)




#
# PB plot
#
p <- ggplot(x, aes(x=log2fc_hist_pb_vs_gf, 
                   y=log2fc_hist_pb_vs_int,
                   color=subtype, 
                   size=subtype)) +
  geom_point(data=none_data, alpha=0.4) + 
  geom_point(data=annot_data, alpha=0.8) +
  geom_text_repel(data=label_data, color="black", size=3, aes(label=gene), max.overlaps=Inf) +
  scale_color_manual(values = color_scale) +
  scale_size_manual(values = size_scale, guide = "none") +
  theme_minimal() +
  xlab("log2(Fold Change) PB vs GF") +
  ylab("log2(Fold Change) PB vs INT") +
  labs(color = "Subtype")
p
f <- file.path(plot_dir, "scatter_de_pb.pdf")
ggsave(f, plot=p, device="pdf", width=8, height=8)


##################################################
#
#
# Plots of individual DE genes
#
#
##################################################

plot_hist_boxplot <- function(s, m, gene_symbol) {
  exprs <- unlist(m[gene_symbol,])
  x <- s %>%
    select(aoi_id, slide_id, segment, histology, path, histpath) %>%
    add_column(gene = exprs)
  p <- ggplot(x, aes(x=histology, y=gene)) + 
#  p <- ggplot(x, aes(x=reorder(histology, gene), y=gene)) + 
    geom_boxplot(aes(fill = histology), outlier.shape=NA) +
    geom_jitter(width=0.1) +
    scale_fill_manual(values=color_scales$histology) +
    xlab("Histology") +
    ylab("log2 Normalized CPM") + 
    labs(fill = "Histology") +
    ggtitle(paste0("Gene: ", gene_symbol)) +
    theme_minimal() + 
    theme(axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.line = element_line(color = "black"))
  return(p)
}

plot_grade_boxplot <- function(s, m, gene_symbol) {
  exprs <- unlist(m[gene_symbol,])
  x <- s %>%
    select(aoi_id, slide_id, segment, histology, path, histpath) %>%
    add_column(gene = exprs)
  p <- ggplot(x, aes(x=histpath, y=gene)) + 
    geom_boxplot(aes(fill = histpath), width=0.5, outlier.shape=NA) +
    geom_jitter(width=0.1) +
    scale_fill_manual(values=color_scales$histpath) +
    xlab("Histopathology") +
    ylab("log2 Normalized CPM") + 
    labs(fill = "Histopath") +
    ggtitle(gene_symbol) +
    theme_minimal() + 
    theme(axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.line = element_line(color = "black"))
  return(p)
}


plot_hgd_boxplot <- function(s, m, gene_symbol) {
  exprs <- unlist(m[gene_symbol,])
  x <- s %>%
    select(aoi_id, slide_id, segment, histology, path, histpath) %>%
    add_column(gene = exprs)
  p <- ggplot(x, aes(x=path, y=gene)) + 
    geom_boxplot(aes(fill = path), outlier.shape=NA, width=0.5) +
    geom_jitter(width=0.1) +
    scale_fill_manual(values=color_scales$path) +
    xlab("Grade") +
    ylab("log2 Normalized CPM") + 
    labs(fill = "Grade") +
    ggtitle(gene_symbol) +
    theme_minimal() + 
    theme(axis.line = element_line(color = "black")) +
    theme(axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))
  return(p)
}


# matrix for plotting
x <- norm_count_lcpm

# plot all de genes in heatmaps
for (g in union(hist_genes, grade_genes)) {
  p <- plot_hist_boxplot(samples, x, g)
  f <- file.path(gene_plot_dir, paste0("boxplot_", g, "_hist.pdf"))
  ggsave(f, p, width=3, height=4)
  p <- plot_grade_boxplot(samples, x, g)
  f <- file.path(gene_plot_dir, paste0("boxplot_", g, "_grade.pdf"))
  ggsave(f, p, width=3, height=4)
}

# manually plot individual genes
p <- plot_hist_boxplot(samples, x, "MUC1")
ggsave(file.path(plot_dir, "hist_boxplot_muc1.pdf"), p, width=3, height=5)
p <- plot_hist_boxplot(samples, x, "MUC4")
ggsave(file.path(plot_dir, "hist_boxplot_muc4.pdf"), p, width=3, height=5)
p <- plot_hist_boxplot(samples, x, "CLU")
ggsave(file.path(plot_dir, "hist_boxplot_clu.pdf"), p, width=3, height=5)
p <- plot_hist_boxplot(samples, x, "RBP4")
ggsave(file.path(plot_dir, "hist_boxplot_rbp4.pdf"), p, width=3, height=5)
p <- plot_hist_boxplot(samples, x, "KRT17")
ggsave(file.path(plot_dir, "hist_boxplot_krt17.pdf"), p, width=3, height=5)


p <- plot_hgd_boxplot(samples, x, "FOSL1")
ggsave(file.path(plot_dir, "hgd_boxplot_fosl1.pdf"), p, width=3, height=5)
p <- plot_hgd_boxplot(samples, x, "AREG")
ggsave(file.path(plot_dir, "hgd_boxplot_areg.pdf"), p, width=3, height=5)
p <- plot_hgd_boxplot(samples, x, "CCL20")
ggsave(file.path(plot_dir, "hgd_boxplot_ccl20.pdf"), p, width=3, height=5)
p <- plot_hgd_boxplot(samples, x, "SERPINB5")
ggsave(file.path(plot_dir, "hgd_boxplot_serpinb5.pdf"), p, width=3, height=5)
p <- plot_hgd_boxplot(samples, x, "ITGA2")
ggsave(file.path(plot_dir, "hgd_boxplot_itga2.pdf"), p, width=3, height=5)
p <- plot_hgd_boxplot(samples, x, "LAMC2")
ggsave(file.path(plot_dir, "hgd_boxplot_lamc2.pdf"), p, width=3, height=5)

p <- plot_grade_boxplot(samples, x, "SMAD4")
ggsave(file.path(plot_dir, "grade_boxplot_smad4.pdf"), p, width=3, height=4)

plot_grade_boxplot(samples, x, "SMAD4")
plot_grade_boxplot(samples, x, "CEACAM6")
plot_grade_boxplot(samples, x, "CD55")
grade_de_merged %>% filter(gene == "CD55") %>% as.data.frame()

# lymphocyte marker
plot_grade_boxplot(samples, x, "PTPRC")
plot_grade_boxplot(samples, x, "CD3D")
plot_grade_boxplot(samples, x, "IL7R")

# EMT-like cells
plot_grade_boxplot(samples, x, "ZEB1")

# ductal cells
plot_grade_boxplot(samples, x, "TACSTD2")
plot_grade_boxplot(samples, x, "KRT7")
# macrophages
plot_grade_boxplot(samples, x, "C1QA")
plot_grade_boxplot(samples, x, "C1QB")

# myeloid-derived suppressor cells
plot_grade_boxplot(samples, x, "S100A9")
plot_grade_boxplot(samples, x, "CCL3")

# KDR and VWF for endothelial cells; and FCER1A and CD1 for dendritic cells
# DC2 dendritic cells
plot_grade_boxplot(samples, x, "THBD")

# fibroblast marker
plot_grade_boxplot(samples, x, "ACTA2")
plot_grade_boxplot(samples, x, "COL1A1")

# inflammatory CAF
plot_grade_boxplot(samples, x, "COL3A1")
plot_grade_boxplot(samples, x, "CXCL12")
plot_grade_boxplot(samples, x, "TTR")

# proliferation markers
plot_grade_boxplot(samples, x, "MKI67")
plot_grade_boxplot(samples, x, "TOP2A")
plot_grade_boxplot(samples, x, "CENPF")
plot_grade_boxplot(samples, x, "UBE2C")

plot_grade_boxplot(samples, x, "KRT17")
plot_grade_boxplot(samples, x, "KRT6B")
plot_grade_boxplot(samples, x, "CXCL8")
plot_grade_boxplot(samples, x, "ITGA2")
plot_grade_boxplot(samples, x, "LAMC2")
plot_grade_boxplot(samples, x, "LAMB3")
plot_grade_boxplot(samples, x, "EPHA2")
plot_grade_boxplot(samples, x, "AREG")

plot_grade_boxplot(samples, x, "LIF")
plot_grade_boxplot(samples, x, "CXCL5")
plot_grade_boxplot(samples, x, "CLU")
plot_grade_boxplot(samples, x, "RBP4")
plot_grade_boxplot(samples, x, "SMAD4")
plot_grade_boxplot(samples, x, "CEACAM6")
plot_grade_boxplot(samples, x, "CXCL14")
plot_grade_boxplot(samples, x, "IL33")
plot_grade_boxplot(samples, x, "HMGA2")
plot_grade_boxplot(samples, x, "CD55")


##################################################
#
#
# GSEA Analysis
#
#
##################################################

get_de_ranks <- function(de, a) {
  ranks <- de %>% 
    filter(analysis == a) %>%
    select(gene, log2fc, padj) %>%
    mutate(rank = log2fc)
  ranks = sort(setNames(ranks$rank, ranks$gene), decreasing = TRUE)
  return(ranks)
}

plot_gsea_enrichment <- function(x, ranks, gs) {
  txt_stat <- paste0("ES=", round(x$ES,2), " NES=", round(x$NES, 2), " padj=", format(x$padj, scientific=TRUE, digits=3))
  txt_title <- paste0("ranks: ", a, " gs: ", x$pathway)
  p <- plotEnrichment(gs[[x$pathway]], ranks) +
    annotate(geom="text", x=150, y=0.1, label=txt_stat, hjust=0) +
    labs(title = txt_title, 
         xlab = "Rank",
         ylab = "Enrichment Score") +
    theme(plot.title = element_text(size=8))
  return(p)
}

run_batch_fgsea <- function(my_analyses, my_de, my_gs, my_prefix, my_plot_dir, my_padj_cutoff=0.01) {
  gsea <- NULL
  for (a in my_analyses) {
    ranks <- get_de_ranks(my_de, a)
    res <- fgsea(pathways = my_gs, stats = ranks, minSize = 10, eps = 0, nPermSimple = 10000)
    res <- mutate(res, analysis = a)
    gsea <- bind_rows(gsea, res)
    print(a)
    for (i in 1:nrow(res)) {
      x <- res[i,]
      if (x$padj >= my_padj_cutoff) { next }
      print(x$pathway)
      p <- plot_gsea_enrichment(x, ranks, my_gs)
      f <- file.path(my_plot_dir, paste0(my_prefix, "_", a, "_gs_", x$pathway, ".pdf"))
      ggsave(f, plot=p, device="pdf", width=5, height=3)
    }
  }
  return(gsea)
}

gmt_to_tbl <- function(gmt) {
  gs_to_tbl <- function(gs_name, gs_genes) {
    bind_cols(gs_name = rep(gs_name, length(gs_genes)),
              gene_symbol = gs_genes)
  }
  tbl <- NULL
  for (i in 1:length(gmt)) {
    tbl <- bind_rows(tbl, gs_to_tbl(names(gmt)[i], gmt[[i]]))
  }
  return(tbl)
}

plot_gsea_volcano <- function(x, padj_cutoff = 0.01, max.overlaps = Inf) {
  p <- ggplot(x, aes(x=NES, y=-log10(padj), color=analysis)) + 
    geom_point(data=subset(x, padj > padj_cutoff), size=1, alpha=0.4) + 
    geom_point(data=subset(x, padj <= padj_cutoff), size=2, alpha=0.8) +
    geom_text_repel(data=subset(x, padj <= padj_cutoff), color="black", size=3, aes(label=pathway), max.overlaps=max.overlaps) +
    ylab("-log10(adjusted p-value)") +
    xlab("NES") +
    theme_minimal() +
    theme(legend.position="bottom") + 
    theme(axis.line = element_line(color = "black")) +
    coord_flip() 
  return(p)
}

plot_gsea_barplot <- function(x, padj_cutoff = 0.01) {
  x <- mutate(x, sig = ifelse(padj < padj_cutoff, ifelse(NES < 0, "dn", "up"), "no"),
              neglog10padj = ifelse(padj < padj_cutoff, -log10(padj), NA))
  p <- ggplot(x, aes(x=reorder(pathway, NES), y=NES, fill=neglog10padj)) +
    geom_col() +
    scale_fill_gradient(low = "blue", high = "cyan", na.value="#aaaaaa") +
    coord_flip() +
    labs(x="Analysis", y="Normalized Enrichment Score") +
    theme_minimal()
  return(p)
}

gsea_leading_edge_matrix <- function(x, my_padj_cutoff) {
  # filter for significant results
  sig_pathways <- x %>%
    filter(padj < my_padj_cutoff) %>% select(pathway) %>% distinct() %>% pull()
  x <- x %>%
    filter(pathway %in% sig_pathways) %>%
    arrange(desc(NES))
  
  # get leading edge genes
  leading_edge_genes <- summarise(x, result = unlist(leadingEdge, use.names=FALSE)) %>% distinct() %>% pull()
  
  # matrix genes (rows) vs pathways (columns)
  m <- matrix(data = 0, nrow=length(leading_edge_genes), ncol=nrow(x))
  m <- as.data.frame(m)
  rownames(m) <- leading_edge_genes
  colnames(m) <- x$pathway
  for (i in 1:nrow(x)) {
    p <- x$pathway[i]
    nes <- x$NES[i]
    padj <- x$padj[i]
    for (j in x$leadingEdge[i]) {
      m[j, p] <- sign(nes)
    }
  }
  # filter genes only found in a single pathway
  m <- m[abs(rowSums(m)) >= 2, ]
  # filter empty columns
  m <- m[, abs(colSums(m)) > 1]
  # order by number of pathways sharing the gene
  m <- m[order(rowSums(m), decreasing=TRUE),]
  m <- m[,order(colSums(m), decreasing=TRUE)]
  return(m)
}


#
# Setup Gene Sets
# 
# cptac gene sets
cptac_gs <- get_cptac_gene_sets(cptac_de_data)
lapply(cptac_gs, length)

# hpa gene sets
hpa_gs <- get_hpa_gene_sets(hpa)
hpa_panc_gs <- hpa_gs[c("HPA_PANCREATIC_CANCER_UNFAV_PROGNOSIS", 
                        "HPA_PANCREATIC_CANCER_FAV_PROGNOSIS")]

# read msigdb data
msigdb_gene_sets = msigdbr(species = "Homo sapiens")
# filter CTA genes
msigdb_gene_sets <- filter(msigdb_gene_sets, gene_symbol %in% cta_genes) 

# filter gruetzmann microarray dataset
gruetzmann_gs_names <- c("GRUETZMANN_PANCREATIC_CANCER_UP",
                         "GRUETZMANN_PANCREATIC_CANCER_DN")
gruetzmann_gs <- msigdb_gene_sets %>%
  filter(gs_name %in% gruetzmann_gs_names)
gruetzmann_gs <- split(x = gruetzmann_gs$gene_symbol,
                       f = gruetzmann_gs$gs_name)

# PDAC enrichment gene sets
gs_pdac <- gruetzmann_gs
gs_pdac <- append(gs_pdac, oncotarget_gs)
gs_pdac <- append(gs_pdac, hpa_panc_gs)
gs_pdac <- append(gs_pdac, cptac_gs)
lapply(gs_pdac, length)

# discovery analysis across msigdb gene sets
gs_hallmark <- filter(msigdb_gene_sets, gs_cat == "H")
gs_hallmark <- split(x = gs_hallmark$gene_symbol, f = gs_hallmark$gs_name)
gs_gobp <- filter(msigdb_gene_sets, gs_cat == "C5", gs_subcat == "GO:BP")
gs_gobp <- split(x = gs_gobp$gene_symbol, f = gs_gobp$gs_name)

#
# Setup analysis
#
de <- bind_rows(hist_de, grade_de)
analyses <- unique(de$analysis)

#
# PDAC GSEA Analysis
#
prefix <- "gsea_pdac"
padj_cutoff <- 0.05
gsea <- run_batch_fgsea(analyses, de, gs_pdac, "gsea_pdac", gsea_plot_dir, padj_cutoff)
gsea_pdac <- as_tibble(gsea)

# write gsea results
f <- file.path(working_dir, paste0(prefix, "_results.tsv"))
data.table::fwrite(gsea, file=f, sep="\t", sep2=c("", " ", ""))
f <- file.path(working_dir, paste0(prefix, "_results.xlsx"))
write_xlsx(gsea, f)


#
# Hallmark gene set analysis
#
prefix <- "gsea_hallmark"
padj_cutoff <- 0.01
gsea <- run_batch_fgsea(analyses, de, gs_hallmark, prefix, gsea_plot_dir, padj_cutoff)
gsea_hallmark <- as_tibble(gsea) %>%
  left_join(select(msigdb_gene_sets, gs_cat, gs_subcat, gs_name) %>% distinct(), 
            by=c("pathway"="gs_name"))
# write gsea results
f <- file.path(working_dir, paste0(prefix, "_results.tsv"))
data.table::fwrite(gsea, file=f, sep="\t", sep2=c("", " ", ""))
f <- file.path(working_dir, paste0(prefix, "_results.xlsx"))
write_xlsx(gsea, f)

#
# Gene Ontology analysis
#
prefix <- "gsea_gobp"
padj_cutoff <- 0.01
gsea <- run_batch_fgsea(analyses, de, gs_gobp, prefix, gsea_plot_dir, padj_cutoff)
gsea_gobp <- as_tibble(gsea) %>%
  left_join(select(msigdb_gene_sets, gs_cat, gs_subcat, gs_name) %>% distinct(), 
            by=c("pathway"="gs_name"))
# write gsea results
f <- file.path(working_dir, paste0(prefix, "_results.tsv"))
data.table::fwrite(gsea, file=f, sep="\t", sep2=c("", " ", ""))
f <- file.path(working_dir, paste0(prefix, "_results.xlsx"))
write_xlsx(gsea, f)


#
# HGD vs LGD gsea plots
#
x <- filter(gsea_hallmark, analysis == "path_hgd_vs_lgd")
x$pathway <- gsub("HALLMARK_", "", x$pathway)
p <- plot_gsea_barplot(x, padj_cutoff = 0.01)
ggsave(file.path(plot_dir, "gsea_hgd_vs_lgd_hallmark_barplot.pdf"), plot=p, width=6, height=6)

x <- filter(gsea_gobp, analysis == "path_hgd_vs_lgd", padj < 0.01)
x$pathway <- gsub("GOBP_", "", x$pathway)
p <- plot_gsea_barplot(x, padj_cutoff = 0.01)

x <- filter(gsea_hallmark, analysis == "path_hgd_vs_lgd")
x$pathway <- gsub("HALLMARK_", "", x$pathway)
p <- plot_gsea_volcano(x, padj_cutoff = 0.01) + 
  ggtitle("MSigDB Hallmark Gene Sets")
p
ggsave(file.path(plot_dir, "gsea_hgd_vs_lgd_hallmark_volcano.pdf"), plot=p, width=5, height=3)

x <- filter(gsea_gobp, analysis == "path_hgd_vs_lgd")
x$pathway <- gsub("GOBP_", "", x$pathway)
p <- plot_gsea_volcano(x, padj_cutoff = 0.01, max.overlaps=Inf) + 
  ggtitle("GO:BP Gene Sets")
p
ggsave(file.path(plot_dir, "gsea_hgd_vs_lgd_gobp_volcano.pdf"), plot=p, width=8, height=3)


#
# PDAC HGD vs LGD plots
#
padj_cutoff <- 0.05
x <- filter(gsea_pdac, analysis == "path_hgd_vs_lgd")
x <- mutate(x, sig = ifelse(padj < padj_cutoff, ifelse(NES < 0, "dn", "up"), "no"),
            neglog10padj = ifelse(padj < padj_cutoff, -log10(padj), NA))
p <- ggplot(x, aes(x=reorder(pathway, NES), y=NES, fill=neglog10padj)) +
  geom_col() +
  scale_fill_gradient(low = "#004488", high = "cyan", na.value="#aaaaaa") +
  labs(x="Analysis", y="Normalized Enrichment Score",
       title="PDAC Gene Sets") +
  coord_flip() +
  theme_minimal() + 
  theme(axis.text.y = element_text(hjust=1, size=4)) +
  theme(legend.position="bottom")
p
ggsave(file.path(plot_dir, "gsea_pdac_hgd_vs_lgd_barplots.pdf"), p, width=6, height=6)


#
# PDAC plots comparing analysis
#
unique(gsea_pdac$analysis)
padj_cutoff <- 0.05
plotanalyses <- c("path_intlgd_vs_pblgd", "path_inthgd_vs_pblgd", "path_inthgd_vs_intlgd",
                  "path_pbhgd_vs_inthgd", "path_pbhgd_vs_intlgd", "path_pbhgd_vs_pblgd",
                  "path_hgd_vs_lgd")
x <- filter(gsea_pdac, analysis %in% plotanalyses)
x$analysis <- factor(x$analysis, levels=plotanalyses)
x <- x %>% mutate( 
  sig = ifelse(padj < padj_cutoff, ifelse(NES < 0, "dn", "up"), "no"),
  neglog10padj = -log10(x$padj),
  score = NES * -log10(x$padj),
  fillneglog10padj = ifelse(padj < padj_cutoff, -log10(padj), NA),
  fillnes = ifelse(padj < padj_cutoff, NES, NA)
)

p <- ggplot(x, aes(x=analysis, y=NES, fill=fillneglog10padj)) +
  geom_bar(stat="identity", position=position_dodge()) +
#  coord_flip() +
#  scale_fill_viridis_c() +
  scale_fill_gradient(low = "#004488", high = "cyan", na.value="#aaaaaa") +
  labs(x="Analysis", y="Normalized Enrichment Score",
       title="PDAC Gene Signatures") + 
  theme_minimal() +
  theme(axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1, size=6)) +
  facet_wrap(~ pathway, nrow=1)
p
ggsave(file.path(plot_dir, "gsea_pdac_subtype_barplots.pdf"), p, width=9, height=5)


#
# Hallmark analysis by subtype
#
padj_cutoff <- 0.01

x <- filter(gsea_hallmark, analysis %in% plotanalyses)
x$analysis <- factor(x$analysis, levels=plotanalyses)
x$abbrev <- gsub("HALLMARK_", "", x$pathway)
p <- ggplot(x, aes(x=NES, y=-log10(padj), color=analysis)) + 
  geom_point(data=subset(x, padj > padj_cutoff), size=1, alpha=0.4) + 
  geom_point(data=subset(x, padj <= padj_cutoff), size=2, alpha=0.8) +
  geom_text_repel(data=subset(x, padj <= padj_cutoff), color="black", size=2, aes(label=abbrev), max.overlaps=Inf) +
  ylab("-log10(adjusted p-value)") +
  xlab("NES") +
  ggtitle("MSigDB Hallmark Gene Sets") +
  theme_minimal() +
  theme(legend.position="bottom") + 
  theme(axis.line = element_line(color = "black")) +
  coord_flip() + 
  facet_wrap(~ analysis, ncol = 1)
p
ggsave(file.path(plot_dir, "gsea_hallmark_subtype_volcano.pdf"), plot=p, width=4, height=10)


# GO
x <- filter(gsea_gobp, analysis %in% plotanalyses)
x$analysis <- factor(x$analysis, levels=plotanalyses)
x$abbrev <- gsub("GOBP_", "", x$pathway)
p <- ggplot(x, aes(x=NES, y=-log10(padj), color=analysis)) + 
  geom_point(data=subset(x, padj > padj_cutoff), size=1, alpha=0.4) + 
  geom_point(data=subset(x, padj <= padj_cutoff), size=2, alpha=0.8) +
  geom_text_repel(data=subset(x, padj <= padj_cutoff), color="black", size=2, aes(label=abbrev), max.overlaps=Inf) +
  ylab("-log10(adjusted p-value)") +
  xlab("NES") +
  ggtitle("GO:BP Gene Sets") +
  theme_minimal() +
  theme(legend.position="bottom") + 
  theme(axis.line = element_line(color = "black")) +
  coord_flip() + 
  facet_wrap(~ analysis, ncol = 1)
p
ggsave(file.path(plot_dir, "gsea_gobp_subtype_volcano.pdf"), plot=p, width=10, height=15)

#
# GSEA Leading Edge
#
padj_cutoff <- 0.01
gsea_merged <- gsea_hallmark
x <- gsea_merged %>% filter(analysis == "path_hgd_vs_lgd")
m <- gsea_leading_edge_matrix(x, padj_cutoff)

# annotate rows
row_annot <- gene_annot %>%
  filter(gene_symbol %in% rownames(m)) %>%
  select(gene_symbol, grade_subtype, prognostic, de_cptac_rna, de_cptac_prot) %>%
  column_to_rownames("gene_symbol") %>%
  as.data.frame()

p <- pheatmap(m, 
              scale = "none",
              color = colorRampPalette(c("blue", "white", "red"))(50),
              breaks = seq(-1, 1, length.out=51),
              cellwidth = 6,
              cellheight = 6,
              border_color = "black",
              cluster_rows = TRUE,
              cluster_cols = TRUE,
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              clustering_method = "ward.D2",
              annotation_row = row_annot,
              annotation_colors = color_scales,
              fontsize_col = 4,
              fontsize_row = 6)
h <- 2 + nrow(m) / 6
w <- 2 + ncol(m) / 2
save_pheatmap_pdf(p, filename=file.path(plot_dir, "gsea_leading_edge_hallmark.pdf"), width=w, height=h)


# GOBP leading edge
gsea_merged <- gsea_gobp
x <- gsea_merged %>% filter(analysis == "path_hgd_vs_lgd")
m <- gsea_leading_edge_matrix(x, padj_cutoff)

# annotate rows
row_annot <- gene_annot %>%
  filter(gene_symbol %in% rownames(m)) %>%
  select(gene_symbol, grade_subtype, prognostic, de_cptac_rna, de_cptac_prot) %>%
  column_to_rownames("gene_symbol") %>%
  as.data.frame()

p <- pheatmap(m, 
              scale = "none",
              color = colorRampPalette(c("blue", "white", "red"))(50),
              breaks = seq(-1, 1, length.out=51),
              cellwidth = 6,
              cellheight = 6,
              border_color = "black",
              cluster_rows = TRUE,
              cluster_cols = TRUE,
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              clustering_method = "ward.D2",
              annotation_row = row_annot,
              annotation_colors = color_scales,
              fontsize_col = 4,
              fontsize_row = 6)
h <- 2 + nrow(m) / 6
w <- 2 + ncol(m) / 2
save_pheatmap_pdf(p, filename=file.path(plot_dir, "gsea_leading_edge_gobp.pdf"), width=w, height=h)


# GOBP leading edge
gsea_merged <- bind_rows(gsea_hallmark, gsea_gobp)
x <- gsea_merged %>% filter(analysis == "path_hgd_vs_lgd")
m <- gsea_leading_edge_matrix(x, padj_cutoff)

# annotate rows
row_annot <- gene_annot %>%
  filter(gene_symbol %in% rownames(m)) %>%
  select(gene_symbol, grade_subtype, prognostic, de_cptac_rna, de_cptac_prot) %>%
  column_to_rownames("gene_symbol") %>%
  as.data.frame()

p <- pheatmap(m, 
              scale = "none",
              color = colorRampPalette(c("blue", "white", "red"))(50),
              breaks = seq(-1, 1, length.out=51),
              cellwidth = 6,
              cellheight = 6,
              border_color = "black",
              cluster_rows = TRUE,
              cluster_cols = TRUE,
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              clustering_method = "ward.D2",
              annotation_row = row_annot,
              annotation_colors = color_scales,
              fontsize_col = 4,
              fontsize_row = 6)
h <- 2 + nrow(m) / 6
w <- 2 + ncol(m) / 2
save_pheatmap_pdf(p, filename=file.path(plot_dir, "gsea_leading_edge_merged.pdf"), width=w, height=h)


#############################################################
#
# correlation analysis
#
#############################################################

cor_compute <- function(x, y, nperms, cor_method = "spearman") {
  # setup
  ngenes_x <- nrow(x)
  ngenes_y <- nrow(y)
  nsamples <- ncol(x)
  pvals <- array(data = 0, dim = c(ngenes_x, ngenes_y), dimnames = list(rownames(x), rownames(y)))

  # convention is rows are genes, columns are samples but 
  # we need to transpose this for correlation analysis
  x <- t(x)
  y <- t(y)
  r <- cor(x, y, method = cor_method)

  # generate null correlation distribution
  for (i in 1:nperms) {
    r_i <- cor(x, y[sample(nsamples), ], method = cor_method)
    pvals <- pvals + (r_i >= r)
    if (i %% 100 == 0) {
      print(i)
      flush.console()
    }
  }
  pvals <- pvals / nperms
  pvals <- 2 * pmin(pvals, 1 - pvals)
  return(list(r = r, pvals = pvals))
}


pivot_matrix <- function(m, 
                         row.name = "source", 
                         col.name = "target", 
                         value.name = "weight") {
  tbl <- as_tibble(m, rownames=row.name) %>%
    pivot_longer(cols=colnames(m),
                 names_to=col.name,
                 values_to=value.name,
                 values_drop_na=TRUE)
  return(tbl)
}


cor_to_edgelist <- function(corobj, 
                            row.name = "source", 
                            col.name = "target", 
                            self = FALSE) {
  r <- corobj$r
  if (self) {
    r[lower.tri(r, diag=TRUE)] <- NA
  }
  redges <- pivot_matrix(r, row.name, col.name, value.name = "cor")
  
  p <- corobj$pvals
  if (self) {
    p[lower.tri(p, diag=TRUE)] <- NA
  }
  pedges <- pivot_matrix(p, row.name, col.name, value.name = "pval")
  edges <- bind_cols(redges, pval = pedges$pval)
  return(edges)
}


subset_paired_segment <- function(meta, normexpr, segment_name) {
  subset_meta <- meta %>% filter(segment == segment_name)
  subset_mat <- normexpr[, subset_meta$aoi_id]
  colnames(subset_mat) <- subset_meta$pair_id
  genes <- rownames(subset_mat)
  rownames(subset_mat) <- paste0(segment_name, "_", genes)
  return(list(meta = subset_meta,
              mat = subset_mat,
              genes = genes))
}


cor_result_to_edges <- function(cor_result, self=FALSE) {
  edges <- cor_to_edgelist(cor_result, self=self)
  edges$padj <- p.adjust(edges$pval, method = "fdr")
  qobj <- qvalue(edges$pval)
  edges$qval <- qobj$qvalues
  return (edges)
}


#
# correlation analysis
#
nperms <- 10000
cor_qval_cutoff <- 0.01
cor_r_cutoff <- 0.7
cor_weight_beta <- 4
graph_min_csize <- 10
graph_clust_resolution <- 0.50

# panck self correlation
panck_cor_mat <- norm_count_lcpm
rownames(panck_cor_mat) <- paste0("panck_", rownames(norm_count_lcpm))
# panck_cor_res <- cor_compute(panck_cor_mat,
#                              panck_cor_mat,
#                              nperms)
# saveRDS(panck_cor_res, file=file.path(working_dir, "cor_panck.rds"))
panck_cor_res <- readRDS(file.path(working_dir, "cor_panck.rds"))

#
# panck network
#
# edges
edges_all <- cor_result_to_edges(panck_cor_res, self=TRUE)
# filter significant edges
edges <- filter(edges_all, qval < cor_qval_cutoff,
                cor >= cor_r_cutoff)
edges$weight <- abs(edges$cor) ** cor_weight_beta
#edges$weight <- (1 + edges$cor)/2 ** cor_weight_beta
#edges$cordir <- ifelse(edges$cor < 0, "neg", "pos")

# nodes
nodes <- bind_cols(label = rownames(panck_cor_mat),
                   orig_gene = rownames(norm_count_lcpm),
                   type = "panck",
                   grade_subtype = gene_annot$grade_subtype,
                   grade_log2fc = grade_de_merged$log2fc_path_hgd_vs_lgd,
                   grade_pb_log2fc = grade_de_merged$log2fc_path_pbhgd_vs_pblgd,
                   grade_int_log2fc = grade_de_merged$log2fc_path_inthgd_vs_intlgd,
                   grade_abs_log2fc = abs(grade_de_merged$log2fc_path_hgd_vs_lgd),
                   grade_abs_pb_log2fc = abs(grade_de_merged$log2fc_path_pbhgd_vs_pblgd),
                   grade_abs_int_log2fc = abs(grade_de_merged$log2fc_path_inthgd_vs_intlgd),
                   inthgd_scaled_mean = gene_annot$inthgd_scaled_mean,
                   pbhgd_scaled_mean = gene_annot$pbhgd_scaled_mean)
nodes$id <- nodes$label
nodes$name <- nodes$id
nodes <- mutate(nodes, 
                cptac_de = case_when(orig_gene %in% cptac_gs$CPTAC_RNA_UP ~ "up",
                                     orig_gene %in% cptac_gs$CPTAC_RNA_DN ~ "dn",
                                     TRUE ~ "no"))
# filter nodes
nodes <- nodes %>% filter(id %in% union(edges$source, edges$target))


#
# build graph object
#
g <- igraph::graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)

igraph::components(g)$csize


res <- NULL
i <- 1
for (r in seq(0.01, 1, by=0.01)) {
  clust <- cluster_leiden(g, 
                          objective_function="modularity",
                          resolution=r,
                          n_iterations = 100)
  res[[i]] <- tibble(m = modularity(g, membership(clust)),
                     n = clust$nb_clusters,
                     q = clust$quality,
                     r = r)
  i <- i+ 1
}
res <- bind_rows(res)
p1 <- ggplot(res, aes(x=r, y=m)) + 
  geom_point() + 
  geom_vline(xintercept=graph_clust_resolution, linetype="dashed", color="red") +
  labs(x="resolution", y="modularity")
p2 <- ggplot(res, aes(x=r, y=n)) + 
  geom_point() +
  geom_vline(xintercept=graph_clust_resolution, linetype="dashed", color="red") +
  labs(x="resolution", y="# of clusters")

library(patchwork)
p <- p1 / p2
f <- file.path(plot_dir, "clust_leiden_modularity.pdf")
ggsave(f, p)

# clustering/community detection
clust <- cluster_leiden(g,
                        objective_function="modularity",
                        resolution=graph_clust_resolution,
                        n_iterations=1000)
graph_clust_resolution
modularity(g, membership(clust))
table(membership(clust))
V(g)$community <- factor(clust$membership)

clust_nodes <- as_tibble(as_data_frame(g, what="vertices"))
clust_edges <- as_tibble(as_data_frame(g, what="edges"))
clust_nodes <- clust_nodes %>%
  dplyr::rename(label = name)
clust_edges <- clust_edges %>%
  dplyr::rename(source = from,
                target = to)

write.table(clust_nodes, file=file.path(working_dir, "graph_clust_nodes.tsv"),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
write.table(clust_edges, file=file.path(working_dir, "graph_clust_edges.tsv"),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

merged_nodes <- clust_nodes %>%
  dplyr::group_by(community) %>%
  dplyr::summarize(size = n()) %>%
  dplyr::mutate(id = paste0("c", community),
                label = id)

nclust <- length(unique(clust_nodes$community))
e <- filter(edges_all, qval < cor_qval_cutoff)

merged_edges <- NULL
for (i in 1:nclust) {
  ni <- clust_nodes %>% dplyr::filter(community == i) %>% dplyr::select(id) %>% dplyr::pull()
  if (length(ni) < graph_min_csize) {
    next
  }
  for (j in (i+1):nclust) {
    nj <- clust_nodes %>% dplyr::filter(community == j) %>% dplyr::select(id) %>% dplyr::pull()
    if (length(nj) < graph_min_csize) {
      next
    }
    eij <- bind_rows(e %>% filter(source %in% ni, target %in% nj),
                     e %>% filter(source %in% nj, target %in% ni))
    if (nrow(eij) == 0) {
      next
    }
    avgcor <- mean(eij$cor)
    tbl <- tibble(source = paste0("c", i), 
                  target = paste0("c", j), 
                  avgcor = avgcor,
                  n = nrow(eij))
    merged_edges <- bind_rows(merged_edges, tbl)
    print(paste0(i, " d ", j))
  }
}
summary(merged_edges$avgcor)
merged_edges <- merged_edges %>%
  mutate(weight = 1 - avgcor,
         cordir = ifelse(avgcor < 0, "neg", "pos"))
# filter nodes
merged_nodes <- merged_nodes %>% filter(id %in% union(merged_edges$source, merged_edges$target))
merged_nodes[,1] <- merged_nodes$id

write.table(merged_nodes, file=file.path(working_dir, "graph_merged_clust_nodes.tsv"),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
write.table(merged_edges, file=file.path(working_dir, "graph_merged_clust_edges.tsv"),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

mg <- igraph::graph_from_data_frame(d=merged_edges, vertices=merged_nodes, directed=FALSE)
plot(mg, 
     layout = layout_with_kk, 
     vertex.size = 20)

run_batch_enricher <- function(my_nodes, my_gs) {
  seg <- "panck"
  yy <- NULL
  for (i in 1:length(unique(my_nodes$community))) {
    x <- dplyr::filter(my_nodes, type == seg, community == i) %>% 
      dplyr::select(orig_gene) %>% pull()
    res <- enricher(gene = x,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    universe = cta_genes,
                    qvalueCutoff = 0.1,
                    TERM2GENE = my_gs)
    if (is.null(res)) { next }
    y <- res@result
    y$cluster_size <- length(res@gene)
    y$community <- i
    yy <- rbind(yy, y)
  }
  return(as_tibble(yy))
}


gs <- msigdb_gene_sets %>%
  filter(gs_cat == "H") %>%
  select(gs_name, gene_symbol)
clust_hallmark <- run_batch_enricher(clust_nodes, gs)

gs <- msigdb_gene_sets %>%
  dplyr::filter(gs_cat == "C5", gs_subcat == "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)
clust_gobp <- run_batch_enricher(clust_nodes, gs)

gs <- gmt_to_tbl(gs_pdac)
clust_pdac <- run_batch_enricher(clust_nodes, gs)
clust_pdac %>% arrange(p.adjust) %>% as.data.frame()

clust_pdac %>% filter(ID == "CPTAC_RNA_UP") %>% arrange(p.adjust) %>% as.data.frame()


gs <- msigdb_gene_sets %>%
  dplyr::select(gs_name, gene_symbol)
clust_all <- run_batch_enricher(clust_nodes, gs)

# combine and filter results
x <- bind_rows(clust_hallmark, clust_gobp, clust_pdac)
x <- x %>% filter(p.adjust < 0.05, cluster_size >= 10)
write_xlsx(
  list(community_enrichment = x),
  file.path(working_dir, community_enrichment_xlsx)
)

# inspect individual results
x <- clust_all
x <- x %>% filter(p.adjust < 0.05, cluster_size >= 10) %>% select(ID, GeneRatio, p.adjust, cluster_size, community)
x %>% as.data.frame()
x %>% filter(community == 5) %>% as.data.frame()

x <- clust_hallmark
x <- x %>% filter(p.adjust < 0.05, cluster_size >= 10) %>% select(ID, GeneRatio, p.adjust, Count, cluster_size, community)
# x %>% filter(community == 4) %>% arrange(p.adjust) %>% as.data.frame()
x <- x %>% as.data.frame()
x

x <- clust_gobp
x <- x %>% filter(p.adjust < 0.05, cluster_size >= 10) %>% select(ID, GeneRatio, p.adjust, Count, cluster_size, community)
# x %>% filter(community == 7) %>% arrange(p.adjust) %>% as.data.frame()
x %>% as.data.frame()

x <- clust_pdac
x
x %>% arrange(p.adjust) %>% as.data.frame()
x <- x %>% filter(p.adjust < 0.05, cluster_size >= 10) %>% select(ID, GeneRatio, p.adjust, Count, cluster_size, community)
x %>% as.data.frame()


clust_hallmark %>% arrange(p.adjust)

x <- bind_cols(clust_hallmark, oddsratio=sapply(clust_hallmark$GeneRatio, function(x) eval(parse(text=x))) / sapply(clust_hallmark$BgRatio, function(x) eval(parse(text=x))))
x %>% arrange(p.adjust)


###################################################################
#
# Address reviewer comments that DE genes may not be expressed at
# levels higher than signal background / limit of detection
#
###################################################################

s <- select(samples, aoi_id, slide_id, histology, path, histpath)
m <- select(counts_probe, -all_of(samples$aoi_id))
x <- select(counts_probe, all_of(samples$aoi_id))
g <- "KRT17"

x <- get_bgsub_qnorm_cpm(counts_probe, samples)
y <- x %>% 
  filter((probe_type == 0) | (gene_symbol == g)) %>%
  pivot_longer(samples$aoi_id, names_to="aoi_id", values_to="cpm") %>%
  inner_join(s, by=c("aoi_id"="aoi_id"), suffix=c("_gene", "_aoi")) %>%
  group_by(aoi_id, probe_type) %>%
  summarise(histology = first(histology),
            path = first(path),
            slide = first(slide_id),
            value = quantile(cpm, qc_negprobe_lod_quantile)) %>%
  ungroup()
y$probe_type <- if_else(y$probe_type == 0, "bg", "gene")

# x <- sweep(x, 2, colSums(x) / 1e6, FUN="/")
# y <- bind_cols(m, x)
# y <- y %>% 
#   filter((probe_type == 0) | (gene_symbol == g)) %>%
#   pivot_longer(colnames(x), names_to="aoi_id", values_to="cpm") %>%
#   inner_join(s, by=c("aoi_id"="aoi_id"), suffix=c("_gene", "_aoi")) %>%
#   group_by(aoi_id, probe_type) %>%
#   summarise(histology = first(histology),
#             path = first(path),
#             slide = first(slide_id),
#             value = quantile(cpm, qc_negprobe_lod_quantile)) %>%
#   ungroup()

p <- ggplot(y, aes(x=reorder(aoi_id, value), y=log2(value), color=factor(probe_type))) + 
  geom_point(size=2, alpha=0.8) +
#  scale_y_log10() +
  scale_color_manual(values = pals::cols25()) +
  labs(color = "LOD", x="AOI", y="log2 Normalized CPM") +
  theme_minimal() + 
  theme(axis.text.x = element_blank()) +
  # theme(axis.text.x = element_blank(),
  #       panel.grid.major.x = element_blank(),
  #       panel.grid.minor.x = element_blank()) +
  facet_grid(~ histology, scales="free")
p
f <- file.path(plot_dir, "qc_bg_detection_krt17.pdf")
ggsave(f, plot=p, device="pdf", width=6, height=3)

y <- y %>% pivot_wider(id_cols = c(aoi_id, slide, path, histology),
                       names_from = probe_type,
                       values_from = value)

table(y$histology, y$gene > y$bg)
table(y$gene > y$bg)
61 / (61+ 22)
median(log2(y$bg))

# p <- ggplot(y, aes(x=bg, y=gene, color=histology)) +
#   geom_point() + 
#   scale_color_manual(values=color_scales$histology) +
#   scale_x_log10() +
#   scale_y_log10() +
#   geom_abline(slope=1, color="red", linetype="dashed")


