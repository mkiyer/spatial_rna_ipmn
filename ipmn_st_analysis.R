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
input_xlsx_file = file.path(data_dir, "ipmn_counts_annot_2022-06-22.xlsx")
# spatial rna sheet names
metadata_sheet <- "SegmentProperties"
metadata_id_key <- "SegmentDisplayName"
count_data_sheet <- "BioProbeCountMatrix"

# external resource files
msigdb_symbol_file <- file.path(data_dir, "msigdb.v7.5.1.symbols.gmt")
hpa_path_file <- file.path(data_dir, "hpa_pathology.tsv")
cptac_panc_file <- file.path(data_dir, "cptac_panc_supp_mmc3.xlsx")
gs_wikipathways_file <- file.path(data_dir, "WikiPathways_2019_Human.txt")
gs_reactome_file <- file.path(data_dir, "Reactome_2016.txt")


# segments included
analysis_segments = c("panck", "cd45", "sma")

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
de_results_hist_xlsx <- "de_results_hist.xlsx"
de_results_grade_xlsx <- "de_results_grade.xlsx"
gsea_results_xlsx <- "gsea_results.xlsx"

# output directories
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
cor_clust_dir <- file.path(working_dir, "cor_clust")
if (!dir.exists(cor_clust_dir)) {
  dir.create(cor_clust_dir)
}

get_color_scales <- function(slide_ids) {
  color_scales = list(
    region = c(epithelia="#00203f", stroma="#adefd1"),
    path = c(lgd="#AAAAAA", hgd="#DD0000"),
    pathbw = c(lgd="#cccccc", hgd="#333333"),
    histology = c(pb = "#32cd32", intestinal = "#670092", gf="#1677c4"),
    histology2 = c(intestinal = "#fff01f", pb = "#1BF118", gf="#360167"),
    histpath = c(panck_intlgd="#f0ceff", panck_inthgd="#b042ff", 
                 panck_pblgd="#9df985", panck_pbhgd="#00ab08",
                 cd45_intlgd="#86c3fa", cd45_inthgd="#1750ac",
                 cd45_pblgd="#fff192", cd45_pbhgd="#ffd400",
                 sma_intlgd="#ff9090", sma_inthgd="#ec0000",
                 sma_pblgd="#ffddb3", sma_pbhgd="#ff8c00"),
    #histpath = c(panck_intlgd="#f0ceff", panck_inthgd="#c175ff", panck_pblgd="#00ffff", panck_pbhgd="#00bb00"),
    hgd_histology = c(pb = "#32cd32", intestinal = "#670092"),
    carcinoma = c(no = "#666666", yes = "#ff0000"),
    auc_quant = c("1"="#0000ff", "2"="#00ff00", "3"="#ff0000"),
    prognostic=c(no="#aaaaaa", fav="#00FFFF", unfav="#FFFF00"),
    de_cptac_rna = c(no="#aaaaaa", dn="#00FFFF", up="#FFFF00", na="#FFFFFF"), 
    de_cptac_prot = c(no="#aaaaaa", dn="#00FFFF", up="#FFFF00", na="#FFFFFF")
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
  x <- x  %>%
    filter(CodeClass != "Negative") %>%
    rename(entrez_id = GeneID) %>%
    group_by(entrez_id) %>%
    summarise(
      across(all_of(sample_cols), geomean)
    )
  return(x)
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

quantile_normalize <- function(x, cpm=TRUE, log=TRUE) {
  xn <- x
  if (cpm) {
    xn <- sweep(xn * 1e6, 2, colSums(xn), "/")
    xn <- apply(xn, 2, pmax, 1)
  }
  if (log) {
    xn <- log2(xn)
  }
  xn <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(xn)))
  rownames(xn) <- rownames(x)
  colnames(xn) <- colnames(x)
  if (log) {
    xn <- 2**(xn)
  }
  xn
}


# read input data
ipmn <- process_input(input_xlsx_file, 
                      negprobe_lod_quantile = qc_negprobe_lod_quantile)

# restrict all downstream analysis to specified segments
metadata <- filter(ipmn$metadata, segment %in% analysis_segments)
counts_probe <- select(ipmn$counts,
                       -all_of(ipmn$metadata$aoi_id),
                       all_of(metadata$aoi_id))

# merge multiple probes per gene to a single mean per gene
counts_gene <- merge_probes_to_genes(counts_probe, metadata$aoi_id)
gene_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                      keys=as.character(counts_gene$entrez_id),
                                      column="SYMBOL", 
                                      keytype="ENTREZID", 
                                      multiVals="first")
counts_gene <- mutate(counts_gene, gene_symbol = gene_symbols, .after = "entrez_id")

# filter aois based on qc parameters
metadata <- metadata %>%
  mutate(qcfail = (num_counts <= qc_min_counts | auc <= qc_min_auc | 
                     negprobe_geomean <= qc_min_negprobe_geomean))

# filter aois
fmetadata <- filter(metadata, !qcfail)
counts_gene <- select(counts_gene,
                      -all_of(metadata$aoi_id), 
                      all_of(fmetadata$aoi_id))

# compute quantiles of important qc metrics
fmetadata <- fmetadata %>%
  mutate(negprobe_quant = factor(ntile(negprobe_geomean, 4)),
         auc_quant = factor(ntile(auc, 4)),
         count_quant = factor(ntile(num_counts, 4)))

# fraction of aois expressing each genes
frac_expr <- calc_frac_expr(select(counts_gene, fmetadata$aoi_id), fmetadata$negprobe_lod)
counts_gene <- mutate(counts_gene, frac_expr = frac_expr, .after = "gene_symbol")

# filter genes with low expression
fcounts_gene <- filter(counts_gene, frac_expr > qc_min_frac_expr)

# background subtraction
gene_meta <- select(fcounts_gene, -fmetadata$aoi_id)
fcounts_gene_mat <- fcounts_gene %>%
  select(gene_symbol, all_of(fmetadata$aoi_id)) %>%
  column_to_rownames("gene_symbol")
fcounts_bgsub_gene_mat <- background_subtract(fcounts_gene_mat, fmetadata$negprobe_geomean)

# quantile normalization
fcounts_qnorm_gene_mat <- round(quantile_normalize(fcounts_bgsub_gene_mat))

# round to nearest count
fcounts_gene_mat <- round(fcounts_gene_mat)
fcounts_bgsub_gene_mat <- round(fcounts_bgsub_gene_mat)

# write processed data
write_xlsx(
  list(unfiltered_metadata = metadata,
       unfiltered_counts_probe = counts_probe,
       metadata = fmetadata,
       gene_meta = gene_meta,
       unfiltered_gene_counts = counts_gene,
       gene_count_matrix = as_tibble(fcounts_gene_mat, rownames="gene_symbol"),
       bgsub_gene_count_matrix = as_tibble(fcounts_bgsub_gene_mat, rownames="gene_symbol"),
       qnorm_gene_count_matrix = as_tibble(fcounts_qnorm_gene_mat, rownames="gene_symbol")),
  file.path(working_dir, processed_xlsx_file)
)


# setup samples
samples <- fmetadata

# extra params
samples <- mutate(samples, 
                  region = ifelse(segment == "panck", "epithelia", "stroma"))
samples$region <- factor(samples$region, levels=c("stroma", "epithelia"))

samples$slide_id <- factor(samples$slide_id)
samples$aoi_id <- factor(samples$aoi_id)
samples$path <- factor(samples$path, levels=c("lgd", "hgd"))
samples$histology <- factor(samples$histology)
samples$hgd_histology <- factor(samples$hgd_histology)
samples$carcinoma <- factor(samples$carcinoma)
samples$histpath <- factor(samples$histpath)
samples$roi_group <- factor(samples$roi_group)

table(samples$segment, samples$histology)
table(samples$histpath)

# color scales for plots
color_scales <- get_color_scales(unique(samples$slide_id))

# panck data
panck = list()
panck$samples <- filter(samples, segment == "panck")
panck$samples$region <- factor(panck$samples$region)
panck$samples$histology <- factor(panck$samples$histology)
panck$samples$carcinoma <- factor(panck$samples$carcinoma)
panck$samples$histpath <- factor(panck$samples$histpath)
panck$bgsub_gene_mat <- fcounts_bgsub_gene_mat[, panck$samples$aoi_id]
panck$qnorm_gene_mat <- fcounts_qnorm_gene_mat[, panck$samples$aoi_id]
panck$gene_meta <- gene_meta

# cd45 data
cd45 = list()
cd45$samples <- filter(samples, segment == "cd45")
cd45$samples$region <- factor(cd45$samples$region)
cd45$samples$histology <- factor(cd45$samples$histology)
cd45$samples$path <- factor(cd45$samples$path)
cd45$samples$carcinoma <- factor(cd45$samples$carcinoma)
cd45$samples$histpath <- factor(cd45$samples$histpath)
cd45$bgsub_gene_mat <- fcounts_bgsub_gene_mat[, cd45$samples$aoi_id]
cd45$qnorm_gene_mat <- fcounts_qnorm_gene_mat[, cd45$samples$aoi_id]
cd45$gene_meta <- gene_meta

# sma data
sma = list()
sma$samples <- filter(samples, segment == "sma")
sma$samples$region <- factor(sma$samples$region)
sma$samples$histology <- factor(sma$samples$histology)
sma$samples$path <- factor(sma$samples$path)
sma$samples$carcinoma <- factor(sma$samples$carcinoma)
sma$samples$histpath <- factor(sma$samples$histpath)
sma$bgsub_gene_mat <- fcounts_bgsub_gene_mat[, sma$samples$aoi_id]
sma$qnorm_gene_mat <- fcounts_qnorm_gene_mat[, sma$samples$aoi_id]
sma$gene_meta <- gene_meta
table(sma$samples$path)
table(sma$samples$histology)
table(sma$samples$histpath)

##################################################
#
# Quality Control Plots
#
##################################################

table(metadata$segment)
summary(metadata$AOINucleiCount)
summary(metadata$num_counts)
summary(metadata$pct_expressed)
table(metadata$qcfail, metadata$segment)
table(counts_gene$frac_expr > qc_min_frac_expr)
summary(filter(metadata, segment == "panck") %>% select(AOINucleiCount))
summary(filter(metadata, segment != "panck") %>% select(AOINucleiCount))

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
  annotate("text", x=10, y=2e3, label = corlab) +
  theme_minimal()
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
  annotate("text", x=1e3, y=2e3, label = corlab) +
  theme_minimal()
f <- file.path(plot_dir, "scatter_area_vs_counts.pdf")
ggsave(f, plot=p, device="pdf", width = 4, height = 3)

# figure: total counts vs snAUC
p <- ggplot(metadata, aes(x=num_counts, y=auc)) +
  geom_point() + 
  geom_hline(yintercept = qc_min_auc, linetype = "dashed", color = "red") +
  geom_vline(xintercept = qc_min_counts, linetype = "dashed", color = "red") +
  scale_x_continuous(trans="log10") +
  xlab("Total Counts") +
  ylab("snAUC") +
  labs(color = "Segment") +
  theme_minimal() +
  facet_wrap(~factor(segment, levels=c("cd45", "sma", "panck")))
p
f <- file.path(plot_dir, "scatter_counts_vs_auc_bysegment.pdf")
ggsave(f, plot=p, device="pdf", width = 8, height = 4)

# figure: boxplot path versus nuclei count by slide
p <- ggplot(filter(fmetadata, segment == "panck"), aes(x=path, y=AOINucleiCount)) +
  geom_boxplot(aes(fill=path), outlier.shape=NA) +
  geom_point(aes(fill=path), color="black", size=2, position=position_jitterdodge(jitter.width=0.15)) +
  scale_fill_manual(values=color_scales$path) +
  scale_color_manual(values=color_scales$path) +
  xlab("Pathology") +
  ylab("Nuclei Count") +
  labs(fill = "Pathology") +
  theme(legend.position="bottom") +
  facet_grid(cols = vars(slide_id)) 
f <- file.path(plot_dir, "boxplot_path_vs_nuclei_byslide.pdf")
ggsave(f, plot=p, device="pdf", width = 8, height = 4)

# figure: boxplot path versus total counts by slide
p <- ggplot(filter(fmetadata, segment == "panck"), aes(x=path, y=num_counts)) +
  geom_boxplot(aes(fill=path), outlier.shape=NA) +
  geom_point(aes(fill=path), color="black", size=2, position=position_jitterdodge(jitter.width=0.15)) +
  scale_y_continuous(trans="log10") +
  scale_fill_manual(values=color_scales$path) +
  scale_color_manual(values=color_scales$path) +
  xlab("Pathology") +
  ylab("Total Counts") +
  labs(fill = "Pathology") +
  theme(legend.position="bottom") +
  facet_grid(cols = vars(slide_id)) 
f <- file.path(plot_dir, "boxplot_path_vs_counts_byslide.pdf")
ggsave(f, plot=p, device="pdf", width = 8, height = 4)

# figure: boxplot path versus auc by slide
p <- ggplot(filter(fmetadata, segment == "panck"), aes(x=path, y=auc)) +
  geom_boxplot(aes(fill=path), outlier.shape=NA) +
  geom_point(aes(fill=path), color="black", size=2, position=position_jitterdodge(jitter.width=0.15)) +
  scale_y_continuous(trans="log10") +
  scale_fill_manual(values=color_scales$path) +
  scale_color_manual(values=color_scales$path) +
  xlab("Pathology") +
  ylab("snAUC") +
  labs(fill = "Pathology") +
  theme(legend.position="bottom") +
  facet_grid(cols = vars(slide_id)) 
f <- file.path(plot_dir, "boxplot_path_vs_auc_byslide.pdf")
ggsave(f, plot=p, device="pdf", width = 8, height = 4)


# figure: gene filtering fraction expressed cutoff
p <- ggplot(counts_gene, aes(x=frac_expr)) +
  geom_histogram(binwidth=0.05, color="black", fill="grey") +
  geom_vline(xintercept = qc_min_frac_expr, linetype="dashed", color="red") +
  theme_minimal() +
  xlab("Fraction Expressed Above Background") +
  ylab("# of Genes")
f <- file.path(plot_dir, "histogram_frac_expr.pdf")
ggsave(f, plot=p, device="pdf", width = 4, height = 4)


# figure: ridge plot unnormalized counts (expressed no/yes)
x <- counts_gene %>%
  select(gene_symbol, all_of(fmetadata$aoi_id)) %>%
  column_to_rownames("gene_symbol")
x$expressed <- ifelse(counts_gene$frac_expr <= qc_min_frac_expr, "no", "yes")
y <- select(fmetadata, aoi_id, slide_id, segment, hgd_histology, path, count_quant, auc_quant)
x <- x %>%
  pivot_longer(fmetadata$aoi_id, names_to="aoi_id", values_to="count") %>%
  inner_join(y, by=c("aoi_id"="aoi_id"))
p <- x %>%
  ggplot(aes(x=count, y=factor(aoi_id), fill=factor(expressed))) +
  stat_density_ridges(alpha=0.7, scale=10, rel_min_height=0.001, quantile_lines=TRUE) + 
  scale_x_continuous(trans="log10") +
  xlab("Raw Counts") +
  ylab("AOIs") +
  labs(fill = "Expressed") +
  theme_ridges() + 
  theme(axis.text = element_text(size=5))
f <- file.path(plot_dir, "ridge_plot_rawcounts_vs_qcfilter.pdf")
ggsave(f, plot=p, device="pdf", width = 8, height = 8)


# gene count ridge plots
x <- counts_gene %>%
  select(gene_symbol, all_of(fmetadata$aoi_id)) %>%
  column_to_rownames("gene_symbol")
x <- round(background_subtract(x, fmetadata$negprobe_geomean))
x <- DGEList(counts = x)
x <- calcNormFactors(x, method="upperquartile")
x <- cpm(x, log=TRUE, prior.count=0)
x <- as_tibble(x, rownames="gene_symbol")
x$expressed <- ifelse(counts_gene$frac_expr <= qc_min_frac_expr, "no", "yes")
y <- select(fmetadata, aoi_id, slide_id, segment, hgd_histology, path, count_quant, auc_quant)
x <- x %>%
  pivot_longer(fmetadata$aoi_id, names_to="aoi_id", values_to="count") %>%
  inner_join(y, by=c("aoi_id"="aoi_id"))

p <- x %>%
  ggplot(aes(x=count, y=factor(aoi_id), fill=factor(expressed))) +
  stat_density_ridges(alpha=0.7, scale=10, rel_min_height=0.001, quantile_lines=TRUE) + 
  xlab("Counts per million reads (CPM)") +
  ylab("AOIs") +
  labs(fill = "Expressed") +
  theme_ridges() + 
  theme(axis.text = element_text(size=5))
f <- file.path(plot_dir, "ridge_plot_aoi_vs_qcfilter.pdf")
ggsave(f, plot=p, device="pdf", width = 8, height = 8)

# filtered gene ridge plots
x <- DGEList(counts = fcounts_bgsub_gene_mat)
x <- calcNormFactors(x, method = "upperquartile")
x <- cpm(x, log=TRUE, prior.count=0)
x <- as_tibble(x, rownames="gene_symbol")
y <- select(fmetadata, aoi_id, slide_id, segment, histology, path, count_quant, auc_quant)
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
f <- file.path(plot_dir, "ridge_plot_aoi_vs_slide.pdf")
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
f <- file.path(plot_dir, "ridge_plot_path.pdf")
ggsave(f, plot=p, device="pdf", width = 4, height = 4)


p <- filter(x, segment == "panck") %>%
  ggplot(aes(x=count, y=histology, fill=histology)) +
  stat_density_ridges(alpha=0.7, scale=3, rel_min_height=0.001, quantile_lines=TRUE) + 
  scale_fill_manual(values = color_scales$histology2) +
  xlab("log2(cpm)") +
  ylab("Histology") +
  labs(fill = "Histology") +
  theme_ridges() + 
  theme(legend.position="bottom")
f <- file.path(plot_dir, "ridge_plot_hist.pdf")
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

normalize_voom_q3 <- function(m, s, g) {
  y <- DGEList(counts = m, samples = s, genes = g)
  y <- calcNormFactors(y, method = "upperquartile")
  v <- voom(y)
  return(v$E)
}

#
# full dataset (epithelia plus stroma)
#
voom_bgsub_q3 <- normalize_voom_q3(fcounts_bgsub_gene_mat, samples, gene_meta)
pca_res <- run_pca(voom_bgsub_q3, samples)
pca_tbl <- pca_res$tbl
tsne_res <- run_tsne(voom_bgsub_q3, samples, perp=11)
tsne_tbl <- tsne_res$tbl
# write_xlsx(tsne_tbl, file.path(working_dir, tsne_results_xlsx_file))
# tsne_tbl <- read_xlsx(file.path(working_dir, tsne_results_xlsx_file), sheet = "Sheet1")

# pca
pr <- pca_res$pr
pct <- round(pr$sdev / sum(pr$sdev) * 100, 2)
pct <- paste(colnames(pr$x), " (", paste(as.character(pct), "%", ")", sep=""), sep="")
p <- ggplot(pca_tbl, aes(x=PC1, y=PC2, color=region)) +
  geom_point(size=3) +
  scale_color_manual(values=color_scales$region) +
  theme_minimal() +
  labs(x=pct[1], y=pct[2])
p
p <- ggplot(pca_tbl, aes(x=PC1, y=PC2, color=histology)) +
  geom_point(size=3) +
  scale_color_manual(values=color_scales$histology) +
  theme_minimal() +
  labs(x=pct[1], y=pct[2])
p


# tsne
p <- ggplot(tsne_tbl, aes(x=tsne1, y=tsne2, color=region)) +
  geom_point(size=2) +
  scale_color_manual(values=color_scales$region) +
  ggtitle('tSNE clustering') +
  theme_minimal()
p
f <- file.path(plot_dir, "tsne_allgenes_region.pdf")
ggsave(f, plot=p, device="pdf", width=5, height=4)

p <- ggplot(tsne_tbl, aes(x=tsne1, y=tsne2, color=histology)) +
  geom_point(size=2) +
  scale_color_manual(values=color_scales$histology) +
  theme_minimal()
p
f <- file.path(plot_dir, "tsne_allgenes_histology.pdf")
ggsave(f, plot=p, device="pdf", width=5, height=4)


#
# panck data (epithelial)
#
panck_tsne_tbl <- filter(tsne_tbl, segment == "panck")

p <- ggplot(panck_tsne_tbl, aes(x=tsne1, y=tsne2, color=histology, shape=path)) +
  geom_point(size=3) +
  scale_color_manual(values=color_scales$histology) +
  theme_minimal()
p
f <- file.path(plot_dir, "tsne_panck_histology.pdf")
ggsave(f, plot=p, device="pdf", width=6, height=5)

# slide
p <- ggplot(panck_tsne_tbl, aes(x=tsne1, y=tsne2, color=slide_id, shape=path)) +
  geom_point(size=3) +
  scale_color_manual(values=color_scales$slide_id) +
  theme_minimal()
p
f <- file.path(plot_dir, "tsne_panck_slide_id.pdf")
ggsave(f, plot=p, device="pdf", width=6, height=5)

# carcinoma
p <- ggplot(panck_tsne_tbl, aes(x=tsne1, y=tsne2, color=carcinoma, shape=path)) +
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

limmaVoom <- function(y, design, contrasts) {
  # limma voom de analysis
  v <- voom(y, design, plot=TRUE)
  fit <- lmFit(v, design)
  fit <- contrasts.fit(fit, contrasts)
  fit <- eBayes(fit)
  return(list(voom=v$E, fit=fit))
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
design <- model.matrix(~0 + histology, data=panck$samples)
colnames(design) <- levels(panck$samples$histology)
colnames(design)
contrasts = makeContrasts(hist_pb_vs_gf = pb - gf,
                          hist_pb_vs_int = pb - intestinal,
                          hist_int_vs_gf = intestinal - gf,
                          hist_int = intestinal - (pb + gf)/2,
                          levels=design)

# voom with scale normalization
method <- "limma_voom_bgsub_q3"
y <- DGEList(counts=panck$bgsub_gene_mat,
             samples=panck$samples,
             group=panck$samples$histology,
             genes=gene_meta)
y <- calcNormFactors(y, method="upperquartile")
x <- limmaVoom(y, design, contrasts)
res <- processLimma(x$fit, contrasts, method, de_padj_cutoff, de_log2fc_cutoff, de_result_dir)
hist_de <- res
hist_voom_lcpm <- as_tibble(x$voom, rownames="gene_symbol") %>%
  column_to_rownames("gene_symbol")
hist_cpm <- cpm(y)

# merged de analyses
hist_de_merged <- select(hist_de, gene, analysis, avgexpr, log2fc, padj, de) %>%
  pivot_wider(names_from = analysis, values_from = c(log2fc, padj, de))

hist_de_merged <- hist_de_merged %>%
  mutate(de_pb = case_when(de_hist_pb_vs_gf == "up" & de_hist_pb_vs_int == "up" ~ "up",
                           de_hist_pb_vs_gf == "dn" & de_hist_pb_vs_int == "dn" ~ "dn",
                           TRUE ~ "no"),
         de_int = case_when(de_hist_int_vs_gf == "up" & de_hist_pb_vs_int == "dn" ~ "up",
                            de_hist_int_vs_gf == "dn" & de_hist_pb_vs_int == "up" ~ "dn",
                            TRUE ~ "no"),
         de_gf = case_when(de_hist_pb_vs_gf == "dn" & de_hist_int_vs_gf == "dn" ~ "up",
                           de_hist_pb_vs_gf == "up" & de_hist_int_vs_gf == "up" ~ "dn",
                           TRUE ~ "no"))

# write results
write_xlsx(list(samples = panck$samples,
                gene_meta = gene_meta,
                de = hist_de,
                merged_de = hist_de_merged,
                voom_lcpm = hist_voom_lcpm,
                cpm = as_tibble(hist_cpm, rownames="gene_symbol")),
           file.path(working_dir, de_results_hist_xlsx))


#
# grade
#
table(panck$samples$histpath)
design <- model.matrix(~0 + histpath, data=panck$samples)
colnames(design)
colnames(design) <- levels(panck$samples$histpath)
colnames(design)
contrasts = makeContrasts(path_pbhgd_vs_intlgd = panck_pbhgd - panck_intlgd,
                          path_pbhgd_vs_pblgd = panck_pbhgd - panck_pblgd,
                          path_inthgd_vs_intlgd = panck_inthgd - panck_intlgd,
                          path_inthgd_vs_pblgd = panck_inthgd - panck_pblgd,
                          path_hgd_vs_lgd = (panck_pbhgd + panck_inthgd)/2 - (panck_intlgd + panck_pblgd)/2,
                          path_pbhgd_vs_inthgd = panck_pbhgd - panck_inthgd,
                          path_pblgd_vs_intlgd = panck_pblgd - panck_intlgd,
                          levels=design)

# voom with scale normalization
method <- "limma_voom_bgsub_q3"
y <- DGEList(counts=panck$bgsub_gene_mat,
             samples=panck$samples,
             group=panck$samples$histpath,
             genes=gene_meta)
y <- calcNormFactors(y, method="upperquartile")
x <- limmaVoom(y, design, contrasts)
res <- processLimma(x$fit, contrasts, method, de_padj_cutoff, de_log2fc_cutoff, de_result_dir)
grade_de <- res
grade_voom_lcpm <- as_tibble(x$voom, rownames="gene_symbol") %>%
  column_to_rownames("gene_symbol")
grade_cpm <- cpm(y)

# merged de analyses
grade_de_merged <- select(grade_de, gene, analysis, avgexpr, log2fc, padj, de) %>%
  pivot_wider(names_from = analysis, values_from = c(log2fc, padj, de))

grade_de_merged <- grade_de_merged %>%
  mutate(subtype = case_when(de_path_pbhgd_vs_pblgd == "up" & de_path_inthgd_vs_intlgd == "up" ~ "hgd",
                             de_path_pbhgd_vs_pblgd == "dn" & de_path_inthgd_vs_intlgd == "dn" ~ "lgd",
                             de_path_pbhgd_vs_pblgd == "up" ~ "pbhgd",
                             de_path_inthgd_vs_intlgd == "up" ~ "inthgd",
                             de_path_pbhgd_vs_pblgd == "dn" ~ "pblgd",
                             de_path_inthgd_vs_intlgd == "dn" ~ "intlgd",
                             TRUE ~ "none"))

# write gene expression data
write_xlsx(list(samples = samples,
                gene_meta = gene_meta,
                de = grade_de,
                merged_de = grade_de_merged,
                voom_lcpm = grade_voom_lcpm,
                cpm = as_tibble(grade_cpm, rownames="gene_symbol")),
           file.path(working_dir, de_results_grade_xlsx))


#
# cd45 data
#
method <- "limma_voom_bgsub_q3"
design <- model.matrix(~0 + histpath, data=cd45$samples)
colnames(design)
colnames(design) <- levels(cd45$samples$histpath)
colnames(design)
contrasts = makeContrasts(cd45_hgd_vs_lgd = (cd45_inthgd + cd45_pbhgd)/2 - (cd45_intlgd + cd45_pblgd)/2,
                          cd45_pbhgd_vs_pblgd = cd45_pbhgd - cd45_pblgd,
                          cd45_inthgd_vs_intlgd = cd45_inthgd - cd45_intlgd,
                          levels=design)

# voom with scale normalization
y <- DGEList(counts=cd45$bgsub_gene_mat,
             samples=cd45$samples,
             genes=gene_meta)
y <- calcNormFactors(y, method="upperquartile")
x <- limmaVoom(y, design, contrasts)
res <- processLimma(x$fit, contrasts, method, de_padj_cutoff, de_log2fc_cutoff, de_result_dir)
cd45_de <- res
cd45_voom_lcpm <- as_tibble(x$voom, rownames="gene_symbol") %>%
  column_to_rownames("gene_symbol")

#
# sma data
#
method <- "limma_voom_bgsub_q3"
design <- model.matrix(~0 + histpath, data=sma$samples)
colnames(design)
colnames(design) <- levels(sma$samples$histpath)
colnames(design)
contrasts = makeContrasts(sma_hgd_vs_lgd = (sma_inthgd + sma_pbhgd)/2 - (sma_intlgd + sma_pblgd)/2,
                          sma_pbhgd_vs_pblgd = sma_pbhgd - sma_pblgd,
                          sma_inthgd_vs_intlgd = sma_inthgd - sma_intlgd,
                          levels=design)

# voom with scale normalization
y <- DGEList(counts=sma$bgsub_gene_mat,
             samples=sma$samples,
             genes=gene_meta)
y <- calcNormFactors(y, method="upperquartile")
x <- limmaVoom(y, design, contrasts)
res <- processLimma(x$fit, contrasts, method, de_padj_cutoff, de_log2fc_cutoff, de_result_dir)
sma_de <- res
sma_voom_lcpm <- as_tibble(x$voom, rownames="gene_symbol") %>%
  column_to_rownames("gene_symbol")

##################################################
#
#
# Heatmaps
#
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
  # TODO: add duct / epithelial filtered data
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
    select(gene_symbol, prognostic, de_cptac_rna, de_cptac_prot) %>%
    column_to_rownames("gene_symbol") %>%
    as.data.frame()
}

# define background gene set
cta_genes <- unique(hist_de$gene)
# write cta genes to file (background gene set)
write.table(cta_genes, file=file.path(working_dir, "cta_gene_symbols.txt"),
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

gene_annot <- gene_meta %>%
  left_join(hpa_pdac, by=c("gene_symbol"="gene_symbol")) %>%
  left_join(cptac_de_merged, by=c("gene_symbol" = "gene_symbol")) %>%
  mutate(de_cptac_rna = replace_na(de_cptac_rna, "na"),
         de_cptac_prot = replace_na(de_cptac_prot, "na"))

# 
# histology analysis heatmap
#
# union of subtype specific genes
hist_pb_genes <- filter(hist_de_merged, de_pb != "no") %>% select(gene) %>% pull()
hist_int_genes <- filter(hist_de_merged, de_int != "no") %>% select(gene) %>% pull()
hist_gf_genes <- filter(hist_de_merged, de_gf != "no") %>% select(gene) %>% pull()
hist_genes <- union(hist_pb_genes, union(hist_gf_genes, hist_int_genes))
x <- hist_voom_lcpm[hist_genes,]

annot_row <- get_row_annot(hist_genes, gene_annot)
annot_col <- column_to_rownames(panck$samples, "aoi_id") %>% 
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
# union of subtype hgd genes
grade_genes <- filter(grade_de_merged, de_path_hgd_vs_lgd != "no") %>% select(gene) %>% pull()
length(grade_genes)

x <- hist_voom_lcpm[grade_genes,]
annot_row <- get_row_annot(grade_genes, gene_annot)
annot_col <- column_to_rownames(panck$samples, "aoi_id") %>% 
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
# Volcano plots
#
#
##################################################

# TODO

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
label_data <- bind_rows(filter(label_data, subtype == "hgd") %>% slice_max(log2fc_path_pbhgd_vs_pblgd, n = 10),
                        filter(label_data, subtype == "lgd") %>% slice_min(log2fc_path_pbhgd_vs_pblgd, n = 10),
                        filter(label_data, subtype == "hgd") %>% slice_max(log2fc_path_inthgd_vs_intlgd, n = 10),
                        filter(label_data, subtype == "lgd") %>% slice_min(log2fc_path_inthgd_vs_intlgd, n = 10),
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
colnames(hist_de_merged)

x <- hist_de_merged %>%
  mutate(subtype = case_when(de_hist_pb_vs_gf == "up" & de_hist_pb_vs_int == "up" ~ "pb",
                             de_hist_pb_vs_gf == "dn" & de_hist_int_vs_gf == "dn" ~ "gf",
                             de_hist_int_vs_gf == "up" & de_hist_pb_vs_int == "dn" ~ "int",
                             de_hist_int_vs_gf == "dn" & de_hist_pb_vs_int == "up" ~ "pbgf",
                             TRUE ~ "none"))

color_scale = c(none = "#888888",
                pb = "#32cd32",
                int = "#670092",
                gf = "#1677c4",
                pbgf = "#ff6700")

size_scale = c(none = 1,
               pb = 2,
               int = 2,
               gf = 2,
               pbgf = 2)

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
                        filter(x, subtype == "pbgf") %>% slice_min(log2fc_hist_int_vs_gf, n = 10))
# manually add MUC1 as this is featured in plot but not in top 10
label_data <- bind_rows(label_data, filter(x, gene == "MUC1"))
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
#  p <- ggplot(x, aes(x=reorder(histpath, gene), y=gene)) + 
  p <- ggplot(x, aes(x=histpath, y=gene)) + 
    geom_boxplot(aes(fill = histpath), outlier.shape=NA) +
    geom_jitter(width=0.1) +
    scale_fill_manual(values=color_scales$histpath) +
    xlab("Histopathology") +
    ylab("log2 Normalized CPM") + 
    labs(fill = "Histopath") +
    ggtitle(paste0("Gene: ", gene_symbol)) +
    theme_minimal() + 
    theme(axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))
  return(p)
}


# plot all de genes in heatmaps
for (g in union(hist_genes, grade_genes)) {
  p <- plot_hist_boxplot(panck$samples, hist_voom_lcpm, g)
  f <- file.path(gene_plot_dir, paste0("boxplot_", g, "_hist.pdf"))
  ggsave(f, p, width=5, height=5)
  p <- plot_grade_boxplot(panck$samples, grade_voom_lcpm, g)
  f <- file.path(gene_plot_dir, paste0("boxplot_", g, "_grade.pdf"))
  ggsave(f, p, width=5, height=5)
}

# manually plot individual genes
p <- plot_hist_boxplot(panck$samples, hist_voom_lcpm, "MUC1")
ggsave(file.path(plot_dir, "hist_boxplot_muc1.pdf"), p, width=3, height=5)
p <- plot_hist_boxplot(panck$samples, hist_voom_lcpm, "MUC4")
ggsave(file.path(plot_dir, "hist_boxplot_muc4.pdf"), p, width=3, height=5)
p <- plot_hist_boxplot(panck$samples, hist_voom_lcpm, "CLU")
ggsave(file.path(plot_dir, "hist_boxplot_clu.pdf"), p, width=3, height=5)
p <- plot_hist_boxplot(panck$samples, hist_voom_lcpm, "RBP4")
ggsave(file.path(plot_dir, "hist_boxplot_rbp4.pdf"), p, width=3, height=5)


# lymphocyte marker
plot_grade_boxplot(samples, voom_bgsub_q3, "PTPRC")
plot_grade_boxplot(samples, voom_bgsub_q3, "CD3D")
plot_grade_boxplot(samples, voom_bgsub_q3, "IL7R")

# EMT-like cells
plot_grade_boxplot(samples, voom_bgsub_q3, "ZEB1")

# ductal cells
plot_grade_boxplot(samples, voom_bgsub_q3, "TACSTD2")
plot_grade_boxplot(samples, voom_bgsub_q3, "KRT7")
# macrophages
plot_grade_boxplot(samples, voom_bgsub_q3, "C1QA")
plot_grade_boxplot(samples, voom_bgsub_q3, "C1QB")

# myeloid-derived suppressor cells
plot_grade_boxplot(samples, voom_bgsub_q3, "S100A9")
plot_grade_boxplot(samples, voom_bgsub_q3, "CCL3")

# KDR and VWF for endothelial cells; and FCER1A and CD1 for dendritic cells
# DC2 dendritic cells
plot_grade_boxplot(samples, voom_bgsub_q3, "THBD")

# fibroblast marker
plot_grade_boxplot(samples, voom_bgsub_q3, "ACTA2")
plot_grade_boxplot(samples, voom_bgsub_q3, "COL1A1")

# inflammatory CAF
plot_grade_boxplot(samples, voom_bgsub_q3, "COL3A1")
plot_grade_boxplot(samples, voom_bgsub_q3, "CXCL12")
plot_grade_boxplot(samples, voom_bgsub_q3, "TTR")

# proliferation markers
plot_grade_boxplot(samples, voom_bgsub_q3, "MKI67")
plot_grade_boxplot(samples, voom_bgsub_q3, "TOP2A")
plot_grade_boxplot(samples, voom_bgsub_q3, "CENPF")


grade_de %>% filter(gene == "CENPF")

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

#
# Setup Gene Sets
# 
# cptac gene sets
cptac_gs <- get_cptac_gene_sets(cptac_de_data)

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

#
# PDAC GSEA Analysis
#
prefix <- "gsea_pdac"
method <- "limma_voom_bgsub_q3"
analyses <- c("path_inthgd_vs_intlgd", "path_pbhgd_vs_pblgd", "path_hgd_vs_lgd")
de <- bind_rows(hist_de, grade_de)
gs <- gruetzmann_gs
gs <- append(gs, hpa_panc_gs)
gs <- append(gs, cptac_gs)
gs_pdac <- gs

gsea = NULL
for (a in analyses) {
  ranks <- get_de_ranks(de, a)
  res <- fgsea(pathways = gs, stats = ranks, minSize = 10, eps = 0, nPermSimple = 10000)
  res <- mutate(res, method = method, analysis = a)
  gsea <- bind_rows(gsea, res)
}
gsea_pdac <- as_tibble(gsea)

# write gsea results
f <- file.path(working_dir, paste0(prefix, "_results.tsv"))
data.table::fwrite(gsea, file=f, sep="\t", sep2=c("", " ", ""))
f <- file.path(working_dir, paste0(prefix, "_results.xlsx"))
write_xlsx(gsea, f)

# construct bar plots comparing analyses
x <- filter(gsea_pdac, analysis %in% analyses)
x$analysis <- factor(x$analysis, levels=analyses)
x <- mutate(x, sig = ifelse(padj < 0.05, ifelse(NES < 0, "dn", "up"), "no"))
p <- ggplot(x, aes(x=analysis, y=NES)) +
  geom_col(aes(fill = sig)) +
  scale_fill_manual(values = c("no"="#aaaaaa", "dn"="#0000cc", "up"="#cc0000")) +
  coord_flip() +
  labs(x="Analysis", y="Normalized Enrichment Score",
       title="PDAC Gene Signatures") + 
  theme_minimal() +
  facet_wrap(~ pathway, ncol = 1)
p
ggsave(file.path(plot_dir, "gsea_pdac_barplots.pdf"), p, width=6, height=8)


#
# discovery analysis across msigdb gene sets
#
de <- bind_rows(hist_de, grade_de)
method <- "limma_voom_bgsub_q3"
#analyses <- c("path_hgd_vs_lgd")
#analyses <- c("path_hgd_vs_lgd", "path_pbhgd_vs_pblgd", "path_inthgd_vs_intlgd")
analyses <- unique(grade_de$analysis)
color_scale <- c("path_hgd_vs_lgd" = "#ff0000",
                 "path_inthgd_vs_intlgd" = "#670092",
                 "path_pbhgd_vs_pblgd" = "#32cd32")

gs_hallmark <- filter(msigdb_gene_sets, gs_cat == "H")
gs_hallmark <- split(x = gs_hallmark$gene_symbol, f = gs_hallmark$gs_name)
gs_gobp <- filter(msigdb_gene_sets, gs_cat == "C5", gs_subcat == "GO:BP")
gs_gobp <- split(x = gs_gobp$gene_symbol, f = gs_gobp$gs_name)

#
# Hallmark gene set analysis
#
prefix <- "gsea_hallmark"
gs <- gs_hallmark
gsea = NULL
for (a in analyses) {
  ranks <- get_de_ranks(de, a)
  res <- fgsea(pathways = gs, stats = ranks, minSize = 10, eps = 0, nPermSimple = 10000)
  res <- mutate(res, method = method, analysis = a)
  gsea <- bind_rows(gsea, res)
  
  for (i in 1:nrow(res)) {
    x <- res[i,]
    if (x$padj >= gsea_padj_cutoff) { next }
    print(x$pathway)
    p <- plot_gsea_enrichment(x, ranks, gs)
    f <- file.path(gsea_plot_dir, paste0(prefix, "_", a, "_gs_", x$pathway, ".pdf"))
    ggsave(f, plot=p, device="pdf", width=5, height=3)
  }
}
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
gs <- gs_gobp
gsea = NULL
for (a in analyses) {
  ranks <- get_de_ranks(de, a)
  res <- fgsea(pathways = gs, stats = ranks, minSize = 10, eps = 0, nPermSimple = 10000)
  res <- mutate(res, method = method, analysis = a)
  gsea <- bind_rows(gsea, res)

  for (i in 1:nrow(res)) {
    x <- res[i,]
    if (x$padj >= gsea_padj_cutoff) { next }
    print(x$pathway)
    p <- plot_gsea_enrichment(x, ranks, gs)
    f <- file.path(gsea_plot_dir, paste0(prefix, "_", a, "_gs_", x$pathway, ".pdf"))
    ggsave(f, plot=p, device="pdf", width=5, height=3)
  }
}
gsea_gobp <- as_tibble(gsea) %>%
  left_join(select(msigdb_gene_sets, gs_cat, gs_subcat, gs_name) %>% distinct(), 
            by=c("pathway"="gs_name"))

# write gsea results
f <- file.path(working_dir, paste0(prefix, "_results.tsv"))
data.table::fwrite(gsea, file=f, sep="\t", sep2=c("", " ", ""))
f <- file.path(working_dir, paste0(prefix, "_results.xlsx"))
write_xlsx(gsea, f)


#
# GSEA Volcano Plots
#

# hallmark
x <- gsea_hallmark %>% filter(analysis == "path_hgd_vs_lgd")
# x <- gsea_hallmark %>% filter(analysis %in% c("path_hgd_vs_lgd", "path_pbhgd_vs_pblgd", "path_inthgd_vs_intlgd"))
x$abbrev <- gsub("HALLMARK_", "", x$pathway)
p <- ggplot(x, aes(x=NES, y=-log10(padj), color=analysis)) + 
  geom_point(data=subset(x, padj > 0.01), size=1, alpha=0.4) + 
  geom_point(data=subset(x, padj <= 0.01), size=2, alpha=0.8) +
  geom_text_repel(data=subset(x, padj <= 0.01), color="black", size=3, aes(label=abbrev), max.overlaps=Inf) +
  scale_color_manual(values = color_scale) +
  ylab("-log10(adjusted p-value)") +
  xlab("NES") +
  ggtitle("MSigDB Hallmark Gene Sets") +
  theme_minimal() +
  theme(legend.position="bottom") + 
  theme(axis.line = element_line(color = "black")) +
  coord_flip() + 
  facet_wrap(~ analysis, ncol = 1)
p
ggsave(file.path(plot_dir, "gsea_hallmark_volcano.pdf"), plot=p, width=5, height=3)


# GO
x <- gsea_gobp %>%
  filter(analysis == "path_hgd_vs_lgd")
#x <- gsea_gobp %>%
#  filter(analysis %in% c("path_hgd_vs_lgd", "path_pbhgd_vs_pblgd", "path_inthgd_vs_intlgd"))
x$abbrev <- gsub("GOBP_", "", x$pathway)
p <- ggplot(x, aes(x=NES, y=-log10(padj), color=analysis)) + 
  geom_point(data=subset(x, padj > 0.01), size=1, alpha=0.4) + 
  geom_point(data=subset(x, padj <= 0.01), size=2, alpha=0.8) +
  geom_text_repel(data=subset(x, padj <= 0.01), color="black", size=3, aes(label=abbrev), max.overlaps=Inf) +
  scale_color_manual(values = color_scale) +
  ylab("-log10(adjusted p-value)") +
  xlab("NES") +
  ggtitle("GOBP Gene Sets") +
  theme_minimal() +
  theme(legend.position="bottom") + 
  theme(axis.line = element_line(color = "black")) +
  coord_flip() +
  facet_wrap(~ analysis, ncol = 1)
p
ggsave(file.path(plot_dir, "gsea_gobp_volcano.pdf"), plot=p, width=5, height=3)

#
# GSEA Leading Edge
#
# merge de analyses
gsea_merged <- bind_rows(gsea_hallmark, gsea_gobp)
# gsea_merged <- gsea_hallmark

# perform leading edge analysis of top concepts
analyses <- c("path_hgd_vs_lgd")
# analyses <- c("path_hgd_vs_lgd", "path_pbhgd_vs_pblgd", "path_inthgd_vs_intlgd")
# analyses <- unique(gsea_merged$analysis)
sig_cutoff <- 0.01

x <- gsea_merged %>%
  filter(analysis %in% analyses)
sig_pathways <- x %>%
  filter(padj < sig_cutoff) %>% select(pathway) %>% distinct() %>% pull()
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
  # if (padj >= 0.05) { next }
  for (j in x$leadingEdge[i]) {
    m[j, p] <- sign(nes)
  }
}

# filter genes only found in a single pathway
m <- m[rowSums(m) >= 2, ]
# filter empty columns
m <- m[, colSums(m) > 1]
# order by number of pathways sharing the gene
m <- m[order(rowSums(m), decreasing=TRUE),]
m <- m[,order(colSums(m), decreasing=TRUE)]

# annotate rows
grade_subtypes <- grade_de_merged %>%
  mutate(grade_subtype = case_when(de_path_pbhgd_vs_pblgd == "up" & de_path_inthgd_vs_intlgd == "up" ~ "hgd",
                                   de_path_pbhgd_vs_pblgd == "dn" & de_path_inthgd_vs_intlgd == "dn" ~ "lgd",
                                   de_path_pbhgd_vs_pblgd == "up" ~ "pbhgd",
                                   de_path_inthgd_vs_intlgd == "up" ~ "inthgd",
                                   de_path_pbhgd_vs_pblgd == "dn" ~ "pblgd",
                                   de_path_inthgd_vs_intlgd == "dn" ~ "intlgd",
                                   TRUE ~ "none")) %>%
  select(gene, grade_subtype)

annot_colors <- color_scales
annot_colors$grade_subtype =  c(none = "#aaaaaa",
                                hgd = "#ff0000",
                                lgd = "#0000ff",
                                pbhgd = "#32cd32",
                                inthgd = "#670092",
                                pblgd = "#1677c4",
                                intlgd = "#ff6700")

row_annot <- gene_annot %>%
  left_join(grade_subtypes, by=c("gene_symbol" = "gene")) %>%
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
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              clustering_method = "ward.D2",
              annotation_row = row_annot,
              annotation_colors = annot_colors,
              fontsize_col = 4,
              fontsize_row = 6)
h <- 2 + nrow(m) / 6
w <- 2 + ncol(m) / 2
save_pheatmap_pdf(p, filename=file.path(plot_dir, "gsea_leading_edge.pdf"), width=w, height=h)


#
# GSEA Heatmap Plots
#
# x <- gsea_hallmark
# sig_cutoff <- 0.01
# sig_pathways <- filter(x, padj < sig_cutoff) %>% select(pathway) %>% distinct() %>% pull()
# x <- x %>%
#   filter(pathway %in% sig_pathways) %>%
#   mutate(sig_nes = ifelse(padj < 0.05, NES, 0)) %>%
#   # mutate(sig_nes = NES) %>%
#   select(pathway, sig_nes, analysis) %>%
#   pivot_wider(names_from = analysis, values_from = sig_nes)
# x <- column_to_rownames(x, "pathway")
# 
# brks <- c(seq(-2.5, -1.25, length.out=25), 0, seq(1.25, 2.5, length.out=25))
# p <- pheatmap(x, 
#               scale = "none",
#               color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#               breaks = brks,
#               border_color = "black",
#               cluster_rows = TRUE,
#               cluster_cols = TRUE,
#               fontsize_col = 6,
#               fontsize_row = 6)
# 
# x <- gsea_gobp
# sig_cutoff <- 0.01
# sig_pathways <- filter(x, padj < sig_cutoff) %>% select(pathway) %>% distinct() %>% pull()
# x <- x %>%
#   filter(pathway %in% sig_pathways) %>%
#   # mutate(sig_nes = ifelse(padj < 0.05, NES, 0)) %>%
#   mutate(sig_nes = NES) %>%
#   select(pathway, sig_nes, analysis) %>%
#   pivot_wider(names_from = analysis, values_from = sig_nes)
# x <- column_to_rownames(x, "pathway")
# 
# brks <- c(seq(-2.5, -1.25, length.out=25), 0, seq(1.25, 2.5, length.out=25))
# p <- pheatmap(x, 
#               scale = "none",
#               color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#               breaks = brks,
#               border_color = "black",
#               cluster_rows = TRUE,
#               cluster_cols = TRUE,
#               fontsize_col = 6,
#               fontsize_row = 6)


#########################################################
# GSEA Result Tables
#########################################################

# write gsea results
x <- bind_rows(gsea_hallmark, gsea_gobp) %>%
  filter(padj < 0.01) %>% arrange(desc(NES))
write_xlsx(list(gsea = x),
           file.path(working_dir, gsea_results_xlsx))

# # enumerate founder gene sets of HALLMARK_TNF gene set
# founder_gs <- pull(read.delim("founder_gs_tnf.txt"))
# nrow(filter(x, pathway %in% founder_gs))
# 
# # merged de analyses
# gsea_merged <- gsea %>%
#   mutate(sig = padj < gsea_padj_cutoff)
# gsea_merged <- select(gsea, pathway, padj, ES, NES, analysis) %>%
#   pivot_wider(names_from = analysis, values_from = c(padj, ES, NES)) %>%
#   mutate(nes_pb_vs_int = NES_path_pbhgd_vs_pblgd - NES_path_inthgd_vs_intlgd,
#          hist = case_when(padj_path_pbhgd_vs_pblgd < gsea_padj_cutoff & padj_path_inthgd_vs_intlgd < gsea_padj_cutoff ~ "pbint",
#                           padj_path_pbhgd_vs_pblgd < gsea_padj_cutoff ~ "pb",
#                           padj_path_inthgd_vs_intlgd < gsea_padj_cutoff ~ "int",
#                           TRUE ~ "none"))
# gsea_merged %>%
#   filter(padj_path_pbhgd_vs_pblgd < 0.01) %>%
#   arrange(desc(nes_pb_vs_int)) %>%
#   as.data.frame()
# 
# x <- gsea_merged
# 
# color_scale = c(none = "#888888",
#                 pbint = "#FF0000",
#                 pb = "#32cd32",
#                 int = "#670092")
# size_scale = c(none = 1,
#                pbint = 2,
#                pb = 2,
#                int = 2)
# 
# label_pathways = c()
# 
# annot_data <- x %>%
#   filter(hist != "none")
# none_data <- x %>%
#   filter(hist == "none")
# sig_data <- x %>%
#   filter(hist != "none")
# 
# p <- ggplot(gsea_merged, aes(x=NES_path_pbhgd_vs_pblgd, 
#                              y=NES_path_inthgd_vs_intlgd,
#                              color = hist,
#                              size = hist)) +
#   geom_point(data=none_data, alpha=0.2) + 
#   geom_point(data=sig_data, alpha=0.8) +
#   scale_color_manual(values = color_scale) +
#   scale_size_manual(values = size_scale, guide = "none") +
#   theme_minimal()
# p

#############################################################
#
# Cluster Labels
#
#############################################################




g <- filter(grade_de_merged, de_path_pbhgd_vs_pblgd != "no") %>% select(gene) %>% pull()
length(g)
g <- union(g, filter(grade_de_merged, de_path_inthgd_vs_intlgd != "no") %>% select(gene) %>% pull())
length(g)
g <- union(g, filter(grade_de_merged, de_path_hgd_vs_lgd != "no") %>% select(gene) %>% pull())
length(g)

d <- dist(t(grade_voom_lcpm[grade_genes, ]), method="euclidean")
# d <- dist(t(grade_voom_lcpm[g, ]), method="euclidean")
clust <- hclust(d, method = "ward.D2")
plot(clust, hang=-1, cex=0.5)
grade_clust_labels <- cutree(clust, k = 2) - 1

panck$samples$clust <- factor(grade_clust_labels)

design <- model.matrix(~0 + clust, data=panck$samples)
colnames(design)
contrasts = makeContrasts(clust_1_vs_0 = clust1 - clust0,
                          levels=design)
# voom with scale normalization
method <- "limma_voom_bgsub_q3"
y <- DGEList(counts=panck$bgsub_gene_mat,
             samples=panck$samples,
             group=panck$samples$clust,
             genes=gene_meta)
y <- calcNormFactors(y, method="upperquartile")
x <- limmaVoom(y, design, contrasts)
res <- processLimma(x$fit, contrasts, method, de_padj_cutoff, de_log2fc_cutoff, de_result_dir)
clust_de <- res
clust_voom_lcpm <- as_tibble(x$voom, rownames="gene_symbol") %>%
  column_to_rownames("gene_symbol")

clust_de %>% arrange(-log2fc) %>% head(20)

# gs <- gruetzmann_gs
# gs <- append(gs, hpa_panc_gs)
# gs <- append(gs, cptac_gs)
gs <- msigdb_gene_sets
gs <- split(x = gs$gene_symbol,
            f = gs$gs_name)


method <- "limma_voom_bgsub_q3"
analysis <- "clust_1_vs_0"
ranks <- get_de_ranks(clust_de, analysis)
res <- fgsea(pathways = gs, 
             stats = ranks,
             minSize = 10,
             eps = 0,
             nPermSimple = 10000)
res <- mutate(res, method = method, analysis = analysis)


filter(res, padj < 0.05) %>% arrange(desc(NES)) %>% select(pathway, pval, padj, ES, NES, size) %>% as.data.frame()

# write gsea results
f <- file.path(working_dir, "gsea_clust_results.tsv")
data.table::fwrite(res, file=f, sep="\t", sep2=c("", " ", ""))
f <- file.path(working_dir, "gsea_clust_results.xlsx")
write_xlsx(res, f)


rf <- randomForest(x=t(grade_voom_lcpm[g, ]), 
                   y=factor(grade_clust_labels),
                   ntree = 10000,
                   importance = TRUE)
rf
varImpPlot(rf)


names(cutree(clust, k = 2)) == colnames(panck$bgsub_gene_mat)
rect.hclust(clust, k=2, border="red")

# Ward Hierarchical Clustering
d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward")
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
rect.hclust(fit, k=5, border="red")

?hclust
?pheatmap
clustering_distance_rows = "euclidean",
clustering_distance_cols = "euclidean", clustering_method = "complete",
clustering_callback = identity2, cutree_rows = NA, cutree_cols = NA,
treeheight_row = ifelse((class(cluster_rows) == "hclust") || cluster_rows,
                        50, 0), treeheight_col = ifelse((class(cluster_cols) == "hclust") ||
                                                          cluster_cols, 50, 0)
hclust


#############################################################
#
# Figure data for manuscript
#
#############################################################

# MUC1 expression associated with histology
as.data.frame(hist_de_merged %>% filter(gene == "CLU"))
as.data.frame(grade_de_merged %>% filter(gene == "SMAD4"))
as.data.frame(filter(grade_de_merged, gene == "CXCL8"))
as.data.frame(filter(grade_de_merged, gene == "SERPINB5"))
as.data.frame(filter(grade_de_merged, gene == "LAMB3"))
as.data.frame(grade_de_merged %>% filter(gene == "SMAD4"))
as.data.frame(grade_de_merged %>% filter(gene == "LIF"))


filter(grade_de_merged, de_path_pbhgd_vs_lgd != "no") %>% select(gene) %>% pull()
as.data.frame(filter(grade_de_merged, de_path_inthgd_vs_lgd != "no"))
as.data.frame(filter(grade_de_merged, de_inthgd != "no"))
colnames(grade_de_merged)


#############################################################
#
# PanCK / CD45 / SMA correlation analysis
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
# extract pairs of colocalized regions
# get normalized count matrices
#
# panck cd45
panck_cd45_meta <- samples %>%
  filter(segment == "panck" | segment == "cd45") %>%
  group_by(slide_id, roi_group) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  mutate(pair_id = paste0("s_", slide_id, "_", roi_group)) %>%
  arrange(pair_id)
table(panck_cd45_meta$segment)

paired_panck_cd45 <- list()
paired_panck_cd45$panck <- subset_paired_segment(panck_cd45_meta, voom_bgsub_q3, "panck")
paired_panck_cd45$cd45 <- subset_paired_segment(panck_cd45_meta, voom_bgsub_q3, "cd45")

# panck sma
panck_sma_meta <- samples %>%
  filter(segment == "panck" | segment == "sma") %>%
  group_by(slide_id, roi_group) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  mutate(pair_id = paste0("s_", slide_id, "_", roi_group)) %>%
  arrange(pair_id)
table(panck_sma_meta$segment)

paired_panck_sma <- list()
paired_panck_sma$panck <- subset_paired_segment(panck_sma_meta, voom_bgsub_q3, "panck")
paired_panck_sma$sma <- subset_paired_segment(panck_sma_meta, voom_bgsub_q3, "sma")


#
# correlation analysis
#
nperms <- 10000
cor_qval_cutoff <- 0.01
cor_r_cutoff <- 0.7
cor_weight_beta <- 4
graph_min_csize <- 10
graph_clust_resolution <- 0.25

# panck-cd45 correlation
# paired_panck_cd45$cor_result <- cor_compute(paired_panck_cd45$panck$mat,
#                                             paired_panck_cd45$cd45$mat,
#                                             nperms)
# saveRDS(paired_panck_cd45$cor_result, file=file.path(working_dir, "cor_panck_cd45.rds"))
paired_panck_cd45$cor_result <- readRDS(file.path(working_dir, "cor_panck_cd45.rds"))

# panck-sma correlation
# paired_panck_sma$cor_result <- cor_compute(paired_panck_sma$panck$mat,
#                                            paired_panck_sma$sma$mat,
#                                            nperms)
# saveRDS(paired_panck_sma$cor_result, file=file.path(working_dir, "cor_panck_sma.rds"))
paired_panck_sma$cor_result <- readRDS(file.path(working_dir, "cor_panck_sma.rds"))

# panck self correlation
panck_cor_mat <- grade_voom_lcpm
rownames(panck_cor_mat) <- paste0("panck_", rownames(grade_voom_lcpm))
# panck_cor_res <- cor_compute(panck_cor_mat, 
#                              panck_cor_mat, 
#                              nperms)
# saveRDS(panck_cor_res, file=file.path(working_dir, "cor_panck.rds"))
panck_cor_res <- readRDS(file.path(working_dir, "cor_panck.rds"))

# cd45 self correlation
cd45_cor_mat <- cd45_voom_lcpm
rownames(cd45_cor_mat) <- paste0("cd45_", rownames(cd45_voom_lcpm))
# cd45_cor_res <- cor_compute(cd45_cor_mat, 
#                             cd45_cor_mat, 
#                             nperms)
# saveRDS(cd45_cor_res, file=file.path(working_dir, "cor_cd45.rds"))
cd45_cor_res <- readRDS(file.path(working_dir, "cor_cd45.rds"))

# sma self correlation
sma_cor_mat <- sma_voom_lcpm
rownames(sma_cor_mat) <- paste0("sma_", rownames(sma_voom_lcpm))
# sma_cor_res <- cor_compute(sma_cor_mat, 
#                            sma_cor_mat, 
#                            nperms)
# saveRDS(sma_cor_res, file=file.path(working_dir, "cor_sma.rds"))
sma_cor_res <- readRDS(file.path(working_dir, "cor_sma.rds"))



#
# construct network
#
# edges
edges_all <- bind_rows(cor_result_to_edges(paired_panck_cd45$cor_result),
                       cor_result_to_edges(paired_panck_sma$cor_result),
                       cor_result_to_edges(panck_cor_res, self=TRUE),
                       cor_result_to_edges(cd45_cor_res, self=TRUE),
                       cor_result_to_edges(sma_cor_res, self=TRUE))

# filter significant edges
edges <- filter(edges_all, qval < cor_qval_cutoff,
                cor >= cor_r_cutoff)
edges$weight <- abs(edges$cor) ** cor_weight_beta
# nodes
nodes_panck <- bind_cols(label = rownames(panck_cor_mat),
                         orig_gene = rownames(grade_voom_lcpm),
                         type = "panck",
                         grade_subtype = grade_de_merged$subtype,
                         grade_log2fc = grade_de_merged$log2fc_path_hgd_vs_lgd)
nodes_cd45 <- bind_cols(label = rownames(cd45_cor_mat),
                        orig_gene = rownames(cd45_voom_lcpm),
                        type = "cd45",
                        grade_subtype = "cd45",
                        grade_log2fc = 0)
nodes_sma <- bind_cols(label = rownames(sma_cor_mat),
                       orig_gene = rownames(sma_voom_lcpm),
                       type = "sma",
                       grade_subtype = "sma",
                       grade_log2fc = 0)
nodes <- bind_rows(nodes_panck, nodes_cd45, nodes_sma)
nodes$id <- nodes$label
# filter nodes
nodes <- nodes %>% filter(id %in% union(edges$source, edges$target))


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
                   orig_gene = rownames(grade_voom_lcpm),
                   type = "panck",
                   grade_subtype = grade_de_merged$subtype,
                   grade_log2fc = grade_de_merged$log2fc_path_hgd_vs_lgd,
                   grade_pb_log2fc = grade_de_merged$log2fc_path_pbhgd_vs_pblgd,
                   grade_int_log2fc = grade_de_merged$log2fc_path_inthgd_vs_intlgd,
                   grade_abs_log2fc = abs(grade_de_merged$log2fc_path_hgd_vs_lgd),
                   grade_abs_pb_log2fc = abs(grade_de_merged$log2fc_path_pbhgd_vs_pblgd),
                   grade_abs_int_log2fc = abs(grade_de_merged$log2fc_path_inthgd_vs_intlgd))
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
ggplot(res, aes(x=r, y=m)) + 
  geom_point()
ggplot(res, aes(x=r, y=n)) + 
  geom_point()

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

# TODO: explore 'nest' and 'unnest'
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
  if (length(ni) < min_csize) {
    next
  }
  for (j in (i+1):nclust) {
    nj <- clust_nodes %>% dplyr::filter(community == j) %>% dplyr::select(id) %>% dplyr::pull()
    if (length(nj) < min_csize) {
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

# plot(clust, g,
#      layout=layout_with_kk,
#      vertex.size=2,
#      vertex.label=NA)

gs <- msigdb_gene_sets %>%
  filter(gs_cat == "H") %>%
  select(gs_name, gene_symbol)
gs <- msigdb_gene_sets %>%
  dplyr::filter(gs_cat == "C5", gs_subcat == "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)

gs <- msigdb_gene_sets %>%
  dplyr::filter(gs_cat == "C7", gs_subcat == "IMMUNESIGDB") %>%
  dplyr::select(gs_name, gene_symbol)

gs <- msigdb_gene_sets %>%
  dplyr::select(gs_name, gene_symbol)


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

gs <- gmt_to_tbl(cptac_gs)
gs <- gmt_to_tbl(gmtPathways(gs_wikipathways_file))
gs <- gmt_to_tbl(gmtPathways(gs_reactome_file))



seg <- "panck"
yy <- NULL
for (i in 1:length(unique(clust_nodes$community))) {
  x <- dplyr::filter(clust_nodes, type == seg, community == i) %>% 
    dplyr::select(orig_gene) %>% pull()
  if (length(x) < 10) {
    yy[[i]] <- NULL
  } else {
    yy[[i]] <- enricher(gene = x,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        universe = cta_genes,
                        qvalueCutoff = 0.1,
                        TERM2GENE = gs)
    print(yy[[i]])
  }
}

yy
edo2 <- gseDO(geneList)


dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")

?enrichplot::dotplot


enrichplot::dotplot(yy[[10]], showCategory = 3, font.size = 6)


str(yy)
str(yy[[2]])
yy[[8]]
yy[[2]]@result$qvalue
yy[[16]]
enrichplot::dotplot(yy[[16]])



enricher(gene = x,
         pvalueCutoff = 0.2,
         pAdjustMethod = "BH",
         universe = cta_genes,
         TERM2GENE = gs)

enrichGO(gene = x,
         universe = cta_genes,
         OrgDb = "org.Hs.eg.db",
         keyType = "SYMBOL",
         ont = "BP",
         pAdjustMethod = "BH",
         pvalueCutoff  = 0.2,
         qvalueCutoff  = 0.2,
         minGSSize = 10,
         maxGSSize = 500,
         readable = TRUE)




y <- dplyr::filter(gene_meta, gene_symbol %in% x) %>% dplyr::select(entrez_id) %>% dplyr::pull()
y <- filter(gene_meta, gene_symbol %in% x) %>% dplyr::select(entrez_id) %>% dplyr::pull()
enricher(x, 
         pvalueCutoff = 0.05,
         pAdjustMethod = "BH",
         universe = cta_genes,
         TERM2GENE = dplyr::select(
           hs_kegg_df,
           gs_name,
           human_entrez_gene
         )
)






e <- cor_result_to_edges(panck_cor_res, self=TRUE)
table(e$qval < 0.01)
e %>% arrange(desc(cor))
plot(grade_voom_lcpm["ACTB",])
plot(t(grade_voom_lcpm["LIF",]), t(grade_voom_lcpm["EPHA2",]), col=color_scales$path[panck$samples$path])
cor.test(t(grade_voom_lcpm["LIF",]), t(grade_voom_lcpm["CXCL8",]), method="pearson")

# find connected components (at least 10 vertices)
# TODO: improve this
# min_csize <- 10
# gcomponents <- decompose(g, mode="weak", min.vertices=min_csize)
# gcomponents
# sort(table(components(g)$membership))
# components(g)$csize



# edges_panck <- cor_to_edgelist(cor_panck, self=TRUE)
# edges_panck$padj <- p.adjust(edges_panck$pval, method = "fdr")
# qobj <- qvalue(edges_panck$pval)
# edges_panck$qval <- qobj$qvalues
# 
# edges_cd45 <- cor_to_edgelist(cor_cd45, self=TRUE)
# edges_cd45$padj <- p.adjust(edges_cd45$pval, method = "fdr")
# qobj <- qvalue(edges_cd45$pval)
# edges_cd45$qval <- qobj$qvalues


#
# process correlation results
#
edges_panck <- cor_to_edgelist(cor_panck, self=TRUE)
edges_panck$padj <- p.adjust(edges_panck$pval, method = "fdr")
qobj <- qvalue(edges_panck$pval)
edges_panck$qval <- qobj$qvalues

edges_cd45 <- cor_to_edgelist(cor_cd45, self=TRUE)
edges_cd45$padj <- p.adjust(edges_cd45$pval, method = "fdr")
qobj <- qvalue(edges_cd45$pval)
edges_cd45$qval <- qobj$qvalues

edges_panck_cd45 <- cor_to_edgelist(cor_panck_cd45, self=FALSE)
edges_panck_cd45$padj <- p.adjust(edges_panck_cd45$pval, method = "fdr")
qobj <- qvalue(edges_panck_cd45$pval)
edges_panck_cd45$qval <- qobj$qvalues
edges_panck_cd45$weight <- abs(edges_panck_cd45$cor)

# select panck-cd45 edges
edges <- edges_panck_cd45 %>% filter(qval < 0.01)
nodes <- nodes %>% filter(id %in% union(edges$source, edges$target))

# build graph object
g <- igraph::graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)

# find connected components (at least 10 vertices)
# TODO: improve this
min_csize <- 10
gcomponents <- decompose(g, mode="weak", min.vertices=min_csize)
table(components(g)$membership)
components(g)$csize
g <- gcomponents[[1]]
# choose largest graph with most components
# g <- igraph::induced_subgraph(g, igraph::V(g)[igraph::components(g)$membership == which.max(igraph::components(g)$csize)])

# clustering/community detection
clust <- cluster_leiden(g, 
                        objective_function="modularity",
                        resolution=0.1,
                        n_iterations = 100)
table(membership(clust))
modularity(g, membership(clust))
plot(clust, g, 
     layout=layout_with_kk, 
     vertex.size=2,
     vertex.label=NA)
V(g)$community <- factor(clust$membership)

clust_nodes <- as_tibble(as_data_frame(g, what="vertices"))
clust_edges <- as_tibble(as_data_frame(g, what="edges"))

clust_nodes <- clust_nodes %>%
  rename(label = name)
clust_edges <- clust_edges %>%
  rename(source = from,
         target = to)
clust_edges$cordir = ifelse(clust_edges$cor >= 0, "pos", "neg")


write.table(clust_nodes, file=file.path(working_dir, "gephi_clust_panck_cd45_nodes.tsv"),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
write.table(clust_edges, file=file.path(working_dir, "gephi_clust_panck_cd45_edges.tsv"),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)


clust_nodes <- as_tibble(as_data_frame(g, what="vertices"))
clust_nodes[1,]
for (i in 1:clust$nb_clusters) {
  comm_nodes <- filter(clust_nodes, community == i)
  panck_genes <- comm_nodes %>% filter(type == "panck") %>% select(orig_gene) %>% pull()
  cd45_genes <- comm_nodes %>% filter(type == "cd45") %>% select(orig_gene) %>% pull()
  if (length(panck_genes) == 0) {
    next
  }
  if (length(cd45_genes) == 0) {
    next
  }
  f = file.path(cor_clust_dir, paste0("clust", i, "_panck.txt"))
  write(panck_genes, file=f, sep="\n")
  f = file.path(cor_clust_dir, paste0("clust", i, "_cd45.txt"))
  write(cd45_genes, file=f, sep="\n")
}


# filter significant edges
edges <- bind_rows(edges_panck, edges_cd45, edges_panck_cd45)
hist(edges$qval)
edges$weight <- abs(edges$cor)
edges <- edges %>% filter(qval < 0.01)
nodes <- nodes %>% filter(gene %in% union(edges$source, edges$target))

# build graph object
g <- igraph::graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)

# find connected components (at least 10 vertices)
# TODO: improve this
min_csize <- 10
gcomponents <- decompose(g, mode="weak", min.vertices=min_csize)
g <- gcomponents[[1]]
# choose largest graph with most components
# g <- igraph::induced_subgraph(g, igraph::V(g)[igraph::components(g)$membership == which.max(igraph::components(g)$csize)])


clust_nodes <- as_tibble(as_data_frame(g, what="vertices"))
clust_edges <- as_tibble(as_data_frame(g, what="edges"))


# TODO: components with all types (panck and cd45)
E(g)$weight <- abs(E(g)$cor) ** 4

# clustering/community detection
leidclust <- cluster_leiden(g, 
                            objective_function="modularity",
                            resolution=1)
clust <- leidclust


modularity(g, membership(clust))
V(g)$community <- factor(clust$membership)
table(membership(clust))

clust_nodes <- as_tibble(as_data_frame(g, what="vertices"))
clust_nodes[1,]
for (i in 1:clust$nb_clusters) {
  comm_nodes <- filter(clust_nodes, community == i)
  panck_genes <- comm_nodes %>% filter(type == "panck") %>% select(orig_gene) %>% pull()
  cd45_genes <- comm_nodes %>% filter(type == "cd45") %>% select(orig_gene) %>% pull()
  if (length(panck_genes) == 0) {
    next
  }
  if (length(cd45_genes) == 0) {
    next
  }
  f = file.path(cor_clust_dir, paste0("clust", i, "_panck.txt"))
  write(panck_genes, file=f, sep="\n")
  f = file.path(cor_clust_dir, paste0("clust", i, "_cd45.txt"))
  write(cd45_genes, file=f, sep="\n")
}




?cluster_leiden
igraph::modularity(com)
igraph::membership(leidclust)

print(leidclust)
?membership


ebc <- cluster_edge_betweenness(g, weights=1-E(g)$weight)
imc <- cluster_infomap(g)
lec <- cluster_leading_eigen(g)
loc <- cluster_louvain(g)
leidc <- cluster_leiden(g, objective_function="modularity")
sgc <- cluster_spinglass(g)
wtc <- cluster_walktrap(g)
scores <- c(ebc = modularity(g, membership(ebc)),
            infomap = modularity(g,membership(imc)),
            eigen = modularity(g,membership(lec)),
            louvain = modularity(g,membership(loc)),
            leiden = modularity(g,membership(leidc)),
            spinglass = modularity(g,membership(sgc)),
            walk = modularity(g,membership(wtc)))


cluster_leiden(g)

com <- igraph::cluster_louvain(g, resolution = 0.25)
l <- igraph::layout_with_fr(g)



?igraph::communities
com <- igraph::cluster_edge_betweenness(g)

V(g)$color <- com$membership+1
g <- set_graph_attr(g, "layout", layout_with_kk(g))
plot(g, vertex.label.dist=1.5)

?igraph::edge.betweenness.community

plot(g, 
     layout=l,
     vertex.size=1, 
     vertex.label.cex=0.5,
     vertex.frame.color=NA)

plot(g, layout=l, edge.color="black", vertex.label="", main=layout)





write.table(n, file=file.path(working_dir, "graph_panck_cd45_nodes.tsv"),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
write.table(e, file=file.path(working_dir, "graph_panck_cd45_edges.tsv"),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

color_scale = c(none = "#888888",
                hgd = "#FF0000",
                lgd = "#0000FF",
                pbhgd = "#32cd32",
                inthgd = "#670092",
                pblgd = "#1677c4",
                intlgd = "#ff6700",
                cd45 = "#ffff00")
V(g)$color <- color_scale[V(g)$grade_subtype]

# E(g)$weight <- E(g)$cor
E(g)$weight <- 500*abs(E(g)$cor)**6

l = layout_with_fr(g)
plot(g, layout=l, 
     vertex.size=3, 
     vertex.label=vertex_attr(g, "orig_gene"),
     vertex.label.family="Arial",
     vertex.label.cex=0.2,
     vertex.label.color="black")


E(g)$cor
edge_attr(g)
vertex_attr(g, "orig_gene")
plot(g4, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5,
     vertex.color=c( "pink", "skyblue")[1+(V(g4)$gender=="male")] ) 






library(RColorBrewer)
colors=brewer.pal(length(com),'Accent') #make a color palette
V(g.sparrow)$color=colors[membership(com)] #assign each vertex a color based on the community assignment

set.seed(2)
plot(g.sparrow, vertex.label="", edge.width=E(g.sparrow)$weight*5)


panck_grade_nodes <- filter(nodes, grade_subtype != "none", grade_subtype != "cd45")
cd45_nodes <- filter(nodes, type == "cd45")
e <- edges %>% 
  filter(qval < 0.01,
         from %in% panck_grade_nodes$gene,
         to %in% cd45_nodes$gene)

rows <- unique(e$from)
cols <- unique(e$to)
m <- cor_panck_cd45$r[rows, cols]
m <- sign(m) * m**4

p <- pheatmap(m, 
              scale = "none",
              color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
              breaks = seq(-0.9, 0.9, length.out=51),
              border_color = NA,
              cluster_rows = TRUE,
              cluster_cols = TRUE,
              clustering_method = "ward.D2",
              fontsize_row = 6,
              fontsize_col = 6)

save_pheatmap_pdf(p, filename=file.path(plot_dir, "cor_panck_cd45.pdf"), width=10, height=10)



#
# filter edges to enable network creation
#
e <- bind_rows(edges_panck, edges_cd45, edges_panck_cd45)
e$weight <- abs(e$cor) ** 3
e <- e %>% 
  filter(qval < 0.01)


ranks <- cor_tbl %>%
  filter(panck == "LIF") %>%
  select(cd45, cor) %>%
  arrange(desc(cor))
ranks = sort(setNames(ranks$cor, ranks$cd45), decreasing = TRUE)

res <- fgsea(pathways = gs_gobp,
             stats = ranks,
             minSize = 10,
             eps = 0,
             nPermSimple = 10000)
res <- as_tibble(res) %>%
  filter(padj < 0.05) %>% arrange(desc(NES))
res$pathway

get_de_ranks <- function(de, a) {
  ranks <- de %>% 
    filter(analysis == a) %>%
    select(gene, log2fc, padj) %>%
    mutate(rank = log2fc)
  ranks = sort(setNames(ranks$rank, ranks$gene), decreasing = TRUE)
  return(ranks)
}

filter(res, padj < 0.05) %>% arrange(desc(NES)) %>% select(pathway, pval, padj, ES, NES, size) %>% as.data.frame()

dim(x)
x <- cor_tbl %>%
  filter(pval < 0.05)
ggplot(x, aes(x = cor, y = -log10(padj))) + 
  geom_point()

x %>% arrange(desc(cor)) %>% head(30) %>% as.data.frame()
x %>% filter(cd45 == "ITGA2" | panck == "ITGA2", padj < 0.05) %>% arrange(cor) %>% as.data.frame()
x %>% filter(panck == "CD55", pval < 0.05) %>% arrange(cor) %>% as.data.frame()



e <- edges_panck_cd45 %>%
  mutate(source_orig_gene = gsub("panck_", "", source),
         target_orig_gene = gsub("cd45_", "", target)) %>%
  mutate(same_gene = source_orig_gene == target_orig_gene)

ggplot(e, aes(x=cor, y=factor(same_gene))) + 
  stat_density_ridges(alpha=0.7, scale=2, rel_min_height=0.0, quantile_lines=TRUE)


plot(paired_panck$mat["panck_LIF",], paired_cd45$mat["cd45_LIF",])



# # filter genes associated with disease (include cd45 nodes)
# nrow(nodes)
# n <- filter(nodes, grade_subtype != "none")
# # get edges
# nrow(e)
# e <- e %>% filter((from %in% n$gene) | (to %in% n$gene))
# nrow(e)
# subset nodes
# sum(nodes$gene %in% union(e$from, e$to))
# n <- filter(nodes, gene %in% union(e$from, e$to))

# 
# paired_cor_test_v1 <- function(x, y, nperms, method = "spearman") {
#   nrows <- nrow(x)
#   ncols <- ncol(x)
#   r <- array(data = 0, dim = c(ncols, ncols), dimnames = list(colnames(x), colnames(y)))
#   pvals <- array(data = 1, dim = c(ncols, ncols), dimnames = list(colnames(x), colnames(y)))
#   set.seed(1)
#   perms <- replicate(nperms, sample(nrows))
#   
#   for (i in 1:ncols) {
#     r[i, ] <- cor(x[, i], y, method = method)
#     rnull <- array(0, dim = c(nperms, ncols))
#     for (p in 1:nperms) {
#       rnull[p, ] <- abs(cor(x[, i], y[perms[, p], ], method = method))
#     }
#     pvals_row <- sweep(rnull, 2, FUN=">=", abs(r[i, ]))
#     pvals[i, ] <- apply(pvals_row, 2, mean)
#     print(i)
#   }
#   return(list(r = r, pvals = pvals))
# }




x <- paired_panck_cd45$panck$mat[paste0("panck_", c("LIF", "CXCL8")),]
y <- paired_panck_cd45$cd45$mat
r <- cor_compute(x, y, 10000)
edges <- cor_result_to_edges(r)


n <- bind_cols(label = rownames(paired_panck_cd45$cd45$mat),
               orig_gene = paired_panck_cd45$cd45$genes,
               type = "cd45")

e <- filter(edges, source == "panck_LIF") %>%
  inner_join(n, by=c("target"="label"))

x <- e %>% filter(cor < 0, pval < 0.01) %>% select(orig_gene) %>% pull()
write(x, file=file.path(working_dir, "nodes_dn.txt"), sep="\n")
x <- e %>% filter(cor >= 0, pval < 0.01) %>% select(orig_gene) %>% pull()
write(x, file=file.path(working_dir, "nodes_up.txt"), sep="\n")



ranks <- cor_tbl %>%
  filter(panck == "LIF") %>%
  select(cd45, cor) %>%
  arrange(desc(cor))
ranks = sort(setNames(ranks$cor, ranks$cd45), decreasing = TRUE)

res <- fgsea(pathways = gs_gobp,
             stats = ranks,
             minSize = 10,
             eps = 0,
             nPermSimple = 10000)
res <- as_tibble(res) %>%
  filter(padj < 0.05) %>% arrange(desc(NES))
res$pathway

get_de_ranks <- function(de, a) {
  ranks <- de %>% 
    filter(analysis == a) %>%
    select(gene, log2fc, padj) %>%
    mutate(rank = log2fc)
  ranks = sort(setNames(ranks$rank, ranks$gene), decreasing = TRUE)
  return(ranks)
}

filter(res, padj < 0.05) %>% arrange(desc(NES)) %>% select(pathway, pval, padj, ES, NES, size) %>% as.data.frame()

dim(x)
x <- cor_tbl %>%
  filter(pval < 0.05)
ggplot(x, aes(x = cor, y = -log10(padj))) + 
  geom_point()

x %>% arrange(desc(cor)) %>% head(30) %>% as.data.frame()
x %>% filter(cd45 == "ITGA2" | panck == "ITGA2", padj < 0.05) %>% arrange(cor) %>% as.data.frame()
x %>% filter(panck == "CD55", pval < 0.05) %>% arrange(cor) %>% as.data.frame()



# #
# # full network (panck + sma + cd45)
# #
# # edges
# edges_all <- bind_rows(cor_result_to_edges(paired_panck_cd45$cor_result),
#                        cor_result_to_edges(paired_panck_sma$cor_result),
#                        cor_result_to_edges(panck_cor_res, self=TRUE),
#                        cor_result_to_edges(cd45_cor_res, self=TRUE),
#                        cor_result_to_edges(sma_cor_res, self=TRUE))
# 
# # filter significant edges
# # edges <- filter(edges_all, qval < cor_qval_cutoff,
# #                 abs(cor) >= cor_r_cutoff)
# edges <- filter(edges_all, qval < cor_qval_cutoff,
#                 cor >= cor_r_cutoff)
# edges$weight <- abs(edges$cor) ** cor_weight_beta
# #edges$cordir <- ifelse(edges$cor >= 0, "pos", "neg")
# #edges$weight <- edges$cor ** cor_weight_beta
# #edges$weight <- ((1 + edges$cor)/2) ** cor_weight_beta
# #edges$weight <- abs(edges$cor) ** cor_weight_beta
# #edges$weight <- sign(edges$cor) * abs(edges$cor) ** cor_weight_beta
# 
# 
# # nodes
# nodes_panck <- bind_cols(label = rownames(panck_cor_mat),
#                          orig_gene = rownames(grade_voom_lcpm),
#                          type = "panck",
#                          grade_subtype = grade_de_merged$subtype,
#                          grade_log2fc = grade_de_merged$log2fc_path_hgd_vs_lgd)
# nodes_cd45 <- bind_cols(label = rownames(cd45_cor_mat),
#                         orig_gene = rownames(cd45_voom_lcpm),
#                         type = "cd45",
#                         grade_subtype = "cd45",
#                         grade_log2fc = 0)
# nodes_sma <- bind_cols(label = rownames(sma_cor_mat),
#                        orig_gene = rownames(sma_voom_lcpm),
#                        type = "sma",
#                        grade_subtype = "sma",
#                        grade_log2fc = 0)
# nodes <- bind_rows(nodes_panck, nodes_cd45, nodes_sma)
# nodes$id <- nodes$label
# 
# # filter nodes
# nodes <- nodes %>% filter(id %in% union(edges$source, edges$target))
# 
# # build graph object
# g <- igraph::graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)
# 
# # choose largest weakly connected graph
# # g <- igraph::induced_subgraph(g, igraph::V(g)[igraph::components(g)$membership == which.max(igraph::components(g)$csize)])
# 
# # # plot graph
# # plot(g,
# #      layout=layout_with_fr,
# #      vertex.size=2,
# #      vertex.label=NA)
# 
# # clustering/community detection
# clust <- cluster_leiden(g,
#                         objective_function="modularity",
#                         resolution=graph_clust_resolution,
#                         n_iterations=1000)
# modularity(g, membership(clust))
# table(membership(clust))
# # plot(clust, g,
# #      layout=layout_with_kk,
# #      vertex.size=2,
# #      vertex.label=NA)
# V(g)$community <- factor(clust$membership)
# 
# clust_nodes <- as_tibble(as_data_frame(g, what="vertices"))
# clust_edges <- as_tibble(as_data_frame(g, what="edges"))
# clust_nodes <- clust_nodes %>%
#   dplyr::rename(label = name)
# clust_edges <- clust_edges %>%
#   dplyr::rename(source = from,
#                 target = to)
# 
# write.table(clust_nodes, file=file.path(working_dir, "graph_clust_nodes.tsv"),
#             quote=FALSE,
#             sep="\t",
#             row.names=FALSE)
# write.table(clust_edges, file=file.path(working_dir, "graph_clust_edges.tsv"),
#             quote=FALSE,
#             sep="\t",
#             row.names=FALSE)
# 
# 
# 
# # TODO: explore 'nest' and 'unnest'
# merged_nodes <- clust_nodes %>%
#   dplyr::group_by(community) %>%
#   dplyr::summarize(size = n()) %>%
#   dplyr::mutate(id = paste0("c", community),
#                 label = id)
# 
# nclust <- length(unique(clust_nodes$community))
# e <- filter(edges_all, qval < cor_qval_cutoff)
# 
# merged_edges <- NULL
# for (i in 1:nclust) {
#   ni <- clust_nodes %>% dplyr::filter(community == i) %>% dplyr::select(id) %>% dplyr::pull()
#   if (length(ni) < min_csize) {
#     next
#   }
#   for (j in (i+1):nclust) {
#     nj <- clust_nodes %>% dplyr::filter(community == j) %>% dplyr::select(id) %>% dplyr::pull()
#     if (length(nj) < min_csize) {
#       next
#     }
#     eij <- bind_rows(edges_all %>% filter(source %in% ni, target %in% nj),
#                      edges_all %>% filter(source %in% nj, target %in% ni))
#     if (nrow(eij) == 0) {
#       next
#     }
#     avgcor <- mean(eij$cor)
#     tbl <- tibble(source = paste0("c", i), 
#                   target = paste0("c", j), 
#                   avgcor = avgcor,
#                   n = nrow(eij))
#     merged_edges <- bind_rows(merged_edges, tbl)
#     print(paste0(i, " d ", j))
#   }
# }
# 
# summary(merged_edges$avgcor)
# merged_edges <- merged_edges %>%
#   mutate(weight = 1 - avgcor,
#          cordir = ifelse(avgcor < 0, "neg", "pos"))
# 
# 
# write.table(merged_nodes, file=file.path(working_dir, "graph_merged_clust_nodes.tsv"),
#             quote=FALSE,
#             sep="\t",
#             row.names=FALSE)
# write.table(merged_edges, file=file.path(working_dir, "graph_merged_clust_edges.tsv"),
#             quote=FALSE,
#             sep="\t",
#             row.names=FALSE)
# 
# mg <- igraph::graph_from_data_frame(d=merged_edges, vertices=merged_nodes, directed=FALSE)
# plot(mg, 
#      vertex.size = 2)
# 
# merged_nodes[1,]
# # plot(clust, g,
# #      layout=layout_with_kk,
# #      vertex.size=2,
# #      vertex.label=NA)
