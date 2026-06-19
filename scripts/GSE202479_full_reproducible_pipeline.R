#!/usr/bin/env Rscript

# Reproduce the final GSE202479 normal vs gallstones analysis:
# 1) download GEO processed raw-count matrix and MSigDB GMT
# 2) extract normal/gallstones raw counts
# 3) compute edgeR CPM and DESeq2 differential expression
# 4) run clusterProfiler GSEA against full MSigDB plus custom marker sets
# 5) draw enrichplot GSEA plots for IL17/IL22, IL17, TCR, and marker gene sets

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)
  project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
} else {
  project_dir <- normalizePath(getwd(), mustWork = FALSE)
}
if (!dir.exists(file.path(project_dir, "data"))) {
  project_dir <- normalizePath(getwd(), mustWork = FALSE)
}
setwd(project_dir)

required_pkgs <- c(
  "data.table", "edgeR", "DESeq2", "clusterProfiler",
  "enrichplot", "ggplot2", "gridExtra"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install them in the bioinfo environment before running this script."
  )
}

suppressPackageStartupMessages(suppressWarnings({
  library(data.table)
  library(edgeR)
  library(DESeq2)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
}))

force_download <- identical(Sys.getenv("FORCE_DOWNLOAD"), "1")
force_deseq <- identical(Sys.getenv("FORCE_DESEQ2"), "1")
force_gsea <- identical(Sys.getenv("FORCE_CLUSTERPROFILER_GSEA"), "1")

raw_dir <- "data/raw"
processed_dir <- "data/processed"
deseq_dir <- "results/deseq2"
gsea_dir <- "results/clusterprofiler_enrichplot_full_msigdb_v2026_1"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deseq_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)

geo_expression_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202479/suppl/",
  "GSE202479_gene_expression_anno.txt.gz"
)
geo_series_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202479/matrix/",
  "GSE202479_series_matrix.txt.gz"
)
msigdb_url <- paste0(
  "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/",
  "msigdb.v2026.1.Hs.symbols.gmt"
)

geo_expression_file <- file.path(raw_dir, "GSE202479_gene_expression_anno.txt.gz")
geo_series_file <- file.path(raw_dir, "GSE202479_series_matrix.txt.gz")
msigdb_file <- file.path(raw_dir, "msigdb.v2026.1.Hs.symbols.gmt")

download_if_needed <- function(url, dest, force = FALSE) {
  if (file.exists(dest) && file.info(dest)$size > 0 && !force) {
    message("Using cached file: ", dest)
    return(invisible(dest))
  }
  message("Downloading: ", url)
  tmp <- paste0(dest, ".tmp")
  if (file.exists(tmp)) unlink(tmp)
  utils::download.file(url, tmp, mode = "wb", method = "libcurl", quiet = FALSE)
  if (!file.exists(tmp) || file.info(tmp)$size == 0) {
    stop("Download failed or produced an empty file: ", url)
  }
  file.rename(tmp, dest)
  invisible(dest)
}

message("Step 1: downloading inputs")
download_if_needed(geo_expression_url, geo_expression_file, force_download)
download_if_needed(geo_series_url, geo_series_file, force_download)
if (!file.exists(msigdb_file) || file.info(msigdb_file)$size == 0 || force_download) {
  local_msigdb <- file.path(Sys.getenv("HOME"), "Downloads", "msigdb.v2026.1.Hs.symbols.gmt.txt")
  if (!force_download && file.exists(local_msigdb) && file.info(local_msigdb)$size > 0) {
    message("Copying local MSigDB file from Downloads: ", local_msigdb)
    file.copy(local_msigdb, msigdb_file, overwrite = TRUE)
  } else {
    download_if_needed(msigdb_url, msigdb_file, force_download)
  }
} else {
  message("Using cached file: ", msigdb_file)
}

samples <- data.frame(
  accession = c(
    "GSM6123008", "GSM6123009", "GSM6123010",
    "GSM6123011", "GSM6123012", "GSM6123013", "GSM6123014"
  ),
  sample = c("N8", "N10", "N20", "Y8", "Y12", "Y13", "Y16"),
  group = c(rep("Normal", 3), rep("Gallstones", 4)),
  count_column = c(
    "N8_count", "N10_count", "N20_count",
    "Y8_count", "Y12_count", "Y13_count", "Y16_count"
  ),
  title = c(
    "normal gallbladder [N8]",
    "normal gallbladder [N10]",
    "normal gallbladder [N20]",
    "gallbladder with chronic inflammation induced by gallstones [Y8]",
    "gallbladder with chronic inflammation induced by gallstones [Y12]",
    "gallbladder with chronic inflammation induced by gallstones [Y13]",
    "gallbladder with chronic inflammation induced by gallstones [Y16]"
  )
)

counts_file <- file.path(processed_dir, "GSE202479_normal_gallstones_raw_counts.csv")
metadata_file <- file.path(processed_dir, "GSE202479_normal_gallstones_sample_metadata.csv")
metadata_detailed_file <- file.path(processed_dir, "GSE202479_normal_gallstones_sample_metadata_detailed.csv")

message("Step 2: extracting normal/gallstones raw count matrix")
raw_anno <- fread(geo_expression_file, quote = "", data.table = FALSE, check.names = FALSE)
required_cols <- c("gene_id", "gene_name", samples$count_column)
missing_cols <- setdiff(required_cols, names(raw_anno))
if (length(missing_cols) > 0) {
  stop("Missing expected columns in GEO matrix: ", paste(missing_cols, collapse = ", "))
}

raw_counts <- raw_anno[, c("gene_id", "gene_name", samples$count_column), drop = FALSE]
names(raw_counts) <- c("gene_id", "gene_name", samples$sample)
for (sample_id in samples$sample) {
  raw_counts[[sample_id]] <- as.integer(round(as.numeric(raw_counts[[sample_id]])))
}
fwrite(raw_counts, counts_file)
fwrite(samples[, c("sample", "group")], metadata_file)
fwrite(samples, metadata_detailed_file)

custom_marker_term2gene <- function() {
  marker_sets <- list(
    "Tuft cell markers" = c(
      "TRPM5", "ETV1", "POU2F3", "SOX4", "DKK3", "TLR4", "TCF4", "WNT5A",
      "GFI1B", "DCLK1", "AVIL", "LRMP", "SH2D6", "IL17RB", "KIT", "CHAT",
      "IL25", "ALOX5AP"
    ),
    "gamma delta T cell markers" = c(
      "TRDC", "TRGC1", "TRGC2", "TRDV1", "TRDV2", "TRGV9",
      "CD3D", "CD3E", "CD3G", "CD247", "CD2", "CD5", "CD7",
      "CCL5", "GNLY", "NKG7", "SLC4A10", "KLRB1", "CCR6",
      "GZMA", "GZMB", "PRF1", "KLRD1", "KLRC1", "PTPRC"
    )
  )
  rbindlist(lapply(names(marker_sets), function(term) {
    data.frame(term = term, gene = unique(toupper(marker_sets[[term]])))
  }))
}

safe_file_part <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

pretty_pathway_label <- function(x, wrap_width = NULL) {
  out <- gsub("_", " ", x)
  out[x == "GSE10240_CTRL_VS_IL17_AND_IL22_STIM_PRIMARY_BRONCHIAL_EPITHELIAL_CELLS_UP"] <-
    "IL-17/IL-22-stimulated epithelial response"
  out <- gsub("^REACTOME ", "Reactome: ", out)
  out <- gsub("^KEGG ", "KEGG: ", out)
  out <- gsub("\\bTCR\\b", "TCR", out)
  out <- gsub("\\bINTERLEUKIN 17\\b", "IL-17", out)
  out <- tools::toTitleCase(tolower(out))
  out <- gsub("\\bTcr\\b", "TCR", out)
  out <- gsub("\\bIl-17\\b", "IL-17", out)
  out <- gsub("\\bKegg\\b", "KEGG", out)
  out <- gsub("\\bReactome\\b", "Reactome", out)
  if (!is.null(wrap_width)) {
    out <- vapply(out, function(label) {
      paste(strwrap(label, width = wrap_width), collapse = "\n")
    }, FUN.VALUE = character(1))
  }
  out
}

collection_name <- function(pathway) {
  out <- sub("^([^_]+).*", "\\1", pathway)
  out[pathway %in% c("Tuft cell markers", "gamma delta T cell markers")] <- "CUSTOM"
  out
}

message("Step 3: edgeR CPM and DESeq2 differential expression")
raw <- fread(counts_file, data.table = FALSE, check.names = FALSE)
meta <- fread(metadata_file, data.table = FALSE)
stopifnot(all(meta$sample %in% names(raw)))

counts <- as.matrix(raw[, meta$sample, drop = FALSE])
storage.mode(counts) <- "integer"
rownames(counts) <- raw$gene_id
gene_map <- data.frame(
  gene_id = raw$gene_id,
  gene_name = toupper(raw$gene_name),
  stringsAsFactors = FALSE
)

cpm_file <- file.path(processed_dir, "GSE202479_edgeR_cpm_from_rawcounts.csv")
log2cpm_file <- file.path(processed_dir, "GSE202479_edgeR_log2cpm_from_rawcounts.csv")
deseq_result_file <- file.path(deseq_dir, "GSE202479_DESeq2_gallstones_vs_normal_from_rawcounts.csv")
normalized_counts_file <- file.path(deseq_dir, "GSE202479_DESeq2_normalized_counts.csv")
rank_file <- file.path(deseq_dir, "GSE202479_DESeq2_stat_rank_for_GSEA.csv")

cpm_mat <- edgeR::cpm(counts, log = FALSE)
log2cpm_mat <- edgeR::cpm(counts, log = TRUE, prior.count = 1)
fwrite(data.frame(gene_map, cpm_mat, check.names = FALSE), cpm_file)
fwrite(data.frame(gene_map, log2cpm_mat, check.names = FALSE), log2cpm_file)

if (file.exists(deseq_result_file) && file.exists(rank_file) && !force_deseq) {
  message("Using cached DESeq2 outputs")
  res_df <- fread(deseq_result_file, data.table = FALSE)
  rank_df <- fread(rank_file, data.table = FALSE)
} else {
  keep <- rowSums(cpm_mat >= 1) >= min(table(meta$group))
  counts_filtered <- counts[keep, , drop = FALSE]
  gene_map_filtered <- gene_map[keep, , drop = FALSE]
  message(
    "Genes before filtering: ", nrow(counts),
    "; after CPM>=1 in >=3 samples: ", nrow(counts_filtered)
  )

  coldata <- data.frame(
    row.names = meta$sample,
    group = factor(meta$group, levels = c("Normal", "Gallstones"))
  )
  dds <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData = coldata,
    design = ~ group
  )
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", "Gallstones", "Normal"))
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- merge(gene_map_filtered, res_df, by = "gene_id", sort = FALSE)

  norm_counts <- counts(dds, normalized = TRUE)
  fwrite(data.frame(gene_map_filtered, norm_counts, check.names = FALSE), normalized_counts_file)

  mean_cpm_normal <- rowMeans(cpm_mat[res_df$gene_id, meta$sample[meta$group == "Normal"], drop = FALSE])
  mean_cpm_gall <- rowMeans(cpm_mat[res_df$gene_id, meta$sample[meta$group == "Gallstones"], drop = FALSE])
  res_df$mean_CPM_normal <- mean_cpm_normal
  res_df$mean_CPM_gallstones <- mean_cpm_gall
  res_df <- res_df[, c(
    "gene_id", "gene_name", "baseMean", "log2FoldChange", "lfcSE",
    "stat", "pvalue", "padj", "mean_CPM_normal", "mean_CPM_gallstones"
  )]
  res_df <- res_df[order(res_df$padj, -abs(res_df$stat)), ]
  fwrite(res_df, deseq_result_file)

  rank_df <- res_df[
    is.finite(res_df$stat) & !is.na(res_df$gene_name) & res_df$gene_name != "",
  ]
  rank_df <- rank_df[order(rank_df$gene_name, -abs(rank_df$stat), -rank_df$baseMean), ]
  rank_df <- rank_df[!duplicated(rank_df$gene_name), ]
  rank_df <- rank_df[order(-rank_df$stat), ]
  fwrite(
    rank_df[, c("gene_id", "gene_name", "stat", "log2FoldChange", "pvalue", "padj", "baseMean")],
    rank_file
  )
}

message("Step 4: clusterProfiler GSEA using full MSigDB")
rank_df <- fread(rank_file, data.table = FALSE)
rank_df <- rank_df[is.finite(rank_df$stat) & !is.na(rank_df$gene_name) & rank_df$gene_name != "", ]
rank_df$gene_name <- toupper(rank_df$gene_name)
rank_df <- rank_df[order(rank_df$gene_name, -abs(rank_df$stat), -rank_df$baseMean), ]
rank_df <- rank_df[!duplicated(rank_df$gene_name), ]
gene_list <- rank_df$stat
names(gene_list) <- rank_df$gene_name
gene_list <- sort(gene_list, decreasing = TRUE)

term2gene <- clusterProfiler::read.gmt(msigdb_file)
term2gene$term <- as.character(term2gene$term)
term2gene$gene <- toupper(as.character(term2gene$gene))
term2gene <- rbind(term2gene, custom_marker_term2gene())
term2name <- unique(term2gene[, c("term", "term")])
colnames(term2name) <- c("term", "name")

gsea_rds <- file.path(gsea_dir, "GSE202479_clusterProfiler_GSEA_gseaResult.rds")
if (file.exists(gsea_rds) && !force_gsea && !force_deseq) {
  message("Using cached clusterProfiler gseaResult: ", gsea_rds)
  gsea_res <- readRDS(gsea_rds)
} else {
  set.seed(20260530)
  gsea_res <- clusterProfiler::GSEA(
    geneList = gene_list,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 0,
    seed = TRUE,
    by = "fgsea",
    verbose = FALSE
  )
  saveRDS(gsea_res, gsea_rds)
}

gsea_df <- as.data.frame(gsea_res)
gsea_df$direction <- ifelse(gsea_df$NES >= 0, "enriched_in_gallstones", "enriched_in_normal")
gsea_df$collection <- collection_name(gsea_df$ID)
gsea_df <- gsea_df[order(gsea_df$p.adjust, -abs(gsea_df$NES)), ]
gsea_res@result$Description <- pretty_pathway_label(gsea_res@result$ID)

gsea_full_csv <- file.path(gsea_dir, "GSE202479_clusterProfiler_GSEA_full_msigdb_DESeq2stat.csv")
fwrite(gsea_df, gsea_full_csv)

target_pathways <- c(
  "GSE10240_CTRL_VS_IL17_AND_IL22_STIM_PRIMARY_BRONCHIAL_EPITHELIAL_CELLS_UP",
  "REACTOME_INTERLEUKIN_17_SIGNALING",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "REACTOME_TCR_SIGNALING",
  "REACTOME_DOWNSTREAM_TCR_SIGNALING",
  "Tuft cell markers",
  "gamma delta T cell markers"
)
target_df <- gsea_df[gsea_df$ID %in% target_pathways, ]
target_csv <- file.path(gsea_dir, "GSE202479_clusterProfiler_GSEA_target_pathways_IL17_TCR_custom.csv")
fwrite(target_df, target_csv)

il17_pattern <- "IL17|IL_17|INTERLEUKIN_17|TH17|T_HELPER_17|RORC|IL23|IL_23"
il17_keep <- grepl(il17_pattern, gsea_df$ID, ignore.case = TRUE) |
  grepl(il17_pattern, gsea_df$Description, ignore.case = TRUE)
il17_df <- gsea_df[il17_keep, ]
il17_df <- il17_df[order(il17_df$p.adjust, -il17_df$NES), ]
il17_csv <- file.path(gsea_dir, "GSE202479_clusterProfiler_GSEA_IL17_TH17_related_candidates.csv")
fwrite(il17_df, il17_csv)

message("Step 5: enrichplot GSEA plots")
plot_one <- function(pathway_id) {
  if (!pathway_id %in% gsea_res@result$ID) {
    warning("Pathway not available after GSEA filtering: ", pathway_id)
    return(invisible(NULL))
  }
  row <- gsea_df[gsea_df$ID == pathway_id, , drop = FALSE]
  direction_label <- ifelse(row$NES[1] >= 0, "Gallstones enriched", "Normal enriched")
  title <- sprintf(
    "%s\nNES = %.2f, FDR = %.3g, %s",
    pretty_pathway_label(pathway_id), row$NES[1], row$p.adjust[1], direction_label
  )
  line_color <- ifelse(row$NES[1] >= 0, "#D55E00", "#0072B2")
  p <- enrichplot::gseaplot2(
    gsea_res,
    geneSetID = pathway_id,
    title = title,
    pvalue_table = TRUE,
    base_size = 11,
    color = line_color,
    ES_geom = "line"
  )
  png_file <- file.path(gsea_dir, paste0("GSE202479_enrichplot_gseaplot2_", safe_file_part(pathway_id), ".png"))
  pdf_file <- file.path(gsea_dir, paste0("GSE202479_enrichplot_gseaplot2_", safe_file_part(pathway_id), ".pdf"))
  ggsave(png_file, p, width = 10.5, height = 6.5, dpi = 300, limitsize = FALSE)
  ggsave(pdf_file, p, width = 10.5, height = 6.5, limitsize = FALSE)
  invisible(p)
}
invisible(lapply(target_pathways, plot_one))

dot_plot <- suppressMessages(
  enrichplot::dotplot(gsea_res, showCategory = 10, split = ".sign") +
    facet_grid(. ~ .sign) +
    scale_y_discrete(labels = function(x) pretty_pathway_label(x, wrap_width = 36)) +
    scale_color_gradient(low = "#2C7BB6", high = "#D7191C", name = "p.adjust") +
    labs(title = "clusterProfiler GSEA: full MSigDB v2026.1", x = "Gene ratio", y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 8, lineheight = 0.9)
    )
)
ggsave(
  file.path(gsea_dir, "GSE202479_enrichplot_dotplot_top10_by_sign.png"),
  dot_plot, width = 12, height = 7, dpi = 300, limitsize = FALSE
)
ggsave(
  file.path(gsea_dir, "GSE202479_enrichplot_dotplot_top10_by_sign.pdf"),
  dot_plot, width = 12, height = 7, limitsize = FALSE
)

summary_df <- data.frame(
  item = c(
    "project_dir", "geo_expression_url", "geo_expression_file", "geo_series_url",
    "geo_series_file", "msigdb_url", "msigdb_file", "counts_file",
    "metadata_file", "comparison", "cpm_method", "deseq2_input",
    "filter", "rank_metric", "gene_sets_plus_custom", "tested_gene_sets",
    "gsea_package", "plot_package", "gsea_full_csv", "target_csv", "il17_csv"
  ),
  value = c(
    project_dir, geo_expression_url, geo_expression_file, geo_series_url,
    geo_series_file, msigdb_url, msigdb_file, counts_file,
    metadata_file, "Gallstones vs Normal", "edgeR::cpm(raw counts)",
    "integer raw counts after CPM-based low-expression filtering",
    "CPM >= 1 in at least 3 samples",
    "DESeq2 Wald statistic", length(unique(term2gene$term)), nrow(gsea_df),
    paste0("clusterProfiler ", as.character(packageVersion("clusterProfiler"))),
    paste0("enrichplot ", as.character(packageVersion("enrichplot"))),
    gsea_full_csv, target_csv, il17_csv
  )
)
summary_file <- file.path(gsea_dir, "GSE202479_full_pipeline_summary.csv")
fwrite(summary_df, summary_file)

message("Done.")
message("Key target pathway results:")
print(target_df[, c("ID", "NES", "p.adjust", "direction", "setSize")], row.names = FALSE)
message("Outputs:")
message("  DESeq2: ", deseq_result_file)
message("  Full GSEA: ", gsea_full_csv)
message("  Target GSEA: ", target_csv)
message("  IL17/TH17 candidates: ", il17_csv)
message("  Summary: ", summary_file)
