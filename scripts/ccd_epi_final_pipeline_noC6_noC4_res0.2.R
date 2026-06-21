#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

set.seed(1234)

input_rds <- "gall0603.rds"
out_dir <- "epi"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

resolution_initial <- 0.5
resolution_final <- 0.2
dims_use <- 1:30
npcs_use <- 50

selected_marker_genes <- list(
  Epithelial = c("Epcam", "Krt8", "Krt18", "Krt19"),
  Cycling = c("Mki67", "Top2a", "Birc5", "Pclaf"),
  Selected = c(
    "Aqp1", "Slc2a2", "Cyp2f2", "Cftr", "Nlrp6",
    "Reg1", "S100g", "Slc4a2", "Klf5",
    "Cldn4", "Mt1", "Mt2", "Tff3",
    "C3", "Lcn2", "Cxcl5"
  )
)

paper_dotplot_genes <- c(
  "Sox4", "Hes1", "Hnf1b", "Onecut3", "Onecut1", "Hhex",
  "Pclaf", "Ccna2", "Cenpf", "Birc5", "Mki67",
  "Ctrb1", "Plet1", "Cxcl5", "C3", "Lcn2", "Chil4",
  "Ets2", "Fosl1", "Mt2", "Cldn4", "Mt1", "Tff3",
  "Klf5", "S100g", "Reg1", "Slc4a2", "Cftr", "Cyp2f2",
  "Nlrp6", "Lamb1", "Slc2a2", "Aqp1",
  "Epcam", "Krt19"
)

celltype_order <- c("Aqp1+", "Reg1+", "Cldn4+", "C3+", "Mki67+")

run_standard_workflow <- function(obj, resolution) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, features = rownames(obj))
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs_use, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims_use)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj, dims = dims_use)
  obj
}

new_object_from_counts <- function(source_obj, keep_cells, project) {
  DefaultAssay(source_obj) <- "RNA"
  counts <- GetAssayData(source_obj, assay = "RNA", layer = "counts")[, keep_cells, drop = FALSE]
  drop_meta <- c(
    "seurat_clusters",
    grep("_snn_res\\.", colnames(source_obj@meta.data), value = TRUE)
  )
  meta <- source_obj@meta.data[
    keep_cells,
    setdiff(colnames(source_obj@meta.data), drop_meta),
    drop = FALSE
  ]
  obj <- CreateSeuratObject(counts = counts, meta.data = meta, project = project)
  obj$orig.ident <- factor(obj$orig.ident, levels = c("WT_CCD", "KO_CCD"))
  obj
}

write_marker_tables <- function(obj, markers_csv, top100_csv) {
  DefaultAssay(obj) <- "RNA"
  markers <- FindAllMarkers(
    obj,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  markers <- markers[order(markers$cluster, markers$p_val_adj, markers$p_val, -markers$avg_log2FC), ]
  write.csv(markers, markers_csv, row.names = FALSE)

  top100 <- do.call(
    rbind,
    lapply(split(markers, markers$cluster), function(x) {
      x <- x[order(-x$avg_log2FC, x$p_val_adj, x$p_val), , drop = FALSE]
      head(x, 100)
    })
  )
  rownames(top100) <- NULL
  write.csv(top100, top100_csv, row.names = FALSE)
}

write_cluster_proportions <- function(obj, group_col, prefix, levels_use = NULL) {
  groups_raw <- as.character(obj[[group_col, drop = TRUE]])
  if (is.null(levels_use)) {
    levels_use <- sort(unique(as.numeric(groups_raw)))
  }
  groups <- factor(groups_raw, levels = levels_use)
  samples <- factor(as.character(obj$orig.ident), levels = c("WT_CCD", "KO_CCD"))
  tab <- table(orig.ident = samples, group = groups)

  counts <- as.data.frame.matrix(tab)
  counts$orig.ident <- rownames(counts)
  counts$total <- rowSums(tab)
  counts <- counts[, c("orig.ident", setdiff(colnames(counts), c("orig.ident", "total")), "total")]
  write.csv(counts, paste0(prefix, "_counts.csv"), row.names = FALSE)

  prop_long <- as.data.frame(tab)
  colnames(prop_long) <- c("orig.ident", "group", "count")
  totals <- aggregate(count ~ orig.ident, prop_long, sum)
  colnames(totals)[2] <- "orig_total"
  prop_long <- merge(prop_long, totals, by = "orig.ident")
  prop_long$proportion <- prop_long$count / prop_long$orig_total
  prop_long$percent <- 100 * prop_long$proportion
  write.csv(prop_long, paste0(prefix, "_proportion_by_orig.csv"), row.names = FALSE)

  percent <- round(100 * as.data.frame.matrix(prop.table(tab, margin = 1)), 2)
  write.csv(percent, paste0(prefix, "_percent_by_orig.csv"))

  p <- ggplot(prop_long, aes(x = orig.ident, y = proportion, fill = group)) +
    geom_col(width = 0.65, color = "white", linewidth = 0.2) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Proportion within orig.ident", fill = group_col) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  ggsave(paste0(prefix, ".pdf"), p, width = 5.5, height = 4.5, useDingbats = FALSE)
  ggsave(paste0(prefix, ".png"), p, width = 5.5, height = 4.5, dpi = 300)
}

plot_basic_umap_and_dotplot <- function(obj, prefix, group_col = "seurat_clusters") {
  p_umap <- DimPlot(obj, reduction = "umap", group.by = group_col, label = TRUE) + NoLegend()
  ggsave(paste0(prefix, "_umap.pdf"), p_umap, width = 7, height = 5.5, useDingbats = FALSE)
  ggsave(paste0(prefix, "_umap.png"), p_umap, width = 7, height = 5.5, dpi = 300)

  missing_features <- setdiff(unlist(selected_marker_genes, use.names = FALSE), rownames(obj))
  if (length(missing_features) > 0) {
    stop("Missing genes in object: ", paste(missing_features, collapse = ", "))
  }

  p_dot <- DotPlot(
    obj,
    features = selected_marker_genes,
    group.by = group_col,
    cols = c("grey90", "#D73027"),
    dot.scale = 6
  ) +
    RotatedAxis() +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title = element_blank()
    )
  ggsave(paste0(prefix, "_dotplot.pdf"), p_dot, width = 12, height = 5.5, useDingbats = FALSE)
  ggsave(paste0(prefix, "_dotplot.png"), p_dot, width = 12, height = 5.5, dpi = 300)
}

plot_final_paper_dotplot <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj$plot_celltype <- factor(as.character(obj$epi_annotation), levels = celltype_order)

  missing_genes <- setdiff(paper_dotplot_genes, rownames(obj))
  if (length(missing_genes) > 0) {
    stop("Missing genes in final object: ", paste(missing_genes, collapse = ", "))
  }

  p0 <- DotPlot(
    obj,
    features = paper_dotplot_genes,
    group.by = "plot_celltype",
    cols = c("#fff5f0", "#b21818"),
    dot.scale = 8
  )

  dot_df <- p0$data
  dot_df$id <- factor(dot_df$id, levels = celltype_order)
  dot_df$features.plot <- factor(dot_df$features.plot, levels = rev(paper_dotplot_genes))

  p_dot <- ggplot(dot_df, aes(x = id, y = features.plot, size = pct.exp, color = avg.exp.scaled)) +
    geom_point() +
    scale_color_gradient2(
      low = "#fff5f0",
      mid = "#fcae91",
      high = "#b21818",
      midpoint = 0,
      limits = c(-1.5, 1.5),
      oob = squish,
      name = "Average Expression"
    ) +
    scale_size(
      range = c(0, 8),
      breaks = c(25, 50, 75, 100),
      limits = c(0, 100),
      name = "Percent Expressed"
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, face = "bold.italic", color = "black", size = 20),
      axis.text.y = element_text(face = "bold.italic", color = "black", size = 20),
      legend.title = element_text(face = "bold", size = 18),
      legend.text = element_text(face = "bold", size = 15),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8)
    )

  ggsave(file.path(out_dir, "dotplot_final_annotation_ordered_with_Epcam_Krt19.pdf"), p_dot, width = 6.8, height = 10.5, useDingbats = FALSE)
  ggsave(file.path(out_dir, "dotplot_final_annotation_ordered_with_Epcam_Krt19.png"), p_dot, width = 6.8, height = 10.5, dpi = 300)
}

plot_final_barplot <- function(obj) {
  obj$plot_celltype <- factor(as.character(obj$epi_annotation), levels = celltype_order)
  obj$plot_group <- recode(
    as.character(obj$orig.ident),
    "WT_CCD" = "WT",
    "KO_CCD" = "Pou2f3-/-"
  )
  obj$plot_group <- factor(obj$plot_group, levels = c("WT", "Pou2f3-/-"))

  bar_df <- obj@meta.data %>%
    count(plot_group, plot_celltype, name = "n") %>%
    group_by(plot_group) %>%
    mutate(
      proportion = n / sum(n),
      label = paste0(round(proportion * 100, 1), "%")
    ) %>%
    ungroup()

  sig_df <- data.frame(
    plot_celltype = factor(celltype_order, levels = celltype_order),
    significance = c("ns", "ns", "*", "*", "*")
  ) %>%
    left_join(
      bar_df %>%
        group_by(plot_celltype) %>%
        summarize(y = max(proportion) + 0.055, .groups = "drop"),
      by = "plot_celltype"
    ) %>%
    mutate(
      x = as.numeric(plot_celltype),
      x_start = x - 0.22,
      x_end = x + 0.22,
      y0 = y - 0.015,
      y_text = y + 0.01
    )

  p_bar <- ggplot(bar_df, aes(x = plot_celltype, y = proportion, fill = plot_group)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65, color = NA) +
    geom_text(
      aes(label = label),
      position = position_dodge(width = 0.7),
      vjust = -0.25,
      size = 5,
      fontface = "bold"
    ) +
    geom_segment(data = sig_df, aes(x = x_start, xend = x_end, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.5, color = "black") +
    geom_segment(data = sig_df, aes(x = x_start, xend = x_start, y = y, yend = y0), inherit.aes = FALSE, linewidth = 0.5, color = "black") +
    geom_segment(data = sig_df, aes(x = x_end, xend = x_end, y = y, yend = y0), inherit.aes = FALSE, linewidth = 0.5, color = "black") +
    geom_text(data = sig_df, aes(x = x, y = y_text, label = significance), inherit.aes = FALSE, size = 5, fontface = "bold") +
    scale_fill_manual(values = c("WT" = "#9EC6DD", "Pou2f3-/-" = "#4F73B6")) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.50), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Proportion", fill = NULL) +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x = element_text(face = "bold.italic", color = "black", size = 18),
      axis.text.y = element_text(face = "bold", color = "black", size = 16),
      axis.title.y = element_text(face = "bold", color = "black", size = 20),
      legend.position = c(0.75, 0.88),
      legend.direction = "horizontal",
      legend.text = element_text(face = "bold", size = 16),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8)
    )

  ggsave(file.path(out_dir, "barplot_celltype_proportion_annotated_sig.pdf"), p_bar, width = 8.6, height = 4.5, useDingbats = FALSE)
  ggsave(file.path(out_dir, "barplot_celltype_proportion_annotated_sig.png"), p_bar, width = 8.6, height = 4.5, dpi = 300)
  write.csv(bar_df, file.path(out_dir, "barplot_celltype_proportion_annotated_sig_data.csv"), row.names = FALSE)
  write.csv(sig_df, file.path(out_dir, "barplot_celltype_proportion_annotated_sig_labels.csv"), row.names = FALSE)
}

message("1. Extract WT_CCD/KO_CCD epithelial cells")
original <- readRDS(input_rds)
DefaultAssay(original) <- "RNA"
required_meta <- c("celltype2", "orig.ident")
missing_meta <- setdiff(required_meta, colnames(original@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing metadata column(s): ", paste(missing_meta, collapse = ", "))
}

ccd_epi_cells <- colnames(original)[
  original$celltype2 == "Epithelial" &
    original$orig.ident %in% c("WT_CCD", "KO_CCD")
]
if (length(ccd_epi_cells) == 0) {
  stop("No cells matched celltype2 == 'Epithelial' and orig.ident in WT_CCD/KO_CCD.")
}

ccdepi <- new_object_from_counts(original, ccd_epi_cells, "CCD_Epithelial")
saveRDS(ccdepi, file.path(out_dir, "ccdepi.rds"))

message("2. Initial clustering at resolution 0.5 and remove contaminating cluster 10")
initial <- run_standard_workflow(ccdepi, resolution_initial)
saveRDS(initial, file.path(out_dir, "ccdepi_initial_res0.5.rds"))
cluster10_cells <- colnames(initial)[as.character(Idents(initial)) == "10"]
if (length(cluster10_cells) == 0) {
  stop("Cluster 10 was not found in initial resolution 0.5 clustering.")
}
writeLines(cluster10_cells, file.path(out_dir, "removed_cluster10_cells.txt"))

no_cluster10_cells <- setdiff(colnames(ccdepi), cluster10_cells)
ccdepi_no_cluster10 <- new_object_from_counts(ccdepi, no_cluster10_cells, "CCD_Epithelial_no_cluster10")
saveRDS(ccdepi_no_cluster10, file.path(out_dir, "ccdepi_without_cluster10.rds"))

message("3. Cluster no-cluster10 object at resolution 0.2 and remove C6")
res0.2 <- run_standard_workflow(ccdepi_no_cluster10, resolution_final)
saveRDS(res0.2, file.path(out_dir, "ccdepi_processed_res0.2.rds"))
write_marker_tables(res0.2, file.path(out_dir, "markers_res0.2.csv"), file.path(out_dir, "top100markers_res0.2.csv"))
plot_basic_umap_and_dotplot(res0.2, file.path(out_dir, "res0.2"), group_col = "seurat_clusters")
write_cluster_proportions(res0.2, "seurat_clusters", file.path(out_dir, "cluster_res0.2"))

c6_cells <- colnames(res0.2)[as.character(Idents(res0.2)) == "6"]
if (length(c6_cells) == 0) {
  stop("C6 was not found in the resolution 0.2 result.")
}
writeLines(c6_cells, file.path(out_dir, "removed_C6_res0.2_cells.txt"))

no_c6_cells <- setdiff(colnames(ccdepi_no_cluster10), c6_cells)
ccdepi_no_c6 <- new_object_from_counts(ccdepi_no_cluster10, no_c6_cells, "CCD_Epithelial_no_C6_res0.2")
saveRDS(ccdepi_no_c6, file.path(out_dir, "ccdepi_without_C6_res0.2.rds"))

message("4. Recluster after removing C6, then remove C4")
no_c6_clustered <- run_standard_workflow(ccdepi_no_c6, resolution_final)
saveRDS(no_c6_clustered, file.path(out_dir, "ccdepi_processed_noC6_res0.2.rds"))
write_marker_tables(no_c6_clustered, file.path(out_dir, "markers_noC6_res0.2.csv"), file.path(out_dir, "top100markers_noC6_res0.2.csv"))
plot_basic_umap_and_dotplot(no_c6_clustered, file.path(out_dir, "noC6_res0.2"), group_col = "seurat_clusters")
write_cluster_proportions(no_c6_clustered, "seurat_clusters", file.path(out_dir, "cluster_noC6_res0.2"))

c4_after_no_c6_cells <- colnames(no_c6_clustered)[as.character(Idents(no_c6_clustered)) == "4"]
if (length(c4_after_no_c6_cells) == 0) {
  stop("C4 was not found after removing C6 and reclustering at resolution 0.2.")
}
writeLines(c4_after_no_c6_cells, file.path(out_dir, "removed_C4_after_noC6_res0.2_cells.txt"))

no_c6_no_c4_cells <- setdiff(colnames(ccdepi_no_c6), c4_after_no_c6_cells)
ccdepi_no_c6_no_c4 <- new_object_from_counts(
  ccdepi_no_c6,
  no_c6_no_c4_cells,
  "CCD_Epithelial_no_C6_no_C4_res0.2"
)
saveRDS(ccdepi_no_c6_no_c4, file.path(out_dir, "ccdepi_without_C4_after_noC6_res0.2.rds"))

message("5. Final clustering after removing C6 and C4 at resolution 0.2")
final <- run_standard_workflow(ccdepi_no_c6_no_c4, resolution_final)
saveRDS(final, file.path(out_dir, "ccdepi_processed_noC6_noC4_res0.2.rds"))
write_marker_tables(
  final,
  file.path(out_dir, "markers_noC6_noC4_res0.2.csv"),
  file.path(out_dir, "top100markers_noC6_noC4_res0.2.csv")
)
plot_basic_umap_and_dotplot(final, file.path(out_dir, "noC6_noC4_res0.2"), group_col = "seurat_clusters")
write_cluster_proportions(final, "seurat_clusters", file.path(out_dir, "cluster_noC6_noC4_res0.2"))

message("6. Add final epithelial annotations")
cluster_to_annotation <- c(
  "0" = "Reg1+",
  "1" = "Aqp1+",
  "2" = "Cldn4+",
  "3" = "Mki67+",
  "4" = "C3+"
)

clusters <- as.character(Idents(final))
missing_clusters <- setdiff(unique(clusters), names(cluster_to_annotation))
if (length(missing_clusters) > 0) {
  stop("Missing annotation for final cluster(s): ", paste(missing_clusters, collapse = ", "))
}

final$seurat_clusters_noC6_noC4_res0.2 <- clusters
final$epi_annotation <- factor(unname(cluster_to_annotation[clusters]), levels = celltype_order)
Idents(final) <- "epi_annotation"

write.csv(
  data.frame(
    cluster = names(cluster_to_annotation),
    epi_annotation = unname(cluster_to_annotation),
    row.names = NULL
  ),
  file.path(out_dir, "cluster_annotation_mapping_noC6_noC4_res0.2.csv"),
  row.names = FALSE
)

saveRDS(final, file.path(out_dir, "ccdepi_processed_noC6_noC4_res0.2_annotated.rds"))
plot_basic_umap_and_dotplot(final, file.path(out_dir, "epi_annotation_noC6_noC4_res0.2"), group_col = "epi_annotation")
write_cluster_proportions(
  final,
  "epi_annotation",
  file.path(out_dir, "epi_annotation_noC6_noC4_res0.2"),
  levels_use = celltype_order
)

message("7. Draw final publication-style dotplot and proportion barplot")
plot_final_paper_dotplot(final)
plot_final_barplot(final)

message("Final annotated object: ", file.path(out_dir, "ccdepi_processed_noC6_noC4_res0.2_annotated.rds"))
message("Done")
