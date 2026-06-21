#!/usr/bin/env Rscript
# Author: LYU WEI
# Date: 2026-05-30
# Version: 1.1
# scRNA-seq cleaning pipeline (mouse, WT vs OLD):
#   QC filtering + doublet removal (scDblFinder) + CCA integration + clustering.
#
# Install once:
#   install.packages(c("Seurat","dplyr","ggplot2","patchwork","Matrix"))
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("scDblFinder","SingleCellExperiment"))

suppressPackageStartupMessages({
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)
setwd("~/Desktop/old")

# ---- 1. parameters ----
samples <- data.frame(
  sample_id = c("WT", "OLD"),     # becomes orig.ident
  dir       = c("WT", "24mo"),    # Read10X folder (filtered matrix)
  stringsAsFactors = FALSE
)

# QC thresholds — inspect the *_prefilter.pdf plots first, then tune
min_features <- 200      # lower bound on detected genes
max_features <- 7500     # upper bound (very high = likely multiplet)
max_counts   <- 60000    # upper bound on UMIs
max_mt       <- 15       # max % mitochondrial
max_hb       <- 5        # max % hemoglobin (RBC contamination)
min_cplx     <- 0.80     # min log10(genes)/log10(UMI)  (low = low-complexity/dying)

# ---- 2. per-sample: Seurat -> QC -> doublets ----
process_sample <- function(dir, sid) {
  message(">> processing ", sid)
  
  raw <- Read10X(data.dir = dir)
  obj <- CreateSeuratObject(counts = raw, project = sid,
                            min.cells = 3, min.features = min_features)
  
  # QC metrics (mouse symbols are lowercase-initial)
  obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")
  obj[["percent.hb"]]   <- PercentageFeatureSet(obj, pattern = "^Hb[ab]")
  obj$log10GenesPerUMI  <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  
  # pre-filter QC plot for threshold tuning
  ggsave(paste0("qc_", sid, "_prefilter.pdf"),
         VlnPlot(obj, c("nFeature_RNA","nCount_RNA","percent.mt","percent.hb"),
                 pt.size = 0, ncol = 4),
         width = 12, height = 4)
  
  # QC filtering
  obj <- subset(obj, subset =
                  nFeature_RNA >= min_features & nFeature_RNA <= max_features &
                  nCount_RNA   <= max_counts &
                  percent.mt   <= max_mt &
                  percent.hb   <= max_hb &
                  log10GenesPerUMI >= min_cplx)
  
  # doublet removal (scDblFinder on counts)
  set.seed(42)
  sce2 <- SingleCellExperiment(assays = list(counts = GetAssayData(obj, layer = "counts")))
  sce2 <- scDblFinder(sce2)
  obj$scDblFinder.class <- sce2$scDblFinder.class
  obj$scDblFinder.score <- sce2$scDblFinder.score
  obj <- subset(obj, subset = scDblFinder.class == "singlet")
  
  message("   kept ", ncol(obj), " cells")
  obj
}

obj.list <- mapply(process_sample, samples$dir, samples$sample_id, SIMPLIFY = FALSE)
names(obj.list) <- samples$sample_id

# ---- 3. CCA integration ----
obj.list <- lapply(obj.list, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  x
})
features <- SelectIntegrationFeatures(obj.list, nfeatures = 3000)
anchors  <- FindIntegrationAnchors(obj.list, anchor.features = features,
                                   reduction = "cca", dims = 1:30)
old <- IntegrateData(anchorset = anchors, dims = 1:30)

# ---- 4. scale / PCA / cluster / UMAP ----
DefaultAssay(old) <- "integrated"
old <- ScaleData(old, vars.to.regress = "percent.mt", verbose = FALSE)
old <- RunPCA(old, npcs = 30, verbose = FALSE)
print(ElbowPlot(old, ndims = 30))
old <- FindNeighbors(old, dims = 1:20)
old <- FindClusters(old, resolution = 0.3)
old <- RunUMAP(old, dims = 1:20)

# ---- 5. verify cleanup ----
DefaultAssay(old) <- "RNA"
old <- JoinLayers(old)
print(VlnPlot(old, features = "Il17a", group.by = "seurat_clusters", pt.size = 0) + NoLegend())

# saveRDS(old, "old_clean.rds")





cols <- c("#DC143C","#20B2AA","#FFA500","#9370DB","#228B22","#1E90FF","#FA8072","#EE82EE","#7B68EE",
          "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887","#32CD32")


celltype=data.frame(ClusterID=0:19,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(1),2]='abT'
celltype[celltype$ClusterID %in% c(0,7,18),2]='Epithelial'
celltype[celltype$ClusterID %in% c(6,17),2]='Fibroblast'
celltype[celltype$ClusterID %in% c(8,11,12,15,16),2]='Mac/DC'
celltype[celltype$ClusterID %in% c(2,9),2]='Endothelial'
celltype[celltype$ClusterID %in% c(4,13),2]='B'
celltype[celltype$ClusterID %in% c(14),2]='Tuft'
celltype[celltype$ClusterID %in% c(5),2]='gdT'
celltype[celltype$ClusterID %in% c(3,10,19),2]='Neutrophil'



head(celltype)
celltype 
table(celltype$celltype)
old@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  old@meta.data[which(old@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(old@meta.data$celltype)










old$orig.ident<-factor(old$orig.ident,levels =c("WT", "OLD"))


# 加载必要的包
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)

# 从 Seurat 对象中提取数据并计算比例
cell_prop <- as.data.frame(old@meta.data) %>%
  group_by(orig.ident, celltype) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(orig.ident) %>%
  mutate(proportion = count/sum(count))

# 创建流图
# 创建优化后的流图
ggplot(cell_prop,
       aes(x = orig.ident, 
           stratum = celltype, 
           alluvium = celltype,
           y = proportion,
           fill = celltype)) +
  geom_stratum(width = 0.5) +  # 增加柱状图宽度
  geom_flow(width = 0.5) +     # 减少连接的长度
  scale_fill_manual(values = cols) +  # 使用自定义颜色
  scale_y_continuous(labels = scales::percent,
                     expand = c(0, 0)) +  # 从0开始的y轴
  scale_x_discrete(expand = c(0.2, 0.2)) +  # 调整x轴间距
  labs(x = "", 
       y = "Cell Type Proportion",
       title = "3ms vs 18ms") +
  theme_classic() +  # 使用经典主题，包含坐标轴
  theme(
    panel.grid = element_blank(),  # 移除网格线
    axis.line = element_line(colour = "black"),  # 添加坐标轴线
    axis.text = element_text(colour = "black"),  # 坐标轴文本颜色
    legend.title = element_text(face = "bold"),  # 图例标题加粗
    plot.title = element_text(hjust = 0.5)       # 标题居中
  )
ggsave("1-C.pdf", width = 4, height = 5)























suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

obj <- readRDS("~/Documents/scRNA/epi/ccdepi_processed_noC6_noC4_res0.2_annotated.rds")
DefaultAssay(obj) <- "RNA"

obj$plot_celltype <- recode(
  as.character(obj$epi_annotation),
  "aqp1+" = "Aqp1+",
  "Reg1+" = "Reg1+",
  "Cldn+" = "Cldn4+",
  "C3+" = "C3+",
  "Mki67+" = "Mki67+"
)

celltype_order <- c("Aqp1+", "Reg1+", "Cldn4+", "C3+", "Mki67+")
obj$plot_celltype <- factor(obj$plot_celltype, levels = celltype_order)

obj$plot_group <- recode(
  as.character(obj$orig.ident),
  "WT_CCD" = "WT",
  "KO_CCD" = "Pou2f3-/-"
)
obj$plot_group <- factor(obj$plot_group, levels = c("WT", "Pou2f3-/-"))

# bar plot: 每个 orig.ident 内部各 celltype 比例
bar_df <- obj@meta.data %>%
  count(plot_group, plot_celltype, name = "n") %>%
  group_by(plot_group) %>%
  mutate(
    proportion = n / sum(n),
    label = paste0(round(proportion * 100, 1), "%")
  ) %>%
  ungroup()

p_bar <- ggplot(bar_df, aes(x = plot_celltype, y = proportion, fill = plot_group)) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.65,
    color = NA
  ) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.7),
    vjust = -0.25,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c("WT" = "#9EC6DD", "Pou2f3-/-" = "#4F73B6")) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 0.55),
    expand = expansion(mult = c(0, 0.03))
  ) +
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

ggsave("./barplot_celltype_proportion_annotation.pdf", p_bar, width = 8.6, height = 4.5)





suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

DefaultAssay(obj) <- "RNA"

obj$plot_celltype <- dplyr::recode(
  as.character(obj$epi_annotation),
  "aqp1+" = "Aqp1+",
  "Reg1+" = "Reg1+",
  "Cldn+" = "Cldn4+",
  "C3+" = "C3+",
  "Mki67+" = "Mki67+"
)

celltype_order <- c("Aqp1+", "Reg1+", "Cldn4+", "C3+", "Mki67+")
obj$plot_celltype <- factor(obj$plot_celltype, levels = celltype_order)

dot_genes <- c(
  "Sox4", "Hes1", "Hnf1b", "Onecut3", "Onecut1", "Hhex",
  "Pclaf", "Ccna2", "Cenpf", "Birc5", "Mki67",
  "Ccl20", "Ccl28", "Ccl2", "Spp1", "Cxcl5", "C3", "Lcn2", "Chil4",
  "Ets2", "Fosl1", "Mt2", "Cldn4", "Mt1", "Tff3",
  "Klf5", "S100g", "Reg1", "Slc4a2", "Cftr", "Cyp2f2",
  "Nlrp6", "Lamb1", "Slc2a2", "Aqp1",
  "Epcam", "Krt19", "Krt7"
)

missing_genes <- setdiff(dot_genes, rownames(obj))
if (length(missing_genes) > 0) {
  stop("Missing genes: ", paste(missing_genes, collapse = ", "))
}

p0 <- DotPlot(
  obj,
  features = dot_genes,
  group.by = "plot_celltype",
  cols = c("#fff5f0", "#b21818"),
  dot.scale = 8
)

dot_df <- p0$data
dot_df$id <- factor(dot_df$id, levels = celltype_order)
dot_df$features.plot <- factor(dot_df$features.plot, levels = rev(dot_genes))

p_dot <- ggplot(
  dot_df,
  aes(x = id, y = features.plot, size = pct.exp, color = avg.exp.scaled)
) +
  geom_point() +
  scale_color_gradient2(
    low = "#fff5f0",
    mid = "#fcae91",
    high = "#b21818",
    midpoint = 0,
    limits = c(-1.5, 1.5),
    oob = scales::squish,
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
    axis.text.x = element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      face = "bold.italic",
      color = "black",
      size = 20
    ),
    axis.text.y = element_text(
      face = "bold.italic",
      color = "black",
      size = 20
    ),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(face = "bold", size = 15),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8)
  )

ggsave(
  "./dotplot_celltype_x_gene_y_with_Epcam_Krt19.png",
  p_dot,
  width = 7,
  height = 12,
  dpi = 300
)

ggsave(
  "./dotplot_celltype_x_gene_y_with_Epcam_Krt19.pdf",
  p_dot,
  width = 7,
  height = 12
)




#统计检验
library(scProportionTest)
prop_test <- sc_utils(obj)
prop_test <- permutation_test(
  prop_test, cluster_identity = "epi_annotation",
  sample_1 = "WT_CCD", sample_2 = "KO_CCD",
  sample_identity = "orig.ident"
)
permutation_plot(prop_test)






















UMAPPlot(obj, group.by = "epi_annotation", label = T, cols=cols, repel = T)+
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))





old$celltype<-factor(old$celltype,levels =c("gdT", "abT", "Mac/DC", "Neutrophil", "B", "Epithelial","Endothelial","Fibroblast", "Tuft"))


UMAPPlot(old, group.by = "celltype", label = T, cols=cols, repel = T)+
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))






#统计检验
library(scProportionTest)
prop_test <- sc_utils(old)
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "WT", sample_2 = "OLD",
  sample_identity = "orig.ident"
)
permutation_plot(prop_test)


pdf("newFig1M.pdf", width = 6, height = 10)
VlnPlot(old, "Il17a", group.by = "celltype", cols = cols, pt.size = 0.05)+  coord_flip()
dev.off()







library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

genes            <- c("Il17a", "Ifng")
cells_of_interest <- c("gdT", "abT")                 
display_name     <- c(gdT = "γδ T", abT = "αβ T")   

## ---- 1. 计算阳性比例 (gene × celltype × 样本) ----
sub <- subset(old, celltype %in% cells_of_interest)

pos_tbl <- FetchData(sub, vars = c("celltype", "orig.ident", genes)) %>%
  group_by(celltype, orig.ident) %>%
  summarise(across(all_of(genes), ~ mean(.x > 0) * 100), .groups = "drop")

## ---- 2. 准备 UMAP 绘图数据（所有细胞都画） ----
emb <- as.data.frame(Embeddings(old, "umap"))
colnames(emb)[1:2] <- c("UMAP1", "UMAP2")
plot_df <- cbind(emb, FetchData(old, vars = c("orig.ident", genes)))
plot_df$orig.ident <- factor(plot_df$orig.ident, levels = c("WT", "OLD"))

theme_umap <- theme_classic(base_size = 12) +
  theme(
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    axis.title       = element_text(size = 10, colour = "grey40"),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 13),
    plot.title       = element_text(face = "bold.italic", hjust = 0.5, size = 14),
    legend.key.height = unit(0.9, "cm"),
    legend.key.width  = unit(0.3, "cm")
  )

## ---- 3. 每个 panel 标注 γδT / αβT 的阳性比例 ----
make_panel <- function(gene) {
  d <- plot_df
  d$expr <- d[[gene]]
  d <- d[order(d$expr), ]                       # 阳性细胞画在上层，不被灰点盖住
  
  lab <- pos_tbl %>%
    mutate(celltype = factor(celltype, levels = cells_of_interest)) %>%
    arrange(orig.ident, celltype) %>%
    mutate(line = paste0(display_name[as.character(celltype)], ": ",
                         sprintf("%.1f", .data[[gene]]), "%")) %>%
    group_by(orig.ident) %>%
    summarise(label = paste(line, collapse = "\n"), .groups = "drop")
  lab$orig.ident <- factor(lab$orig.ident, levels = c("WT", "OLD"))
  lab$x <- min(d$UMAP1)
  lab$y <- max(d$UMAP2)
  
  ggplot(d, aes(UMAP1, UMAP2)) +
    geom_point(aes(colour = expr), size = 0.5, stroke = 0) +
    facet_wrap(~ orig.ident, nrow = 1) +
    scale_colour_gradientn(
      colours = c("grey88", "#fcbba1", "#fb6a4a", "#cb181d", "#67000d"),
      name = "Expr") +
    geom_text(data = lab, aes(x, y, label = label),
              hjust = 0, vjust = 1, size = 3.3, lineheight = 0.95,
              colour = "black", inherit.aes = FALSE) +
    labs(title = gene, x = "UMAP_1", y = "UMAP_2") +
    theme_umap
}

## ---- 4. 拼图：上 Il17a，下 Ifng ----
p <- make_panel("Il17a") / make_panel("Ifng")
p

ggsave("FeaturePlot_split_with_pct.pdf", p, width = 8, height = 7)

