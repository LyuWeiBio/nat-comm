#!/usr/bin/env bash
set -euo pipefail

# Minimal CellChat pipeline:
#   Seurat RDS -> per-condition CellChat -> merged CellChat ->
#   CellChat rankNet information-flow plots for TopN + condition-specific pathways.
#
# Default output is limited to:
#   results/cellchat/plots/information_flow_top5_plus_specific/
#   results/cellchat/tables/cellchat_information_flow_top5_plus_specific_*.csv
# plus the CellChat RDS objects needed for reproducibility.

usage() {
  cat <<'USAGE'
Usage:
  ./cellchat_top5_specific_flow_pipeline.sh [options]

Input and grouping:
  --input FILE              Seurat RDS input. Default: gall0603.rds
  --outdir DIR              CellChat output directory. Default: results/cellchat
  --conditions LIST         Comma-separated orig.ident values.
                            Default: WT_ND,WT_CCD,KO_ND,KO_CCD
  --comparisons LIST        Comma-separated compare_vs_base pairs.
                            Default: WT_CCD_vs_WT_ND,KO_ND_vs_WT_ND
  --condition-col COL       Metadata column for conditions. Default: orig.ident
  --group-col COL           Metadata column for cell types. Default: celltype2
  --celltype-rename MAP     Comma-separated old=new celltype rename map.
                            Default: CD4+ T=ab-T,STromal=Fibroblasts,Stromal=Fibroblasts
  --assay NAME              Seurat assay. Default: RNA
  --layer NAME              Seurat v5 layer or v4 slot. Default: data

CellChat:
  --workers N               Future workers. Default: 4
  --nboot N                 computeCommunProb nboot. Default: 100
  --min-cells N             Minimum cells per cell type. Default: 10
  --force                   Recompute per-condition CellChat objects.
  --skip-analysis           Reuse existing CellChat objects and only remake flow plots.

Information flow plot:
  --top-n N                 Top differential shared pathways to include. Default: 5
  --plot-width N            Plot width in inches. Default: 11
  --dpi N                   PNG DPI. Default: 300

Environment:
  ENV_PREFIX                Micromamba env prefix.
                            Default: /Users/lyuwei/micromamba/envs/bioinfo

Examples:
  ./cellchat_top5_specific_flow_pipeline.sh --skip-analysis
  ./cellchat_top5_specific_flow_pipeline.sh --force
  ./cellchat_top5_specific_flow_pipeline.sh --top-n 5 --skip-analysis
USAGE
}

need_value() {
  if [[ $# -lt 2 || "${2:-}" == --* ]]; then
    echo "Missing value for $1" >&2
    exit 2
  fi
}

INPUT="gall0603.rds"
OUTDIR="results/cellchat"
CONDITIONS="WT_ND,WT_CCD,KO_ND,KO_CCD"
COMPARISONS="WT_CCD_vs_WT_ND,KO_ND_vs_WT_ND"
CONDITION_COL="orig.ident"
GROUP_COL="celltype2"
CELLTYPE_RENAME="CD4+ T=ab-T,STromal=Fibroblasts,Stromal=Fibroblasts"
ASSAY="RNA"
LAYER="data"
WORKERS="4"
NBOOT="100"
MIN_CELLS="10"
TOP_N="5"
PLOT_WIDTH="11"
DPI="300"
FORCE=false
RUN_ANALYSIS=true
DRY_RUN=false
FAST_WILCOX="true"
DB_CATEGORY="all"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) need_value "$@"; INPUT="$2"; shift 2 ;;
    --outdir) need_value "$@"; OUTDIR="$2"; shift 2 ;;
    --conditions) need_value "$@"; CONDITIONS="$2"; shift 2 ;;
    --comparisons) need_value "$@"; COMPARISONS="$2"; shift 2 ;;
    --condition-col) need_value "$@"; CONDITION_COL="$2"; shift 2 ;;
    --group-col) need_value "$@"; GROUP_COL="$2"; shift 2 ;;
    --celltype-rename) need_value "$@"; CELLTYPE_RENAME="$2"; shift 2 ;;
    --assay) need_value "$@"; ASSAY="$2"; shift 2 ;;
    --layer) need_value "$@"; LAYER="$2"; shift 2 ;;
    --workers) need_value "$@"; WORKERS="$2"; shift 2 ;;
    --nboot) need_value "$@"; NBOOT="$2"; shift 2 ;;
    --min-cells) need_value "$@"; MIN_CELLS="$2"; shift 2 ;;
    --top-n) need_value "$@"; TOP_N="$2"; shift 2 ;;
    --plot-width) need_value "$@"; PLOT_WIDTH="$2"; shift 2 ;;
    --dpi) need_value "$@"; DPI="$2"; shift 2 ;;
    --force) FORCE=true; shift ;;
    --skip-analysis) RUN_ANALYSIS=false; shift ;;
    --fast-wilcox) need_value "$@"; FAST_WILCOX="$2"; shift 2 ;;
    --db-category) need_value "$@"; DB_CATEGORY="$2"; shift 2 ;;
    --dry-run) DRY_RUN=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

ENV_PREFIX="${ENV_PREFIX:-/Users/lyuwei/micromamba/envs/bioinfo}"
RSCRIPT_BIN="${ENV_PREFIX}/bin/Rscript"
LOG_DIR="${LOG_DIR:-logs}"
export PATH="${ENV_PREFIX}/bin:${PATH}"

print_settings() {
  cat <<SETTINGS
CellChat TopN + condition-specific information-flow settings
  input:              ${INPUT}
  outdir:             ${OUTDIR}
  conditions:         ${CONDITIONS}
  comparisons:        ${COMPARISONS}
  condition_col:      ${CONDITION_COL}
  group_col:          ${GROUP_COL}
  celltype_rename:    ${CELLTYPE_RENAME}
  assay/layer:        ${ASSAY}/${LAYER}
  workers/nboot:      ${WORKERS}/${NBOOT}
  min_cells:          ${MIN_CELLS}
  top_n:              ${TOP_N}
  run_analysis:       ${RUN_ANALYSIS}
  force:              ${FORCE}
  env_prefix:         ${ENV_PREFIX}
SETTINGS
}

print_settings
if [[ "$DRY_RUN" == true ]]; then
  exit 0
fi

if [[ ! -x "$RSCRIPT_BIN" ]]; then
  echo "Rscript not found or not executable: ${RSCRIPT_BIN}" >&2
  exit 1
fi

mkdir -p "$LOG_DIR" "$OUTDIR"

export FLOW_INPUT="$INPUT"
export FLOW_OUTDIR="$OUTDIR"
export FLOW_CONDITIONS="$CONDITIONS"
export FLOW_COMPARISONS="$COMPARISONS"
export FLOW_CONDITION_COL="$CONDITION_COL"
export FLOW_GROUP_COL="$GROUP_COL"
export FLOW_CELLTYPE_RENAME="$CELLTYPE_RENAME"
export FLOW_ASSAY="$ASSAY"
export FLOW_LAYER="$LAYER"
export FLOW_WORKERS="$WORKERS"
export FLOW_NBOOT="$NBOOT"
export FLOW_MIN_CELLS="$MIN_CELLS"
export FLOW_TOP_N="$TOP_N"
export FLOW_PLOT_WIDTH="$PLOT_WIDTH"
export FLOW_DPI="$DPI"
export FLOW_FORCE="$FORCE"
export FLOW_RUN_ANALYSIS="$RUN_ANALYSIS"
export FLOW_FAST_WILCOX="$FAST_WILCOX"
export FLOW_DB_CATEGORY="$DB_CATEGORY"

"$RSCRIPT_BIN" - <<'RS_FLOW' 2>&1 | tee "${LOG_DIR}/cellchat_top5_specific_flow_pipeline.log"
suppressPackageStartupMessages({
  library(CellChat)
  library(SeuratObject)
  library(Matrix)
  library(future)
  library(ggplot2)
  library(dplyr)
})

env_default <- function(name, default = "") {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

as_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
sanitize_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

cfg <- list(
  input = env_default("FLOW_INPUT", "gall0603.rds"),
  outdir = env_default("FLOW_OUTDIR", "results/cellchat"),
  assay = env_default("FLOW_ASSAY", "RNA"),
  layer = env_default("FLOW_LAYER", "data"),
  condition_col = env_default("FLOW_CONDITION_COL", "orig.ident"),
  group_col = env_default("FLOW_GROUP_COL", "celltype2"),
  celltype_rename = env_default("FLOW_CELLTYPE_RENAME", ""),
  conditions = trimws(strsplit(env_default("FLOW_CONDITIONS", "WT_ND,WT_CCD,KO_ND,KO_CCD"), ",", fixed = TRUE)[[1]]),
  comparisons = trimws(strsplit(env_default("FLOW_COMPARISONS", "WT_CCD_vs_WT_ND,KO_ND_vs_WT_ND"), ",", fixed = TRUE)[[1]]),
  min_cells = as.integer(env_default("FLOW_MIN_CELLS", "10")),
  workers = as.integer(env_default("FLOW_WORKERS", "4")),
  nboot = as.integer(env_default("FLOW_NBOOT", "100")),
  top_n = as.integer(env_default("FLOW_TOP_N", "5")),
  plot_width = as.numeric(env_default("FLOW_PLOT_WIDTH", "11")),
  dpi = as.integer(env_default("FLOW_DPI", "300")),
  force = as_bool(env_default("FLOW_FORCE", "false")),
  run_analysis = as_bool(env_default("FLOW_RUN_ANALYSIS", "true")),
  fast_wilcox = as_bool(env_default("FLOW_FAST_WILCOX", "true")),
  db_category = env_default("FLOW_DB_CATEGORY", "all")
)

object_dir <- file.path(cfg$outdir, "objects")
table_dir <- file.path(cfg$outdir, "tables")
plot_dir <- file.path(cfg$outdir, "plots", paste0("information_flow_top", cfg$top_n, "_plus_specific"))
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

get_assay_layer <- function(obj, assay_name, layer_name) {
  if (!assay_name %in% names(obj@assays)) {
    stop("Assay not found: ", assay_name, call. = FALSE)
  }
  assay <- obj@assays[[assay_name]]
  if ("layers" %in% slotNames(assay)) {
    if (!layer_name %in% names(assay@layers)) {
      stop("Layer '", layer_name, "' not found in assay '", assay_name, "'.", call. = FALSE)
    }
    mat <- assay@layers[[layer_name]]
  } else {
    if (!layer_name %in% slotNames(assay)) {
      stop("Slot '", layer_name, "' not found in assay '", assay_name, "'.", call. = FALSE)
    }
    mat <- slot(assay, layer_name)
  }
  if (is.null(rownames(mat))) rownames(mat) <- rownames(assay)
  if (is.null(colnames(mat))) colnames(mat) <- colnames(obj)
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("Expression matrix needs gene and cell names.", call. = FALSE)
  }
  mat
}

parse_rename_map <- function(rename_map) {
  if (!nzchar(rename_map)) {
    return(data.frame(old = character(), new = character(), stringsAsFactors = FALSE))
  }
  pairs <- trimws(strsplit(rename_map, ",", fixed = TRUE)[[1]])
  pairs <- pairs[nzchar(pairs)]
  if (length(pairs) == 0) {
    return(data.frame(old = character(), new = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, lapply(pairs, function(pair) {
    kv <- strsplit(pair, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2 || !nzchar(trimws(kv[[1]])) || !nzchar(trimws(kv[[2]]))) {
      stop("Invalid --celltype-rename entry: ", pair, call. = FALSE)
    }
    data.frame(old = trimws(kv[[1]]), new = trimws(kv[[2]]), stringsAsFactors = FALSE)
  }))
}

apply_celltype_rename <- function(meta, cfg) {
  rename_df <- parse_rename_map(cfg$celltype_rename)
  if (nrow(rename_df) == 0) return(meta)

  before <- as.character(meta[[cfg$group_col]])
  after <- before
  for (i in seq_len(nrow(rename_df))) {
    after[after == rename_df$old[[i]]] <- rename_df$new[[i]]
  }
  meta[[cfg$group_col]] <- factor(after)

  summary_df <- do.call(rbind, lapply(seq_len(nrow(rename_df)), function(i) {
    renamed <- before == rename_df$old[[i]]
    data.frame(
      old = rename_df$old[[i]],
      new = rename_df$new[[i]],
      n_cells_matching_old_before = sum(renamed, na.rm = TRUE),
      n_cells_renamed = sum(renamed & after == rename_df$new[[i]], na.rm = TRUE),
      n_cells_with_new_label_after = sum(after == rename_df$new[[i]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(summary_df, file.path(table_dir, "cellchat_celltype_rename_summary.csv"), row.names = FALSE)
  message("Applied celltype rename map: ", cfg$celltype_rename)
  print(summary_df)
  meta
}

run_one_condition <- function(condition, data_mat, meta, cfg, db_use) {
  out_rds <- file.path(object_dir, paste0(sanitize_name(condition), ".cellchat.rds"))
  if (file.exists(out_rds) && !cfg$force) {
    message("Reusing existing CellChat object: ", out_rds)
    return(readRDS(out_rds))
  }

  cells <- rownames(meta)[meta[[cfg$condition_col]] == condition]
  if (length(cells) == 0) stop("No cells found for condition: ", condition, call. = FALSE)

  meta_sub <- meta[cells, , drop = FALSE]
  meta_sub[[cfg$group_col]] <- droplevels(factor(meta_sub[[cfg$group_col]]))
  celltype_counts <- table(meta_sub[[cfg$group_col]])
  keep_types <- names(celltype_counts[celltype_counts >= cfg$min_cells])
  keep_cells <- rownames(meta_sub)[meta_sub[[cfg$group_col]] %in% keep_types]
  meta_sub <- meta_sub[keep_cells, , drop = FALSE]
  meta_sub[[cfg$group_col]] <- droplevels(factor(meta_sub[[cfg$group_col]]))
  data_sub <- data_mat[, rownames(meta_sub), drop = FALSE]

  message(
    "Running ", condition, ": ", ncol(data_sub), " cells, ",
    length(levels(meta_sub[[cfg$group_col]])), " cell types"
  )
  cellchat <- createCellChat(object = data_sub, meta = meta_sub, group.by = cfg$group_col, do.sparse = TRUE)
  cellchat@DB <- db_use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(
    cellchat,
    do.fast = cfg$fast_wilcox,
    min.cells = cfg$min_cells
  )
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean",
    nboot = cfg$nboot,
    seed.use = 1L,
    raw.use = TRUE,
    population.size = FALSE
  )
  cellchat <- filterCommunication(cellchat, min.cells = cfg$min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  saveRDS(cellchat, out_rds)
  message("Saved ", out_rds)
  cellchat
}

read_existing_objects <- function(conditions) {
  objects <- setNames(vector("list", length(conditions)), conditions)
  for (condition in conditions) {
    rds <- file.path(object_dir, paste0(sanitize_name(condition), ".cellchat.rds"))
    if (!file.exists(rds)) stop("Missing CellChat object: ", rds, call. = FALSE)
    objects[[condition]] <- readRDS(rds)
  }
  objects
}

if (cfg$run_analysis) {
  if (cfg$workers > 1) {
    options(future.globals.maxSize = 32 * 1024^3)
    future::plan("multisession", workers = cfg$workers)
  } else {
    future::plan("sequential")
  }

  message("Reading ", cfg$input)
  obj <- readRDS(cfg$input)
  meta <- obj@meta.data
  if (!identical(rownames(meta), colnames(obj))) {
    stop("meta.data rownames are not identical to Seurat cell names.", call. = FALSE)
  }
  required_cols <- c(cfg$condition_col, cfg$group_col)
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  meta <- apply_celltype_rename(meta, cfg)
  data_mat <- get_assay_layer(obj, cfg$assay, cfg$layer)

  missing_conditions <- setdiff(cfg$conditions, unique(as.character(meta[[cfg$condition_col]])))
  if (length(missing_conditions) > 0) {
    stop(
      "Missing conditions in ", cfg$condition_col, ": ",
      paste(missing_conditions, collapse = ", "),
      call. = FALSE
    )
  }

  data(CellChatDB.mouse)
  db_use <- CellChatDB.mouse
  if (!identical(tolower(cfg$db_category), "all")) {
    db_use <- subsetDB(CellChatDB.mouse, search = cfg$db_category)
  }

  objects <- setNames(vector("list", length(cfg$conditions)), cfg$conditions)
  for (condition in cfg$conditions) {
    objects[[condition]] <- run_one_condition(condition, data_mat, meta, cfg, db_use)
  }
  future::plan("sequential")
} else {
  objects <- read_existing_objects(cfg$conditions)
}

all_celltypes <- unique(unlist(lapply(objects, function(x) levels(x@idents))))
need_lift <- any(!vapply(objects, function(x) identical(levels(x@idents), all_celltypes), logical(1)))
if (need_lift) {
  message("Lifting CellChat objects to a shared cell type set: ", paste(all_celltypes, collapse = ", "))
  objects <- lapply(objects, function(x) liftCellChat(x, group.new = all_celltypes))
}
for (condition in cfg$conditions) {
  saveRDS(objects[[condition]], file.path(object_dir, paste0(sanitize_name(condition), ".cellchat.rds")))
}
merged <- mergeCellChat(objects, add.names = cfg$conditions, cell.prefix = TRUE)
saveRDS(merged, file.path(object_dir, "merged.cellchat.rds"))
saveRDS(objects, file.path(object_dir, "cellchat_object_list.rds"))

pathway_strength <- function(obj, condition) {
  arr <- obj@netP$prob
  if (is.null(arr) || length(arr) == 0 || length(dim(arr)) != 3) {
    return(data.frame(condition = condition, pathway_name = character(), strength = numeric()))
  }
  pathways <- dimnames(arr)[[3]]
  data.frame(
    condition = condition,
    pathway_name = pathways,
    strength = vapply(pathways, function(p) sum(arr[, , p], na.rm = TRUE), numeric(1)),
    stringsAsFactors = FALSE
  )
}

strength_long <- bind_rows(lapply(names(objects), function(condition) {
  pathway_strength(objects[[condition]], condition)
}))
write.csv(strength_long, file.path(table_dir, "cellchat_information_flow_all_pathway_strength_long.csv"), row.names = FALSE)

make_pairwise <- function(base_condition, compare_condition) {
  base <- strength_long %>%
    filter(condition == base_condition) %>%
    select(pathway_name, strength_base = strength)
  compare <- strength_long %>%
    filter(condition == compare_condition) %>%
    select(pathway_name, strength_compare = strength)
  full_join(base, compare, by = "pathway_name") %>%
    mutate(
      strength_base = ifelse(is.na(strength_base), 0, strength_base),
      strength_compare = ifelse(is.na(strength_compare), 0, strength_compare),
      base_condition = base_condition,
      compare_condition = compare_condition,
      comparison = paste(compare_condition, base_condition, sep = "_vs_"),
      delta_strength = strength_compare - strength_base,
      abs_delta_strength = abs(delta_strength),
      present_base = strength_base > 0,
      present_compare = strength_compare > 0,
      presence_class = case_when(
        present_base & present_compare ~ "shared",
        present_base & !present_compare ~ paste0(base_condition, "_only"),
        !present_base & present_compare ~ paste0(compare_condition, "_only"),
        TRUE ~ "absent"
      )
    ) %>%
    arrange(desc(abs_delta_strength), pathway_name)
}

parse_comparison <- function(comparison) {
  parts <- strsplit(comparison, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    stop("Invalid comparison. Expected compare_vs_base, got: ", comparison, call. = FALSE)
  }
  c(compare = parts[[1]], base = parts[[2]])
}

condition_index <- setNames(seq_along(cfg$conditions), cfg$conditions)
summary_rows <- list()

for (comparison_name in cfg$comparisons) {
  pair <- parse_comparison(comparison_name)
  compare_condition <- pair[["compare"]]
  base_condition <- pair[["base"]]
  if (!all(c(compare_condition, base_condition) %in% cfg$conditions)) {
    stop("Comparison contains condition not in --conditions: ", comparison_name, call. = FALSE)
  }

  pair_df <- make_pairwise(base_condition, compare_condition)
  top_pathways <- pair_df %>%
    filter(abs_delta_strength > 0) %>%
    slice_head(n = cfg$top_n) %>%
    pull(pathway_name)
  specific_pathways <- pair_df %>%
    filter(presence_class != "shared", presence_class != "absent") %>%
    pull(pathway_name)
  plot_pathways <- unique(c(top_pathways, specific_pathways))
  if (length(plot_pathways) == 0) {
    warning("No pathways selected for ", comparison_name)
    next
  }

  pathway_table <- pair_df %>%
    mutate(
      selected_top_delta = pathway_name %in% top_pathways,
      selected_condition_specific = pathway_name %in% specific_pathways,
      selected_for_plot = pathway_name %in% plot_pathways
    ) %>%
    arrange(
      desc(selected_for_plot),
      desc(selected_condition_specific),
      desc(selected_top_delta),
      desc(abs_delta_strength),
      pathway_name
    )
  write.csv(pathway_table, file.path(table_dir, paste0("cellchat_information_flow_top", cfg$top_n, "_plus_specific_", comparison_name, ".csv")), row.names = FALSE)
  write.csv(
    pathway_table %>%
      filter(selected_for_plot) %>%
      select(pathway_name, presence_class, strength_base, strength_compare, delta_strength, selected_top_delta, selected_condition_specific),
    file.path(table_dir, paste0("cellchat_information_flow_top", cfg$top_n, "_plus_specific_selected_pathways_", comparison_name, ".csv")),
    row.names = FALSE
  )

  plot_height <- max(7, 3.0 + 0.32 * length(plot_pathways))
  # Reversed CellChat colors: base group is teal, compare group is coral.
  p <- rankNet(
    merged,
    mode = "comparison",
    comparison = c(condition_index[[base_condition]], condition_index[[compare_condition]]),
    signaling = plot_pathways,
    stacked = TRUE,
    do.stat = FALSE,
    measure = "weight",
    color.use = c("#00BFC4", "#F8766D"),
    title = paste0(compare_condition, " vs ", base_condition, ": information flow Top", cfg$top_n, " + condition-specific")
  ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
      axis.title.x = element_text(face = "bold", size = 18),
      axis.title.y = element_blank(),
      axis.text.x = element_text(face = "bold", size = 15, color = "black"),
      axis.text.y = element_text(face = "bold", size = 15, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(face = "bold", size = 17, color = "black"),
      legend.key.size = unit(0.8, "cm"),
      plot.margin = margin(10, 18, 10, 14)
    )

  file_base <- file.path(plot_dir, paste0("cellchat_information_flow_top", cfg$top_n, "_plus_specific_", comparison_name, "_reversed_bold"))
  ggsave(paste0(file_base, ".pdf"), p, width = cfg$plot_width, height = plot_height, useDingbats = FALSE, limitsize = FALSE)
  ggsave(paste0(file_base, ".png"), p, width = cfg$plot_width, height = plot_height, dpi = cfg$dpi, limitsize = FALSE)

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    comparison = comparison_name,
    top_n = cfg$top_n,
    n_top_pathways = length(top_pathways),
    n_condition_specific = length(specific_pathways),
    n_selected_pathways = length(plot_pathways),
    includes_TGFb = "TGFb" %in% plot_pathways,
    plot_height_in = plot_height,
    stringsAsFactors = FALSE
  )
}

summary_df <- if (length(summary_rows) > 0) bind_rows(summary_rows) else data.frame()
write.csv(summary_df, file.path(table_dir, paste0("cellchat_information_flow_top", cfg$top_n, "_plus_specific_summary.csv")), row.names = FALSE)
plot_files <- list.files(plot_dir, pattern = "\\.(pdf|png)$", full.names = TRUE)
write.csv(
  data.frame(file = basename(plot_files), size_bytes = file.info(plot_files)$size, stringsAsFactors = FALSE),
  file.path(plot_dir, paste0("cellchat_information_flow_top", cfg$top_n, "_plus_specific_plot_files.csv")),
  row.names = FALSE
)

print(summary_df)
message("Wrote Top", cfg$top_n, " + condition-specific information-flow plots to ", plot_dir)
RS_FLOW
