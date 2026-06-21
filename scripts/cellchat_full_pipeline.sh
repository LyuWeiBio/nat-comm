#!/usr/bin/env bash
set -euo pipefail

# One-file CellChat workflow for gall0603.rds.
# It can optionally install dependencies, then run CellChat per orig.ident,
# summarize pairwise communication strength, and draw compact IL17 plots.

usage() {
  cat <<'USAGE'
Usage:
  ./cellchat_full_pipeline.sh [options]

Common options:
  --input FILE              Seurat RDS input. Default: gall0603.rds
  --outdir DIR              Output directory. Default: results/cellchat
  --conditions LIST         Comma-separated orig.ident values.
                            Default: WT_ND,WT_CCD,KO_ND,KO_CCD
  --condition-col COL       Metadata column for conditions. Default: orig.ident
  --group-col COL           Metadata column for cell types. Default: celltype2
  --celltype-rename MAP     Comma-separated old=new celltype rename map.
                            Default: CD4+ T=ab-T,STromal=Fibroblasts,Stromal=Fibroblasts
  --assay NAME              Seurat assay. Default: RNA
  --layer NAME              Seurat v5 layer or v4 slot. Default: data
  --workers N               CellChat future workers. Default: 4
  --nboot N                 CellChat computeCommunProb nboot. Default: 100
  --min-cells N             Minimum cells per celltype in each condition. Default: 10
  --focus NAME              Pathway/focus for summary and plots. Default: IL17

Run control:
  --install                 Install/verify dependencies in the bioinfo env first.
  --only-install            Install dependencies and exit.
  --force                   Recompute per-condition CellChat objects.
  --skip-analysis           Do not run/reuse CellChat analysis objects.
  --skip-summary            Do not generate pairwise CSV summaries.
  --skip-plot               Do not generate CellChat pathway plots.
  --skip-overview           Do not generate four-group overview tables/plots.
  --only-plot               Only generate plots from existing objects/tables.
  --only-overview           Only generate four-group overview tables/plots.
  --dry-run                 Print resolved settings and exit.

Environment:
  ENV_PREFIX                Micromamba env prefix.
                            Default: /Users/lyuwei/micromamba/envs/bioinfo
  ENV_NAME                  Micromamba env name for install. Default: bioinfo

Examples:
  ./cellchat_full_pipeline.sh
  ./cellchat_full_pipeline.sh --force
  ./cellchat_full_pipeline.sh --install
  ./cellchat_full_pipeline.sh --only-plot --focus IL17
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
CONDITION_COL="orig.ident"
GROUP_COL="celltype2"
CELLTYPE_RENAME="CD4+ T=ab-T,STromal=Fibroblasts,Stromal=Fibroblasts"
ASSAY="RNA"
LAYER="data"
WORKERS="4"
NBOOT="100"
MIN_CELLS="10"
FOCUS="IL17"
PVAL_THRESHOLD="0.05"
EPS="1e-12"
FAST_WILCOX="true"
DB_CATEGORY="all"
SOURCE=""
TARGET=""
BUBBLE_WIDTH="5.2"
BUBBLE_HEIGHT="3.8"
AGGREGATE_WIDTH="4.8"
AGGREGATE_HEIGHT="4.8"
CONTRIBUTION_WIDTH="4.8"
CONTRIBUTION_HEIGHT="3.0"
OVERVIEW_WIDTH="14"
OVERVIEW_HEIGHT="10"
OVERVIEW_DPI="300"
DPI="300"
INSTALL=false
RUN_ANALYSIS=true
RUN_SUMMARY=true
RUN_PLOT=true
RUN_OVERVIEW=true
FORCE=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) need_value "$@"; INPUT="$2"; shift 2 ;;
    --outdir) need_value "$@"; OUTDIR="$2"; shift 2 ;;
    --conditions) need_value "$@"; CONDITIONS="$2"; shift 2 ;;
    --condition-col) need_value "$@"; CONDITION_COL="$2"; shift 2 ;;
    --group-col) need_value "$@"; GROUP_COL="$2"; shift 2 ;;
    --celltype-rename) need_value "$@"; CELLTYPE_RENAME="$2"; shift 2 ;;
    --assay) need_value "$@"; ASSAY="$2"; shift 2 ;;
    --layer) need_value "$@"; LAYER="$2"; shift 2 ;;
    --workers) need_value "$@"; WORKERS="$2"; shift 2 ;;
    --nboot) need_value "$@"; NBOOT="$2"; shift 2 ;;
    --min-cells) need_value "$@"; MIN_CELLS="$2"; shift 2 ;;
    --focus) need_value "$@"; FOCUS="$2"; shift 2 ;;
    --pvalue-threshold) need_value "$@"; PVAL_THRESHOLD="$2"; shift 2 ;;
    --eps) need_value "$@"; EPS="$2"; shift 2 ;;
    --fast-wilcox) need_value "$@"; FAST_WILCOX="$2"; shift 2 ;;
    --db-category) need_value "$@"; DB_CATEGORY="$2"; shift 2 ;;
    --source) need_value "$@"; SOURCE="$2"; shift 2 ;;
    --target) need_value "$@"; TARGET="$2"; shift 2 ;;
    --bubble-width) need_value "$@"; BUBBLE_WIDTH="$2"; shift 2 ;;
    --bubble-height) need_value "$@"; BUBBLE_HEIGHT="$2"; shift 2 ;;
    --aggregate-width) need_value "$@"; AGGREGATE_WIDTH="$2"; shift 2 ;;
    --aggregate-height) need_value "$@"; AGGREGATE_HEIGHT="$2"; shift 2 ;;
    --contribution-width) need_value "$@"; CONTRIBUTION_WIDTH="$2"; shift 2 ;;
    --contribution-height) need_value "$@"; CONTRIBUTION_HEIGHT="$2"; shift 2 ;;
    --overview-width) need_value "$@"; OVERVIEW_WIDTH="$2"; shift 2 ;;
    --overview-height) need_value "$@"; OVERVIEW_HEIGHT="$2"; shift 2 ;;
    --overview-dpi) need_value "$@"; OVERVIEW_DPI="$2"; shift 2 ;;
    --dpi) need_value "$@"; DPI="$2"; shift 2 ;;
    --install) INSTALL=true; shift ;;
    --only-install) INSTALL=true; RUN_ANALYSIS=false; RUN_SUMMARY=false; RUN_PLOT=false; RUN_OVERVIEW=false; shift ;;
    --force) FORCE=true; shift ;;
    --skip-analysis) RUN_ANALYSIS=false; shift ;;
    --skip-summary) RUN_SUMMARY=false; shift ;;
    --skip-plot) RUN_PLOT=false; shift ;;
    --skip-overview) RUN_OVERVIEW=false; shift ;;
    --only-plot) RUN_ANALYSIS=false; RUN_SUMMARY=false; RUN_PLOT=true; RUN_OVERVIEW=false; shift ;;
    --only-overview) RUN_ANALYSIS=false; RUN_SUMMARY=false; RUN_PLOT=false; RUN_OVERVIEW=true; shift ;;
    --dry-run) DRY_RUN=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

ENV_PREFIX="${ENV_PREFIX:-/Users/lyuwei/micromamba/envs/bioinfo}"
ENV_NAME="${ENV_NAME:-bioinfo}"
R_BIN="${ENV_PREFIX}/bin/R"
RSCRIPT_BIN="${ENV_PREFIX}/bin/Rscript"
PYTHON_BIN="${ENV_PREFIX}/bin/python"
LOG_DIR="${LOG_DIR:-logs}"

export PATH="${ENV_PREFIX}/bin:${PATH}"

print_settings() {
  cat <<SETTINGS
CellChat full pipeline settings
  input:              ${INPUT}
  outdir:             ${OUTDIR}
  conditions:         ${CONDITIONS}
  condition_col:      ${CONDITION_COL}
  group_col:          ${GROUP_COL}
  celltype_rename:    ${CELLTYPE_RENAME}
  assay/layer:        ${ASSAY}/${LAYER}
  workers/nboot:      ${WORKERS}/${NBOOT}
  min_cells:          ${MIN_CELLS}
  focus:              ${FOCUS}
  install:            ${INSTALL}
  run_analysis:       ${RUN_ANALYSIS}
  run_summary:        ${RUN_SUMMARY}
  run_plot:           ${RUN_PLOT}
  run_overview:       ${RUN_OVERVIEW}
  force:              ${FORCE}
  env_prefix:         ${ENV_PREFIX}
SETTINGS
}

download_tarball() {
  local label="$1"
  local url="$2"
  local outfile="$3"

  rm -f "$outfile"
  echo "Downloading ${label}..."
  curl -L --fail --connect-timeout 20 --max-time 300 \
    -o "$outfile" \
    "https://gh-proxy.com/${url}" || \
  curl -L --fail --connect-timeout 20 --max-time 300 \
    -o "$outfile" \
    "$url"
}

package_installed() {
  local pkg="$1"
  "$RSCRIPT_BIN" -e "quit(status = ifelse(requireNamespace('${pkg}', quietly = TRUE), 0, 1))" >/dev/null 2>&1
}

install_dependencies() {
  local mamba_bin
  mamba_bin="${MAMBA_BIN:-$(command -v micromamba || true)}"
  if [[ -z "$mamba_bin" ]]; then
    echo "micromamba was not found in PATH. Set MAMBA_BIN or activate a shell with micromamba available." >&2
    exit 1
  fi

  "$mamba_bin" install -y -n "$ENV_NAME" -c conda-forge \
    r-seuratobject r-matrix r-future r-dplyr r-ggplot2

  "$mamba_bin" install -y -n "$ENV_NAME" -c conda-forge \
    numba llvmlite scanpy anndata h5py statsmodels seaborn plotnine pynndescent \
    umap-learn natsort

  "$PYTHON_BIN" -m pip install cellphonedb

  "$mamba_bin" install -y -n "$ENV_NAME" -c conda-forge -c bioconda \
    r-circlize bioconductor-complexheatmap bioconductor-biocneighbors \
    r-pbapply r-irlba r-ggalluvial r-svglite r-rspectra r-reticulate \
    r-sna r-fnn r-shape r-ggpubr r-ggnetwork r-plotly r-collapse \
    r-remotes r-biocmanager

  "$RSCRIPT_BIN" - <<'RS_INSTALL_CRAN'
options(repos = c(CRAN = "https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
if (!requireNamespace("NMF", quietly = TRUE)) {
  install.packages("NMF", Ncpus = 4)
}
RS_INSTALL_CRAN

  if ! package_installed CellChat; then
    rm -rf /tmp/CellChat-src /tmp/CellChat.tar.gz
    mkdir -p /tmp/CellChat-src
    download_tarball \
      "CellChat" \
      "https://github.com/jinworks/CellChat/archive/refs/heads/main.tar.gz" \
      "/tmp/CellChat.tar.gz"
    tar -xzf /tmp/CellChat.tar.gz -C /tmp/CellChat-src --strip-components=1
    "$R_BIN" CMD INSTALL /tmp/CellChat-src
  fi

  if ! package_installed presto; then
    rm -rf /tmp/presto-src /tmp/presto.tar.gz
    mkdir -p /tmp/presto-src
    download_tarball \
      "presto" \
      "https://github.com/immunogenomics/presto/archive/refs/heads/master.tar.gz" \
      "/tmp/presto.tar.gz"
    tar -xzf /tmp/presto.tar.gz -C /tmp/presto-src --strip-components=1
    "$R_BIN" CMD INSTALL /tmp/presto-src
  fi

  "$RSCRIPT_BIN" - <<'RS_VERIFY_INSTALL'
suppressPackageStartupMessages(library(CellChat))
suppressPackageStartupMessages(library(presto))
data(CellChatDB.mouse)
cat("CellChat ", as.character(packageVersion("CellChat")), "\n", sep = "")
cat("presto ", as.character(packageVersion("presto")), "\n", sep = "")
cat("CellChatDB.mouse interactions: ", nrow(CellChatDB.mouse$interaction), "\n", sep = "")
RS_VERIFY_INSTALL
}

print_settings

if [[ "$DRY_RUN" == true ]]; then
  exit 0
fi

mkdir -p "$LOG_DIR" "$OUTDIR"

if [[ "$INSTALL" == true ]]; then
  install_dependencies 2>&1 | tee "${LOG_DIR}/full_pipeline_install.log"
fi

if [[ ! -x "$RSCRIPT_BIN" ]]; then
  echo "Rscript not found or not executable: ${RSCRIPT_BIN}" >&2
  exit 1
fi

export PIPE_INPUT="$INPUT"
export PIPE_OUTDIR="$OUTDIR"
export PIPE_CONDITIONS="$CONDITIONS"
export PIPE_CONDITION_COL="$CONDITION_COL"
export PIPE_GROUP_COL="$GROUP_COL"
export PIPE_CELLTYPE_RENAME="$CELLTYPE_RENAME"
export PIPE_ASSAY="$ASSAY"
export PIPE_LAYER="$LAYER"
export PIPE_WORKERS="$WORKERS"
export PIPE_NBOOT="$NBOOT"
export PIPE_MIN_CELLS="$MIN_CELLS"
export PIPE_FOCUS="$FOCUS"
export PIPE_PVALUE_THRESHOLD="$PVAL_THRESHOLD"
export PIPE_EPS="$EPS"
export PIPE_FAST_WILCOX="$FAST_WILCOX"
export PIPE_DB_CATEGORY="$DB_CATEGORY"
export PIPE_FORCE="$FORCE"
export PIPE_SOURCE="$SOURCE"
export PIPE_TARGET="$TARGET"
export PIPE_BUBBLE_WIDTH="$BUBBLE_WIDTH"
export PIPE_BUBBLE_HEIGHT="$BUBBLE_HEIGHT"
export PIPE_AGGREGATE_WIDTH="$AGGREGATE_WIDTH"
export PIPE_AGGREGATE_HEIGHT="$AGGREGATE_HEIGHT"
export PIPE_CONTRIBUTION_WIDTH="$CONTRIBUTION_WIDTH"
export PIPE_CONTRIBUTION_HEIGHT="$CONTRIBUTION_HEIGHT"
export PIPE_OVERVIEW_WIDTH="$OVERVIEW_WIDTH"
export PIPE_OVERVIEW_HEIGHT="$OVERVIEW_HEIGHT"
export PIPE_OVERVIEW_DPI="$OVERVIEW_DPI"
export PIPE_DPI="$DPI"

if [[ "$RUN_ANALYSIS" == true ]]; then
  "$RSCRIPT_BIN" - <<'RS_ANALYSIS' 2>&1 | tee "${LOG_DIR}/full_pipeline_analysis.log"
suppressPackageStartupMessages({
  library(CellChat)
  library(SeuratObject)
  library(Matrix)
  library(future)
})

env_default <- function(name, default = "") {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

as_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
sanitize_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

cfg <- list(
  input = env_default("PIPE_INPUT", "gall0603.rds"),
  outdir = env_default("PIPE_OUTDIR", "results/cellchat"),
  assay = env_default("PIPE_ASSAY", "RNA"),
  layer = env_default("PIPE_LAYER", "data"),
  condition_col = env_default("PIPE_CONDITION_COL", "orig.ident"),
  group_col = env_default("PIPE_GROUP_COL", "celltype2"),
  celltype_rename = env_default("PIPE_CELLTYPE_RENAME", ""),
  conditions = env_default("PIPE_CONDITIONS", "WT_ND,WT_CCD,KO_ND,KO_CCD"),
  min_cells = as.integer(env_default("PIPE_MIN_CELLS", "10")),
  workers = as.integer(env_default("PIPE_WORKERS", "4")),
  nboot = as.integer(env_default("PIPE_NBOOT", "100")),
  force = as_bool(env_default("PIPE_FORCE", "false")),
  fast_wilcox = as_bool(env_default("PIPE_FAST_WILCOX", "true")),
  db_category = env_default("PIPE_DB_CATEGORY", "all")
)

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

run_one_condition <- function(condition, data_mat, meta, cfg, db_use, object_dir) {
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

parse_rename_map <- function(rename_map) {
  if (!nzchar(rename_map)) {
    return(data.frame(old = character(), new = character(), stringsAsFactors = FALSE))
  }
  pairs <- trimws(strsplit(rename_map, ",", fixed = TRUE)[[1]])
  pairs <- pairs[nzchar(pairs)]
  if (length(pairs) == 0) {
    return(data.frame(old = character(), new = character(), stringsAsFactors = FALSE))
  }
  parsed <- lapply(pairs, function(pair) {
    kv <- strsplit(pair, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2 || !nzchar(trimws(kv[[1]])) || !nzchar(trimws(kv[[2]]))) {
      stop("Invalid --celltype-rename entry: ", pair, call. = FALSE)
    }
    data.frame(old = trimws(kv[[1]]), new = trimws(kv[[2]]), stringsAsFactors = FALSE)
  })
  do.call(rbind, parsed)
}

apply_celltype_rename <- function(meta, cfg, table_dir) {
  rename_df <- parse_rename_map(cfg$celltype_rename)
  if (nrow(rename_df) == 0) {
    return(meta)
  }

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

  condition_df <- do.call(rbind, lapply(seq_len(nrow(rename_df)), function(i) {
    old <- rename_df$old[[i]]
    new <- rename_df$new[[i]]
    conditions <- sort(unique(as.character(meta[[cfg$condition_col]])))
    data.frame(
      old = old,
      new = new,
      condition = conditions,
      n_cells_matching_old_before = as.integer(table(factor(as.character(meta[[cfg$condition_col]])[before == old], levels = conditions))),
      n_cells_renamed = as.integer(table(factor(as.character(meta[[cfg$condition_col]])[before == old & after == new], levels = conditions))),
      n_cells_with_new_label_after = as.integer(table(factor(as.character(meta[[cfg$condition_col]])[after == new], levels = conditions))),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(condition_df, file.path(table_dir, "cellchat_celltype_rename_by_condition.csv"), row.names = FALSE)

  message("Applied celltype rename map: ", cfg$celltype_rename)
  print(summary_df)
  meta
}

object_dir <- file.path(cfg$outdir, "objects")
table_dir <- file.path(cfg$outdir, "tables")
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

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

meta <- apply_celltype_rename(meta, cfg, table_dir)
data_mat <- get_assay_layer(obj, cfg$assay, cfg$layer)
conditions <- trimws(strsplit(cfg$conditions, ",", fixed = TRUE)[[1]])
missing_conditions <- setdiff(conditions, unique(as.character(meta[[cfg$condition_col]])))
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

run_summary <- do.call(rbind, lapply(conditions, function(condition) {
  cells <- rownames(meta)[meta[[cfg$condition_col]] == condition]
  tab_all <- table(droplevels(factor(meta[cells, cfg$group_col])))
  tab_all <- sort(tab_all[tab_all > 0], decreasing = TRUE)
  tab_keep <- tab_all[tab_all >= cfg$min_cells]
  tab_drop <- tab_all[tab_all < cfg$min_cells]
  data.frame(
    condition = condition,
    n_cells = length(cells),
    min_cells = cfg$min_cells,
    n_celltypes_raw = length(tab_all),
    n_celltypes_analyzed = length(tab_keep),
    celltype_counts_raw = paste(paste(names(tab_all), as.integer(tab_all), sep = ":"), collapse = ";"),
    celltype_counts_analyzed = paste(paste(names(tab_keep), as.integer(tab_keep), sep = ":"), collapse = ";"),
    celltype_counts_dropped = paste(paste(names(tab_drop), as.integer(tab_drop), sep = ":"), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
write.csv(run_summary, file.path(table_dir, "cellchat_run_input_summary.csv"), row.names = FALSE)

object_list <- setNames(vector("list", length(conditions)), conditions)
for (condition in conditions) {
  object_list[[condition]] <- run_one_condition(condition, data_mat, meta, cfg, db_use, object_dir)
}

all_celltypes <- unique(unlist(lapply(object_list, function(x) levels(x@idents))))
need_lift <- any(!vapply(object_list, function(x) identical(levels(x@idents), all_celltypes), logical(1)))
if (need_lift) {
  message("Lifting CellChat objects to a shared cell type set: ", paste(all_celltypes, collapse = ", "))
  object_list <- lapply(object_list, function(x) liftCellChat(x, group.new = all_celltypes))
  for (condition in conditions) {
    saveRDS(object_list[[condition]], file.path(object_dir, paste0(sanitize_name(condition), ".cellchat.rds")))
  }
}

merged <- mergeCellChat(object_list, add.names = conditions, cell.prefix = TRUE)
saveRDS(merged, file.path(object_dir, "merged.cellchat.rds"))
saveRDS(object_list, file.path(object_dir, "cellchat_object_list.rds"))
message("Saved merged CellChat object and object list under ", object_dir)

session <- capture.output(sessionInfo())
writeLines(session, file.path(cfg$outdir, "sessionInfo.txt"))
future::plan("sequential")
message("CellChat analysis done.")
RS_ANALYSIS
fi

if [[ "$RUN_SUMMARY" == true ]]; then
  "$RSCRIPT_BIN" - <<'RS_SUMMARY' 2>&1 | tee "${LOG_DIR}/full_pipeline_summary.log"
suppressPackageStartupMessages({
  library(CellChat)
  library(dplyr)
})

env_default <- function(name, default = "") {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

sanitize_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

array_to_long <- function(arr, value_name, third_name) {
  if (length(dim(arr)) != 3) {
    stop("Expected a 3D CellChat probability array.", call. = FALSE)
  }
  out <- as.data.frame(as.table(arr), stringsAsFactors = FALSE)
  names(out) <- c("source", "target", third_name, value_name)
  out[[value_name]] <- as.numeric(out[[value_name]])
  out
}

extract_lr <- function(cellchat, condition, pvalue_threshold) {
  prob <- array_to_long(cellchat@net$prob, "prob", "interaction_name")
  pval <- array_to_long(cellchat@net$pval, "pval", "interaction_name")
  lr <- prob %>%
    left_join(pval, by = c("source", "target", "interaction_name")) %>%
    mutate(condition = condition, .before = 1)

  lr_meta <- cellchat@LR$LRsig
  lr_meta$interaction_name <- rownames(lr_meta)
  keep_meta <- intersect(
    c("interaction_name", "pathway_name", "ligand", "receptor", "annotation", "evidence"),
    colnames(lr_meta)
  )
  lr <- lr %>% left_join(lr_meta[, keep_meta, drop = FALSE], by = "interaction_name")
  lr %>% filter(prob > 0)
}

extract_pathway <- function(cellchat, condition) {
  if (is.null(cellchat@netP$prob) || length(cellchat@netP$prob) == 0) {
    return(data.frame())
  }
  path <- array_to_long(cellchat@netP$prob, "prob", "pathway_name")
  path %>%
    mutate(condition = condition, .before = 1) %>%
    filter(prob > 0)
}

make_pairwise <- function(df, conditions, keys, eps) {
  if (nrow(df) == 0) {
    return(data.frame())
  }
  pair_index <- combn(conditions, 2, simplify = FALSE)
  out <- lapply(pair_index, function(pair) {
    a <- pair[[1]]
    b <- pair[[2]]
    da <- df %>%
      filter(condition == a) %>%
      select(all_of(keys), prob_a = prob, any_of("pval"))
    if ("pval" %in% colnames(da)) da <- rename(da, pval_a = pval)
    db <- df %>%
      filter(condition == b) %>%
      select(all_of(keys), prob_b = prob, any_of("pval"))
    if ("pval" %in% colnames(db)) db <- rename(db, pval_b = pval)
    full_join(da, db, by = keys) %>%
      mutate(
        condition_a = a,
        condition_b = b,
        comparison = paste(a, b, sep = "_vs_"),
        prob_a = ifelse(is.na(prob_a), 0, prob_a),
        prob_b = ifelse(is.na(prob_b), 0, prob_b),
        delta_prob = prob_b - prob_a,
        ratio = (prob_b + eps) / (prob_a + eps),
        log2_ratio = log2(ratio),
        abs_delta_prob = abs(delta_prob),
        .before = 1
      ) %>%
      arrange(desc(abs_delta_prob))
  })
  bind_rows(out)
}

focus_filter <- function(df, focus, columns) {
  if (nrow(df) == 0) return(rep(FALSE, 0))
  keep <- rep(FALSE, nrow(df))
  for (col in intersect(columns, colnames(df))) {
    keep <- keep | grepl(focus, df[[col]], ignore.case = TRUE)
  }
  keep
}

outdir <- file.path(env_default("PIPE_OUTDIR", "results/cellchat"), "tables")
objects_dir <- file.path(env_default("PIPE_OUTDIR", "results/cellchat"), "objects")
conditions <- trimws(strsplit(env_default("PIPE_CONDITIONS", "WT_ND,WT_CCD,KO_ND,KO_CCD"), ",", fixed = TRUE)[[1]])
focus <- env_default("PIPE_FOCUS", "IL17")
pvalue_threshold <- as.numeric(env_default("PIPE_PVALUE_THRESHOLD", "0.05"))
eps <- as.numeric(env_default("PIPE_EPS", "1e-12"))
focus_file <- sanitize_name(focus)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

objects <- setNames(vector("list", length(conditions)), conditions)
for (condition in conditions) {
  rds <- file.path(objects_dir, paste0(sanitize_name(condition), ".cellchat.rds"))
  if (!file.exists(rds)) stop("Missing CellChat object: ", rds, call. = FALSE)
  objects[[condition]] <- readRDS(rds)
}

lr_long <- bind_rows(lapply(names(objects), function(condition) {
  extract_lr(objects[[condition]], condition, pvalue_threshold)
}))
pathway_long <- bind_rows(lapply(names(objects), function(condition) {
  extract_pathway(objects[[condition]], condition)
}))

write.csv(lr_long, file.path(outdir, "cellchat_lr_strength_long.csv"), row.names = FALSE)
write.csv(pathway_long, file.path(outdir, "cellchat_pathway_strength_long.csv"), row.names = FALSE)

lr_keys <- c("source", "target", "interaction_name", "pathway_name", "ligand", "receptor", "annotation")
lr_keys <- intersect(lr_keys, colnames(lr_long))
pathway_keys <- c("source", "target", "pathway_name")

lr_pairwise <- make_pairwise(lr_long, conditions, lr_keys, eps)
pathway_pairwise <- make_pairwise(pathway_long, conditions, pathway_keys, eps)

write.csv(lr_pairwise, file.path(outdir, "cellchat_lr_pairwise_comparisons.csv"), row.names = FALSE)
write.csv(pathway_pairwise, file.path(outdir, "cellchat_pathway_pairwise_comparisons.csv"), row.names = FALSE)

lr_focus <- lr_pairwise[focus_filter(
  lr_pairwise,
  focus,
  c("pathway_name", "interaction_name", "ligand", "receptor")
), , drop = FALSE]
pathway_focus <- pathway_pairwise[focus_filter(
  pathway_pairwise,
  focus,
  c("pathway_name")
), , drop = FALSE]

write.csv(lr_focus, file.path(outdir, paste0("cellchat_", focus_file, "_lr_pairwise_comparisons.csv")), row.names = FALSE)
write.csv(pathway_focus, file.path(outdir, paste0("cellchat_", focus_file, "_pathway_pairwise_comparisons.csv")), row.names = FALSE)

summary_df <- data.frame(
  table = c("lr_long", "pathway_long", "lr_pairwise", "pathway_pairwise", "lr_focus", "pathway_focus"),
  rows = c(nrow(lr_long), nrow(pathway_long), nrow(lr_pairwise), nrow(pathway_pairwise), nrow(lr_focus), nrow(pathway_focus)),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(outdir, "cellchat_comparison_table_summary.csv"), row.names = FALSE)
print(summary_df)
RS_SUMMARY
fi

if [[ "$RUN_PLOT" == true ]]; then
  "$RSCRIPT_BIN" - <<'RS_PLOT' 2>&1 | tee "${LOG_DIR}/full_pipeline_plot.log"
suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(grid)
})

env_default <- function(name, default = "") {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

sanitize_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

plot_log <- data.frame(
  plot = character(),
  status = character(),
  message = character(),
  stringsAsFactors = FALSE
)

log_plot <- function(plot, status, message = "") {
  assign(
    "plot_log",
    rbind(get("plot_log", envir = .GlobalEnv), data.frame(plot = plot, status = status, message = message)),
    envir = .GlobalEnv
  )
}

save_gg_safe <- function(plot_expr, filename, width, height, dpi) {
  label <- basename(filename)
  p <- try(force(plot_expr), silent = TRUE)
  if (inherits(p, "try-error")) {
    log_plot(label, "skipped", as.character(p))
    return(invisible(FALSE))
  }
  ggsave(paste0(filename, ".pdf"), plot = p, width = width, height = height, useDingbats = FALSE)
  ggsave(paste0(filename, ".png"), plot = p, width = width, height = height, dpi = dpi)
  rows <- if (!is.null(p$data)) paste0("; rows=", nrow(p$data)) else ""
  log_plot(label, "written", rows)
  invisible(TRUE)
}

draw_cellchat_plot <- function(x) {
  if (inherits(x, "ggplot")) {
    print(x)
  } else if (inherits(x, "Heatmap") || inherits(x, "HeatmapList")) {
    ComplexHeatmap::draw(x)
  } else if (inherits(x, "grob") || inherits(x, "gTree")) {
    grid::grid.draw(x)
  }
}

save_base_safe <- function(plot_expr, filename, width, height, dpi) {
  label <- basename(filename)
  expr <- substitute(plot_expr)
  env <- parent.frame()
  ok <- TRUE
  err_msg <- ""
  pdf_file <- paste0(filename, ".pdf")
  png_file <- paste0(filename, ".png")

  pdf(pdf_file, width = width, height = height, useDingbats = FALSE)
  res <- try(eval(expr, envir = env), silent = TRUE)
  if (inherits(res, "try-error")) {
    ok <- FALSE
    err_msg <- as.character(res)
  }
  invisible(dev.off())

  if (ok) {
    png(png_file, width = width, height = height, units = "in", res = dpi)
    res <- try(eval(expr, envir = env), silent = TRUE)
    if (inherits(res, "try-error")) {
      ok <- FALSE
      err_msg <- as.character(res)
    }
    invisible(dev.off())
  }

  if (ok) {
    log_plot(label, "written", "")
  } else {
    if (file.exists(pdf_file)) unlink(pdf_file)
    if (file.exists(png_file)) unlink(png_file)
    log_plot(label, "skipped", err_msg)
  }
  invisible(ok)
}

base_outdir <- env_default("PIPE_OUTDIR", "results/cellchat")
objects_dir <- file.path(base_outdir, "objects")
tables_dir <- file.path(base_outdir, "tables")
conditions <- trimws(strsplit(env_default("PIPE_CONDITIONS", "WT_ND,WT_CCD,KO_ND,KO_CCD"), ",", fixed = TRUE)[[1]])
signaling <- env_default("PIPE_FOCUS", "IL17")
signaling_file <- sanitize_name(signaling)
outdir <- file.path(base_outdir, "plots", signaling_file)
width <- as.numeric(env_default("PIPE_BUBBLE_WIDTH", "5.2"))
height <- as.numeric(env_default("PIPE_BUBBLE_HEIGHT", "3.8"))
aggregate_width <- as.numeric(env_default("PIPE_AGGREGATE_WIDTH", "4.8"))
aggregate_height <- as.numeric(env_default("PIPE_AGGREGATE_HEIGHT", "4.8"))
contribution_width <- as.numeric(env_default("PIPE_CONTRIBUTION_WIDTH", "4.8"))
contribution_height <- as.numeric(env_default("PIPE_CONTRIBUTION_HEIGHT", "3.0"))
dpi <- as.integer(env_default("PIPE_DPI", "300"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

compact_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 10),
  axis.text.x = element_text(size = 8, hjust = 1, vjust = 1),
  axis.text.y = element_text(size = 8),
  axis.title = element_text(size = 9),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 7),
  plot.margin = margin(3, 3, 3, 3)
)

objects <- setNames(vector("list", length(conditions)), conditions)
for (condition in conditions) {
  rds <- file.path(objects_dir, paste0(sanitize_name(condition), ".cellchat.rds"))
  if (!file.exists(rds)) stop("Missing CellChat object: ", rds, call. = FALSE)
  objects[[condition]] <- readRDS(rds)
}

merged_file <- file.path(objects_dir, "merged.cellchat.rds")
if (!file.exists(merged_file)) stop("Missing merged CellChat object: ", merged_file, call. = FALSE)
merged <- readRDS(merged_file)

focus_sources <- NULL
focus_targets <- NULL
source_arg <- env_default("PIPE_SOURCE", "")
target_arg <- env_default("PIPE_TARGET", "")
if (nzchar(source_arg)) focus_sources <- trimws(strsplit(source_arg, ",", fixed = TRUE)[[1]])
if (nzchar(target_arg)) focus_targets <- trimws(strsplit(target_arg, ",", fixed = TRUE)[[1]])

focus_table <- file.path(tables_dir, paste0("cellchat_", signaling_file, "_pathway_pairwise_comparisons.csv"))
if (file.exists(focus_table) && (is.null(focus_sources) || is.null(focus_targets))) {
  focus_df <- read.csv(focus_table, check.names = FALSE)
  if (nrow(focus_df) > 0) {
    if (is.null(focus_sources) && "source" %in% colnames(focus_df)) {
      focus_sources <- unique(focus_df$source)
    }
    if (is.null(focus_targets) && "target" %in% colnames(focus_df)) {
      focus_targets <- unique(focus_df$target)
    }
  }
}

celltypes <- unique(unlist(lapply(objects, function(x) levels(x@idents))))
if (is.null(focus_sources)) focus_sources <- celltypes
if (is.null(focus_targets)) focus_targets <- celltypes

message("Plotting ", signaling, " sources: ", paste(focus_sources, collapse = ", "))
message("Plotting ", signaling, " targets: ", paste(focus_targets, collapse = ", "))

pathway_strength_rows <- lapply(conditions, function(condition) {
  obj <- objects[[condition]]
  pathways <- dimnames(obj@netP$prob)[[3]]
  if (is.null(pathways) || !signaling %in% pathways) {
    return(data.frame(
      condition = condition,
      total_strength = 0,
      source_target = "none",
      n_nonzero_pairs = 0,
      stringsAsFactors = FALSE
    ))
  }
  mat <- obj@netP$prob[, , signaling, drop = FALSE][, , 1]
  idx <- which(mat > 0, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(data.frame(
      condition = condition,
      total_strength = 0,
      source_target = "none",
      n_nonzero_pairs = 0,
      stringsAsFactors = FALSE
    ))
  }
  pairs <- paste(rownames(mat)[idx[, 1]], colnames(mat)[idx[, 2]], sep = " -> ")
  data.frame(
    condition = condition,
    total_strength = sum(mat[idx]),
    source_target = paste(unique(pairs), collapse = "; "),
    n_nonzero_pairs = nrow(idx),
    stringsAsFactors = FALSE
  )
})
pathway_strength_df <- do.call(rbind, pathway_strength_rows)
pathway_strength_df$condition <- factor(pathway_strength_df$condition, levels = conditions)
write.csv(
  pathway_strength_df,
  file.path(tables_dir, paste0("cellchat_", signaling_file, "_four_group_pathway_strength.csv")),
  row.names = FALSE
)

subtitle_pairs <- unique(unlist(strsplit(
  paste(pathway_strength_df$source_target[pathway_strength_df$source_target != "none"], collapse = "; "),
  "; ",
  fixed = TRUE
)))
subtitle_pairs <- subtitle_pairs[nzchar(subtitle_pairs)]
subtitle_text <- if (length(subtitle_pairs) == 0) {
  "CellChat pathway probability, summed across source-target pairs"
} else {
  paste0(
    "CellChat pathway probability, summed across source-target pairs; nonzero pair",
    ifelse(length(subtitle_pairs) > 1, "s: ", ": "),
    paste(subtitle_pairs, collapse = "; ")
  )
}
if (nchar(subtitle_text) > 140) {
  subtitle_text <- "CellChat pathway probability, summed across source-target pairs"
}
max_strength <- max(pathway_strength_df$total_strength, na.rm = TRUE)
plot_scale <- ifelse(max_strength > 0, max_strength, 1)
pathway_strength_df$label_y <- ifelse(
  pathway_strength_df$total_strength > 0,
  pathway_strength_df$total_strength + plot_scale * 0.06,
  plot_scale * 0.04
)

save_gg_safe(
  ggplot(pathway_strength_df, aes(x = condition, y = total_strength, fill = condition)) +
    geom_col(width = 0.62, color = "grey25", linewidth = 0.25) +
    geom_point(size = 2.6, color = "black") +
    geom_text(aes(y = label_y, label = sprintf("%.4g", total_strength)), size = 4) +
    scale_fill_manual(values = c(WT_ND = "#7A7A7A", WT_CCD = "#D95F02", KO_ND = "#1B9E77", KO_CCD = "#7570B3")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = paste0(signaling, " pathway signaling strength across orig.ident groups"),
      subtitle = subtitle_text,
      x = NULL,
      y = paste0(signaling, " pathway strength")
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title.y = element_text(size = 12),
      plot.margin = margin(8, 12, 8, 12)
    ),
  file.path(outdir, paste0(signaling_file, "_four_group_pathway_strength")),
  7,
  4.5,
  dpi
)

pair_index <- combn(seq_along(conditions), 2, simplify = FALSE)
for (pair in pair_index) {
  a <- pair[[1]]
  b <- pair[[2]]
  comp_name <- paste0(conditions[[a]], "_vs_", conditions[[b]])

  save_gg_safe(
    netVisual_bubble(
      merged,
      sources.use = focus_sources,
      targets.use = focus_targets,
      signaling = signaling,
      comparison = c(a, b),
      max.dataset = b,
      remove.isolate = TRUE,
      angle.x = 45,
      font.size = 8,
      font.size.title = 10,
      dot.size.min = 2,
      dot.size.max = 5,
      title.name = paste0(signaling, " increased in ", conditions[[b]], " vs ", conditions[[a]])
    ) + compact_theme,
    file.path(outdir, paste0(signaling_file, "_", comp_name, "_increased_in_", conditions[[b]])),
    width,
    height,
    dpi
  )

  save_gg_safe(
    netVisual_bubble(
      merged,
      sources.use = focus_sources,
      targets.use = focus_targets,
      signaling = signaling,
      comparison = c(a, b),
      max.dataset = a,
      remove.isolate = TRUE,
      angle.x = 45,
      font.size = 8,
      font.size.title = 10,
      dot.size.min = 2,
      dot.size.max = 5,
      title.name = paste0(signaling, " decreased in ", conditions[[b]], " vs ", conditions[[a]])
    ) + compact_theme,
    file.path(outdir, paste0(signaling_file, "_", comp_name, "_decreased_in_", conditions[[b]])),
    width,
    height,
    dpi
  )
}

for (condition in conditions) {
  obj <- objects[[condition]]
  group_size <- as.numeric(table(obj@idents))
  names(group_size) <- names(table(obj@idents))

  save_base_safe(
    {
      p <- netVisual_aggregate(
        obj,
        signaling = signaling,
        layout = "circle",
        remove.isolate = FALSE,
        vertex.weight = group_size,
        vertex.size.max = 16,
        edge.width.max = 8,
        pt.title = 8,
        title.space = 2,
        vertex.label.cex = 0.55,
        signaling.name = paste0(signaling, " - ", condition)
      )
      draw_cellchat_plot(p)
    },
    file.path(outdir, paste0(signaling_file, "_", condition, "_aggregate_circle")),
    aggregate_width,
    aggregate_height,
    dpi
  )

  save_base_safe(
    {
      p <- netAnalysis_contribution(
        obj,
        signaling = signaling,
        title = paste0(signaling, " LR contribution - ", condition)
      )
      draw_cellchat_plot(p)
    },
    file.path(outdir, paste0(signaling_file, "_", condition, "_LR_contribution")),
    contribution_width,
    contribution_height,
    dpi
  )
}

plot_files <- list.files(outdir, pattern = paste0("^", signaling_file, ".*\\.(pdf|png)$"), full.names = TRUE)
file_summary <- data.frame(
  file = basename(plot_files),
  size_bytes = file.info(plot_files)$size,
  stringsAsFactors = FALSE
)
write.csv(file_summary, file.path(outdir, paste0(signaling_file, "_plot_files.csv")), row.names = FALSE)
write.csv(plot_log, file.path(outdir, paste0(signaling_file, "_plot_log.csv")), row.names = FALSE)
message("Wrote ", nrow(file_summary), " plot files to ", outdir)
message("Plot log written to ", file.path(outdir, paste0(signaling_file, "_plot_log.csv")))
RS_PLOT
fi

if [[ "$RUN_OVERVIEW" == true ]]; then
  "$RSCRIPT_BIN" - <<'RS_OVERVIEW' 2>&1 | tee "${LOG_DIR}/full_pipeline_overview.log"
suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(grid)
})

env_default <- function(name, default = "") {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

sanitize_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

plot_log <- data.frame(
  plot = character(),
  status = character(),
  message = character(),
  stringsAsFactors = FALSE
)

log_plot <- function(plot, status, message = "") {
  assign(
    "plot_log",
    rbind(get("plot_log", envir = .GlobalEnv), data.frame(plot = plot, status = status, message = message)),
    envir = .GlobalEnv
  )
}

save_gg_safe <- function(plot_expr, filename, width, height, dpi) {
  label <- basename(filename)
  p <- try(force(plot_expr), silent = TRUE)
  if (inherits(p, "try-error")) {
    log_plot(label, "skipped", as.character(p))
    return(invisible(FALSE))
  }
  ggsave(paste0(filename, ".pdf"), plot = p, width = width, height = height, useDingbats = FALSE, limitsize = FALSE)
  ggsave(paste0(filename, ".png"), plot = p, width = width, height = height, dpi = dpi, limitsize = FALSE)
  log_plot(label, "written", "")
  invisible(TRUE)
}

save_base_safe <- function(plot_expr, filename, width, height, dpi) {
  label <- basename(filename)
  expr <- substitute(plot_expr)
  env <- parent.frame()
  ok <- TRUE
  err_msg <- ""
  pdf_file <- paste0(filename, ".pdf")
  png_file <- paste0(filename, ".png")

  pdf(pdf_file, width = width, height = height, useDingbats = FALSE)
  res <- try(eval(expr, envir = env), silent = TRUE)
  if (inherits(res, "try-error")) {
    ok <- FALSE
    err_msg <- as.character(res)
  }
  invisible(dev.off())

  if (ok) {
    png(png_file, width = width, height = height, units = "in", res = dpi)
    res <- try(eval(expr, envir = env), silent = TRUE)
    if (inherits(res, "try-error")) {
      ok <- FALSE
      err_msg <- as.character(res)
    }
    invisible(dev.off())
  }

  if (ok) {
    log_plot(label, "written", "")
  } else {
    if (file.exists(pdf_file)) unlink(pdf_file)
    if (file.exists(png_file)) unlink(png_file)
    log_plot(label, "skipped", err_msg)
  }
  invisible(ok)
}

mat_to_long <- function(mat, condition, measure, celltypes = NULL) {
  if (!is.null(celltypes)) {
    grid <- expand.grid(source = celltypes, target = celltypes, stringsAsFactors = FALSE)
    value_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    names(value_df) <- c("source", "target", "value")
    out <- merge(grid, value_df, by = c("source", "target"), all.x = TRUE, sort = FALSE)
    out$value[is.na(out$value)] <- 0
  } else {
    out <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    names(out) <- c("source", "target", "value")
  }
  out$condition <- condition
  out$measure <- measure
  out$value <- as.numeric(out$value)
  out[, c("condition", "measure", "source", "target", "value")]
}

array_to_long <- function(arr, condition) {
  if (is.null(arr) || length(arr) == 0 || length(dim(arr)) != 3) {
    return(data.frame())
  }
  out <- as.data.frame(as.table(arr), stringsAsFactors = FALSE)
  names(out) <- c("source", "target", "pathway_name", "value")
  out$condition <- condition
  out$value <- as.numeric(out$value)
  out[, c("condition", "source", "target", "pathway_name", "value")]
}

make_wide_difference <- function(df, keys, conditions, value_name = "value") {
  if (nrow(df) == 0) return(data.frame())
  key_df <- unique(df[, keys, drop = FALSE])
  out <- key_df
  for (condition in conditions) {
    sub <- df[df$condition == condition, c(keys, value_name), drop = FALSE]
    names(sub)[ncol(sub)] <- condition
    out <- merge(out, sub, by = keys, all.x = TRUE, sort = FALSE)
  }
  for (condition in conditions) {
    if (!condition %in% colnames(out)) out[[condition]] <- 0
    out[[condition]][is.na(out[[condition]])] <- 0
  }
  values <- as.matrix(out[, conditions, drop = FALSE])
  max_idx <- max.col(values, ties.method = "first")
  min_idx <- max.col(-values, ties.method = "first")
  out$max_condition <- conditions[max_idx]
  out$max_value <- values[cbind(seq_len(nrow(values)), max_idx)]
  out$min_condition <- conditions[min_idx]
  out$min_value <- values[cbind(seq_len(nrow(values)), min_idx)]
  out$range_value <- out$max_value - out$min_value
  out <- out[order(-out$range_value, out$max_condition), , drop = FALSE]
  rownames(out) <- NULL
  out
}

base_outdir <- env_default("PIPE_OUTDIR", "results/cellchat")
objects_dir <- file.path(base_outdir, "objects")
tables_dir <- file.path(base_outdir, "tables")
outdir <- file.path(base_outdir, "plots", "overview")
conditions <- trimws(strsplit(env_default("PIPE_CONDITIONS", "WT_ND,WT_CCD,KO_ND,KO_CCD"), ",", fixed = TRUE)[[1]])
width <- as.numeric(env_default("PIPE_OVERVIEW_WIDTH", "14"))
height <- as.numeric(env_default("PIPE_OVERVIEW_HEIGHT", "10"))
dpi <- as.integer(env_default("PIPE_OVERVIEW_DPI", "300"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

objects <- setNames(vector("list", length(conditions)), conditions)
for (condition in conditions) {
  rds <- file.path(objects_dir, paste0(sanitize_name(condition), ".cellchat.rds"))
  if (!file.exists(rds)) stop("Missing CellChat object: ", rds, call. = FALSE)
  objects[[condition]] <- readRDS(rds)
}

merged_file <- file.path(objects_dir, "merged.cellchat.rds")
if (!file.exists(merged_file)) stop("Missing merged CellChat object: ", merged_file, call. = FALSE)
merged <- readRDS(merged_file)

celltypes <- unique(unlist(lapply(objects, function(x) levels(x@idents))))
need_lift <- any(!vapply(objects, function(x) identical(levels(x@idents), celltypes), logical(1)))
if (need_lift) {
  message("Lifting overview objects to a shared cell type set: ", paste(celltypes, collapse = ", "))
  objects <- lapply(objects, function(x) liftCellChat(x, group.new = celltypes))
  merged <- mergeCellChat(objects, add.names = conditions, cell.prefix = TRUE)
}

global_summary <- do.call(rbind, lapply(names(objects), function(condition) {
  obj <- objects[[condition]]
  data.frame(
    condition = condition,
    n_celltypes = length(levels(obj@idents)),
    total_count = sum(obj@net$count, na.rm = TRUE),
    total_weight = sum(obj@net$weight, na.rm = TRUE),
    nonzero_count_edges = sum(obj@net$count > 0, na.rm = TRUE),
    nonzero_weight_edges = sum(obj@net$weight > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.csv(global_summary, file.path(tables_dir, "cellchat_four_group_global_summary.csv"), row.names = FALSE)

network_long <- do.call(rbind, lapply(names(objects), function(condition) {
  obj <- objects[[condition]]
  rbind(
    mat_to_long(obj@net$count, condition, "count", celltypes),
    mat_to_long(obj@net$weight, condition, "weight", celltypes)
  )
}))
network_wide <- do.call(rbind, lapply(c("count", "weight"), function(measure) {
  tmp <- make_wide_difference(
    network_long[network_long$measure == measure, , drop = FALSE],
    keys = c("measure", "source", "target"),
    conditions = conditions,
    value_name = "value"
  )
  tmp
}))
write.csv(network_wide, file.path(tables_dir, "cellchat_four_group_network_difference_summary.csv"), row.names = FALSE)

pathway_long <- do.call(rbind, lapply(names(objects), function(condition) {
  array_to_long(objects[[condition]]@netP$prob, condition)
}))
pathway_wide <- make_wide_difference(
  pathway_long,
  keys = c("source", "target", "pathway_name"),
  conditions = conditions,
  value_name = "value"
)
write.csv(pathway_wide, file.path(tables_dir, "cellchat_four_group_pathway_difference_summary.csv"), row.names = FALSE)

write.csv(
  data.frame(
    table = c(
      "cellchat_four_group_global_summary",
      "cellchat_four_group_network_difference_summary",
      "cellchat_four_group_pathway_difference_summary"
    ),
    rows = c(nrow(global_summary), nrow(network_wide), nrow(pathway_wide)),
    stringsAsFactors = FALSE
  ),
  file.path(tables_dir, "cellchat_four_group_overview_table_summary.csv"),
  row.names = FALSE
)

bar_theme <- theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "none"
  )

save_gg_safe(
  compareInteractions(
    merged,
    group = seq_along(conditions),
    measure = "count",
    show.legend = FALSE,
    title.name = "Number of interactions"
  ) + bar_theme,
  file.path(outdir, "cellchat_overview_interaction_count_bar"),
  width * 0.7,
  height * 0.6,
  dpi
)

save_gg_safe(
  compareInteractions(
    merged,
    group = seq_along(conditions),
    measure = "weight",
    show.legend = FALSE,
    title.name = "Interaction strength"
  ) + bar_theme,
  file.path(outdir, "cellchat_overview_interaction_weight_bar"),
  width * 0.7,
  height * 0.6,
  dpi
)

network_heatmap <- function(measure) {
  df <- network_long[network_long$measure == measure, , drop = FALSE]
  df$source <- factor(df$source, levels = celltypes)
  df$target <- factor(df$target, levels = rev(celltypes))
  ggplot(df, aes(x = source, y = target, fill = value)) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_wrap(~ condition, nrow = 2) +
    scale_fill_gradient(low = "white", high = if (measure == "count") "#B2182B" else "#2166AC") +
    coord_equal() +
    labs(x = "Source", y = "Target", fill = measure, title = paste("CellChat", measure, "by source-target pair")) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      strip.text = element_text(size = 13, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16)
    )
}

save_gg_safe(
  network_heatmap("count"),
  file.path(outdir, "cellchat_overview_source_target_count_heatmap"),
  width,
  height,
  dpi
)

save_gg_safe(
  network_heatmap("weight"),
  file.path(outdir, "cellchat_overview_source_target_weight_heatmap"),
  width,
  height,
  dpi
)

save_gg_safe(
  rankNet(
    merged,
    mode = "comparison",
    stacked = TRUE,
    do.stat = FALSE,
    measure = "weight",
    title = "Overall signaling information flow"
  ) + theme_classic(base_size = 11),
  file.path(outdir, "cellchat_overview_pathway_information_flow_rank_stacked"),
  width,
  height * 0.75,
  dpi
)

count_edge_max <- max(unlist(lapply(objects, function(x) max(x@net$count, na.rm = TRUE))), na.rm = TRUE)
weight_edge_max <- max(unlist(lapply(objects, function(x) max(x@net$weight, na.rm = TRUE))), na.rm = TRUE)
pair_width <- max(width, 16)
pair_height <- max(height, 18)

save_base_safe(
  {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mfrow = c(2, 2), mar = c(1, 1, 4, 1), xpd = TRUE)
    for (condition in conditions) {
      obj <- objects[[condition]]
      group_size <- as.numeric(table(obj@idents))
      netVisual_circle(
        obj@net$count,
        title.name = paste0(condition, " - interaction count"),
        vertex.weight = group_size,
        weight.scale = TRUE,
        edge.weight.max = count_edge_max,
        edge.width.max = 8,
        vertex.size.max = 18,
        vertex.label.cex = 0.75,
        label.edge = FALSE
      )
    }
  },
  file.path(outdir, "cellchat_overview_network_count_circle_all_groups"),
  width,
  height,
  dpi
)

save_base_safe(
  {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mfrow = c(2, 2), mar = c(1, 1, 4, 1), xpd = TRUE)
    for (condition in conditions) {
      obj <- objects[[condition]]
      group_size <- as.numeric(table(obj@idents))
      netVisual_circle(
        obj@net$weight,
        title.name = paste0(condition, " - interaction strength"),
        vertex.weight = group_size,
        weight.scale = TRUE,
        edge.weight.max = weight_edge_max,
        edge.width.max = 8,
        vertex.size.max = 18,
        vertex.label.cex = 0.75,
        label.edge = FALSE
      )
    }
  },
  file.path(outdir, "cellchat_overview_network_weight_circle_all_groups"),
  width,
  height,
  dpi
)

pair_index <- combn(seq_along(conditions), 2, simplify = FALSE)

save_base_safe(
  {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mfrow = c(3, 2), mar = c(1, 1, 4, 1), xpd = TRUE)
    for (pair in pair_index) {
      a <- pair[[1]]
      b <- pair[[2]]
      netVisual_diffInteraction(
        merged,
        comparison = c(a, b),
        measure = "count",
        weight.scale = TRUE,
        vertex.size.max = 18,
        vertex.label.cex = 0.75,
        edge.width.max = 8,
        title.name = paste0(conditions[[b]], " vs ", conditions[[a]], " - count")
      )
    }
  },
  file.path(outdir, "cellchat_overview_pairwise_diff_count_circle"),
  pair_width,
  pair_height,
  dpi
)

save_base_safe(
  {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mfrow = c(3, 2), mar = c(1, 1, 4, 1), xpd = TRUE)
    for (pair in pair_index) {
      a <- pair[[1]]
      b <- pair[[2]]
      netVisual_diffInteraction(
        merged,
        comparison = c(a, b),
        measure = "weight",
        weight.scale = TRUE,
        vertex.size.max = 18,
        vertex.label.cex = 0.75,
        edge.width.max = 8,
        title.name = paste0(conditions[[b]], " vs ", conditions[[a]], " - weight")
      )
    }
  },
  file.path(outdir, "cellchat_overview_pairwise_diff_weight_circle"),
  pair_width,
  pair_height,
  dpi
)

plot_files <- list.files(outdir, pattern = "\\.(pdf|png)$", full.names = TRUE)
write.csv(
  data.frame(file = basename(plot_files), size_bytes = file.info(plot_files)$size, stringsAsFactors = FALSE),
  file.path(outdir, "cellchat_overview_plot_files.csv"),
  row.names = FALSE
)
write.csv(plot_log, file.path(outdir, "cellchat_overview_plot_log.csv"), row.names = FALSE)
message("Wrote overview tables to ", tables_dir)
message("Wrote ", length(plot_files), " overview plot files to ", outdir)
RS_OVERVIEW
fi

echo "Full CellChat pipeline finished."
