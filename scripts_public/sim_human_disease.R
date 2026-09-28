#!/usr/bin/env Rscript

# ==============================================================================
# 使用 SIMER 模拟猪群数据:
# - 总输出个体数: 1,000,000 (10 个世代, 每代 100,000)
# - 性状数: 7
# - 遗传力: 0.10 ~ 0.80
# - SNP 面板: 500k / 1M / 2M
#
# 关键重构:
# 1. 原脚本的 20 头/母本 + 0.8 母本选择率会按 8 倍/代爆炸，10 代根本不可控。
# 2. SIMER 内部用 big.matrix(type = "char") 存双行单倍型，不能指数扩群后再上 2M SNP。
# 3. 这里改为"每代固定 100k + 单世代写盘 + 写完立即瘦身"，这是唯一还算能落地的方案。
#
# 运行方式:
# Rscript 50m100n.R
#
# 可选参数:
# --markers=500000,1000000,2000000
# --pop-per-gen=100000
# --generations=10
# --ncpus=16
# --output-dir=./output/simer_1m_10gen_7trait
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("simer", quietly = TRUE)) {
    stop("未检测到 simer 包。先把环境配好，别在批处理作业里现装包。")
  }
  if (!requireNamespace("bigmemory", quietly = TRUE)) {
    stop("未检测到 bigmemory 包。这个脚本现在依赖 file-backed big.matrix。")
  }
  library(simer)
})

parse_args <- function(args) {
  parsed <- list()
  if (length(args) == 0L) {
    return(parsed)
  }

  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      stop("非法参数: ", arg, call. = FALSE)
    }
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    value <- kv[2]
    parsed[[key]] <- value
  }
  parsed
}

to_int_vec <- function(x) {
  as.integer(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
}

format_bytes <- function(bytes) {
  units <- c("B", "KB", "MB", "GB", "TB", "PB")
  value <- as.numeric(bytes)
  idx <- 1L
  while (value >= 1024 && idx < length(units)) {
    value <- value / 1024
    idx <- idx + 1L
  }
  sprintf("%.2f %s", value, units[idx])
}

toeplitz_cor <- function(n_trait, rho) {
  outer(
    seq_len(n_trait),
    seq_len(n_trait),
    FUN = function(i, j) rho ^ abs(i - j)
  )
}

build_named_list <- function(values, trait_names) {
  setNames(as.list(values), trait_names)
}

build_phe_models <- function(trait_names) {
  setNames(
    lapply(seq_along(trait_names), function(i) {
      sprintf("T%d = A + Sex + Year + Cov1 + E", i)
    }),
    trait_names
  )
}

build_pop_env <- function(trait_names) {
  sex_effect <- setNames(
    lapply(seq_along(trait_names), function(i) c(0.40 + 0.03 * i, -0.40 - 0.03 * i)),
    trait_names
  )
  year_effect <- setNames(
    lapply(seq_along(trait_names), function(i) c(0, 0.50 + 0.10 * i, 1.10 + 0.15 * i)),
    trait_names
  )
  cov1_slope <- setNames(
    as.list(seq(0.30, 0.90, length.out = length(trait_names))),
    trait_names
  )

  list(
    Sex = list(
      level = c("Male", "Female"),
      effect = sex_effect
    ),
    Year = list(
      level = c("2022", "2023", "2024"),
      effect = year_effect
    ),
    Cov1 = list(
      level = c(10, 20, 30),
      slope = cov1_slope
    )
  )
}

build_qtn_num <- function(marker_n, trait_names) {
  base_counts <- c(120L, 140L, 160L, 180L, 200L, 220L, 240L)
  scale_factor <- marker_n / 500000
  qtn_counts <- pmax(50L, as.integer(round(base_counts * scale_factor)))
  setNames(as.list(qtn_counts), trait_names)
}

append_log <- function(log_file, ...) {
  msg <- paste0(format(Sys.time(), "%F %T"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

trim_to_latest_generation <- function(SP) {
  latest_geno <- SP$geno$pop.geno[[length(SP$geno$pop.geno)]]
  latest_pop <- SP$pheno$pop[[length(SP$pheno$pop)]]

  SP$geno$pop.geno <- list(latest_geno)
  names(SP$geno$pop.geno) <- "gen1"

  SP$pheno$pop <- list(latest_pop)
  names(SP$pheno$pop) <- "gen1"

  invisible(gc())
  SP
}

build_base_geno_fb <- function(marker_n, pop_ind, backing_dir, file_prefix, seed) {
  nrow_fb <- 2L * pop_ind
  desc_file <- paste0(file_prefix, ".desc")
  bin_file <- paste0(file_prefix, ".bin")

  if (!dir.exists(backing_dir)) {
    dir.create(backing_dir, recursive = TRUE, showWarnings = FALSE)
  }

  unlink(file.path(backing_dir, desc_file))
  unlink(file.path(backing_dir, bin_file))

  bm <- bigmemory::filebacked.big.matrix(
    nrow = nrow_fb,
    ncol = marker_n,
    type = "char",
    init = 0,
    backingpath = backing_dir,
    backingfile = bin_file,
    descriptorfile = desc_file
  )

  # 控制单次随机块大小，避开 R 的 32 位整数 sample 上限。
  target_elements <- 2e8
  block_rows <- max(32L, as.integer(floor(target_elements / marker_n)))
  block_rows <- min(block_rows, nrow_fb)

  set.seed(seed)
  for (row_start in seq.int(1L, nrow_fb, by = block_rows)) {
    row_end <- min(row_start + block_rows - 1L, nrow_fb)
    n_block_rows <- row_end - row_start + 1L
    n_values <- as.double(n_block_rows) * as.double(marker_n)
    block <- matrix(
      sample.int(2L, size = n_values, replace = TRUE) - 1L,
      nrow = n_block_rows,
      ncol = marker_n
    )
    bm[row_start:row_end, ] <- block
    rm(block)
    invisible(gc())
  }

  bm
}

format_marker_label <- function(marker_n) {
  if (marker_n >= 1000000L && marker_n %% 1000000L == 0L) {
    return(paste0(marker_n %/% 1000000L, "m"))
  }
  if (marker_n >= 1000L && marker_n %% 1000L == 0L) {
    return(paste0(marker_n %/% 1000L, "k"))
  }
  as.character(marker_n)
}

run_panel <- function(marker_n, cfg, trait_names, phe_h2a, phe_var, phe_models, pop_env, cor_a, cor_e) {
  panel_label <- format_marker_label(marker_n)
  panel_dir <- file.path(cfg$base_output_directory, paste0("panel_", panel_label))
  if (!dir.exists(panel_dir)) {
    dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
  }

  log_file <- file.path(panel_dir, paste0("simulation_", panel_label, ".log"))
  qtn_num <- build_qtn_num(marker_n, trait_names)
  panel_seed <- cfg$seed_base + match(marker_n, cfg$marker_panels) * 1000L
  base_geno_dir <- file.path(panel_dir, "base_geno")

  current_gen_bigmat_bytes <- 2 * cfg$pop_per_gen * marker_n
  write_stage_bytes <- 1 * cfg$pop_per_gen * marker_n
  conservative_peak_bytes <- 4 * cfg$pop_per_gen * marker_n

  append_log(log_file, "开始面板 ", panel_label,
             " | markers=", marker_n,
             " | 每代个体=", cfg$pop_per_gen,
             " | 世代数=", cfg$generations,
             " | 当前代双行基因型估计=", format_bytes(current_gen_bigmat_bytes),
             " | 写盘阶段估计=", format_bytes(write_stage_bytes),
             " | 保守峰值估计=", format_bytes(conservative_peak_bytes))

  append_log(log_file, "分块构建基础 file-backed 基因型开始")
  base_geno <- build_base_geno_fb(
    marker_n = marker_n,
    pop_ind = cfg$pop_per_gen,
    backing_dir = base_geno_dir,
    file_prefix = paste0("basegeno_", panel_label),
    seed = panel_seed
  )
  append_log(log_file, "基础 file-backed 基因型完成")

  SP <- param.simer(
    out = paste0("pig_1m_7trait_", panel_label, "_gen01"),
    outpath = panel_dir,
    out.format = "plink",
    pop.ind = cfg$pop_per_gen,
    pop.gen = 2,
    prog = cfg$prog,
    ps = cfg$ps,
    reprod.way = "randmate",
    num.chr = 18,
    len.chr = 1.5e8,
    pop.geno = base_geno,
    inrows = 2,
    pop.marker = marker_n,
    index.wt = rep(1 / cfg$n_traits, cfg$n_traits),
    pop.env = pop_env,
    qtn.model = "A",
    qtn.num = qtn_num,
    phe.model = phe_models,
    phe.h2A = phe_h2a,
    phe.var = phe_var,
    phe.corA = cor_a,
    phe.corE = cor_e,
    rate.mut = NULL,
    out.geno.gen = 1,
    out.pheno.gen = 1,
    seed.sim = panel_seed,
    ncpus = cfg$ncpus,
    verbose = TRUE
  )

  append_log(log_file, "注释与基础群体初始化开始")
  SP <- annotation(SP = SP, verbose = TRUE)
  SP <- genotype(SP = SP, ncpus = cfg$ncpus, verbose = TRUE)
  SP <- phenotype(SP = SP, ncpus = cfg$ncpus, verbose = TRUE)

  total_written <- 0L

  for (gen_idx in seq_len(cfg$generations)) {
    current_pop_n <- nrow(SP$pheno$pop[[1]])
    total_written <- total_written + current_pop_n

    SP$global$replication <- gen_idx
    SP$global$out <- sprintf("pig_1m_7trait_%s_gen%02d", panel_label, gen_idx)
    SP$global$outpath <- panel_dir
    SP$global$out.geno.gen <- 1
    SP$global$out.pheno.gen <- 1

    append_log(log_file, "写出第 ", gen_idx, " 代 | 个体数=", current_pop_n)
    SP <- write.file(SP)

    if (gen_idx < cfg$generations) {
      append_log(log_file, "生成第 ", gen_idx + 1L, " 代")
      SP <- simer:::mate.randmate(SP = SP, ncpus = cfg$ncpus, verbose = TRUE)
      SP <- trim_to_latest_generation(SP)
    }
  }

  meta_lines <- c(
    "==================================================",
    paste0("面板: ", panel_label),
    paste0("创建时间: ", Sys.time()),
    paste0("总输出个体数: ", total_written),
    paste0("每代个体数: ", cfg$pop_per_gen),
    paste0("总世代数: ", cfg$generations),
    paste0("SNP 数: ", marker_n),
    paste0("线程数: ", cfg$ncpus),
    paste0("遗传力: ", paste(unlist(phe_h2a), collapse = ", ")),
    paste0("方差: ", paste(unlist(phe_var), collapse = ", ")),
    paste0("QTN 数: ", paste(unlist(qtn_num), collapse = ", ")),
    paste0("输出目录: ", panel_dir),
    "=================================================="
  )
  writeLines(meta_lines, con = file.path(panel_dir, paste0("panel_", panel_label, "_summary.txt")))

  append_log(log_file, "面板 ", panel_label, " 完成 | 总输出个体数=", total_written)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

cfg <- list(
  base_output_directory = if (!is.null(args[["output-dir"]])) {
    args[["output-dir"]]
  } else {
    "./output/simer_1m_10gen_7trait"
  },
  marker_panels = if (!is.null(args[["markers"]])) {
    to_int_vec(args[["markers"]])
  } else {
    c(500000L, 1000000L, 2000000L)
  },
  pop_per_gen = if (!is.null(args[["pop-per-gen"]])) as.integer(args[["pop-per-gen"]]) else 100000L,
  generations = if (!is.null(args[["generations"]])) as.integer(args[["generations"]]) else 10L,
  ncpus = if (!is.null(args[["ncpus"]])) as.integer(args[["ncpus"]]) else 16L,
  prog = 2L,
  ps = c(1, 1),
  n_traits = 7L,
  seed_base = 20260623L
)

if (!dir.exists(cfg$base_output_directory)) {
  dir.create(cfg$base_output_directory, recursive = TRUE, showWarnings = FALSE)
}

trait_names <- paste0("tr", seq_len(cfg$n_traits))
phe_h2a <- build_named_list(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.65, 0.80), trait_names)
phe_var <- build_named_list(c(8, 9, 10, 11, 12, 13, 14), trait_names)
phe_models <- build_phe_models(trait_names)
pop_env <- build_pop_env(trait_names)
cor_a <- toeplitz_cor(cfg$n_traits, rho = 0.60)
cor_e <- toeplitz_cor(cfg$n_traits, rho = 0.25)

cat("==================================================\n")
cat("SIMER 大规模模拟任务启动\n")
cat("输出目录:", cfg$base_output_directory, "\n")
cat("每代个体数:", cfg$pop_per_gen, "\n")
cat("世代数:", cfg$generations, "\n")
cat("总输出个体数:", cfg$pop_per_gen * cfg$generations, "\n")
cat("性状数:", cfg$n_traits, "\n")
cat("遗传力:", paste(unlist(phe_h2a), collapse = ", "), "\n")
cat("SNP 面板:", paste(cfg$marker_panels, collapse = ", "), "\n")
cat("线程数:", cfg$ncpus, "\n")
cat("说明: 为压住 2M SNP 的内存峰值，这里固定每代 100k 个体，写完当前代立刻瘦身。\n")
cat("==================================================\n")

for (marker_n in cfg$marker_panels) {
  run_panel(
    marker_n = marker_n,
    cfg = cfg,
    trait_names = trait_names,
    phe_h2a = phe_h2a,
    phe_var = phe_var,
    phe_models = phe_models,
    pop_env = pop_env,
    cor_a = cor_a,
    cor_e = cor_e
  )
}

cat("全部面板任务完成。\n")
