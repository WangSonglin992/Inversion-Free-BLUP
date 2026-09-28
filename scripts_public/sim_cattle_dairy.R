# ==============================================================================
# SIMER 模拟全流程: 500万系谱 (50万/代 x 10代) -> 表型掩盖 -> PLINK提取子集
# ==============================================================================

library(simer)
library(data.table)

# --- 0. 基础设置 ---
# 输出路径 (修改文件夹名称以区分)
base_output_directory <- "./output/cattle_5m"
if (!dir.exists(base_output_directory)) dir.create(base_output_directory, recursive = TRUE)

# 基础文件名
file_prefix <- "cattle_5M_pedigree"

# --- 1. 模拟参数修改 ---
# 目标: 总共 500 万个体
# 逻辑: randmate 模式下群体规模恒定，10代 x 50万 = 500万
p_ind <- 500000   # 从 100万 降至 50万
p_gen <- 10       # 保持 10 代
p_prog <- 2       # 保持扩繁比例

# --- 2. 执行 SIMER 模拟 ---
cat(paste0("\n>>> [Step 1] 开始模拟 ", p_ind * p_gen / 10000, " 万个体 (50万/代 x 10代)...\n"))

SP <- param.simer(
  out = file_prefix,
  outpath = base_output_directory,
  out.format = "plink", 
  
  pop.ind = p_ind,
  pop.gen = p_gen,
  prog = p_prog,
  reprod.way = "randmate",
  species = "cattle", 
  
  pop.env = NULL,
  qtn.model = "A",
  phe.model = list(tr1 = "T1 = A + E"),
  phe.h2A = list(tr1 = 0.35),
  phe.var = list(tr1 = 100),
  
  seed.sim = 20241010,
  ncpus = 32
)

# 运行模拟
# 注意: SIMER 的 simer() 函数主要进行表型模拟。
# 如果报错提示 G 矩阵内存不足，说明 SIMER 尝试在内存中进行某种关联分析，
# 但对于纯粹的数据生成(Data Generation)，通常不会构建完整 G 矩阵。
simer_results <- simer(SP)
cat(">>> 模拟完成，完整文件已生成。\n")

# --- 3. 表型掩盖 (前 250 万设为 NA) ---
# 阈值调整为总群体的一半
mask_threshold <- 2500000 

cat(paste0("\n>>> [Step 2] 处理表型缺失 (掩盖前 ", mask_threshold/10000, " 万)...\n"))
phe_file <- file.path(base_output_directory, paste0(file_prefix, ".phe"))

if (file.exists(phe_file)) {
  dt_pheno <- fread(phe_file, header = TRUE)
  total_rows <- nrow(dt_pheno)
  
  cat(paste(">>> 当前总行数:", total_rows, "\n"))
  
  if (total_rows >= mask_threshold) {
    target_col <- grep("^T1$|^tr1$", names(dt_pheno), value = TRUE)[1]
    
    if (!is.na(target_col)) {
      # 原地修改前 250 万行
      dt_pheno[1:mask_threshold, (target_col) := NA]
      
      fwrite(dt_pheno, file = phe_file, sep = "\t", quote = FALSE, na = "NA")
      cat(paste(">>> 表型文件更新成功: 前", mask_threshold, "个个体表型已设为 NA。\n"))
    } else {
      cat(">>> 警告: 未找到表型列，跳过掩盖。\n")
    }
  } else {
    cat(">>> 警告: 生成的个体数少于掩盖阈值，未执行掩盖。\n")
  }
} else {
  cat(">>> 错误: 未找到生成的表型文件 (.phe)。\n")
}

# --- 4. 自动调用 PLINK 提取基因型 (最后 160 万) ---
cat("\n>>> [Step 3] 正在调用 PLINK 提取最后 160 万个体的基因型...\n")

fam_file <- file.path(base_output_directory, paste0(file_prefix, ".fam"))
keep_file <- file.path(base_output_directory, "keep_list.txt")
subset_prefix <- "cattle_1.6M_subset"

if (file.exists(fam_file)) {
  dt_fam <- fread(fam_file, header = FALSE, select = c(1, 2)) 
  n_total <- nrow(dt_fam)
  # 依然提取最后 160 万用于验证/计算，这个数量级在 2TB 内存上计算 G 矩阵是安全的
  # 160万^2 * 8 bytes ≈ 20 TB (依然很大，如果是单精度或者通过分块计算是可行的)
  # 但如果你只是提取数据，这步没问题。
  n_subset <- 1600000 
  
  if (n_total >= n_subset) {
    start_idx <- n_total - n_subset + 1
    dt_keep <- dt_fam[start_idx:n_total]
    
    fwrite(dt_keep, file = keep_file, sep = "\t", col.names = FALSE)
    cat(paste(">>> 已生成提取列表:", keep_file, "(包含", nrow(dt_keep), "个个体)\n"))
    
    input_bfile <- file.path(base_output_directory, file_prefix)
    output_bfile <- file.path(base_output_directory, subset_prefix)
    
    plink_cmd <- paste(
      "plink",
      "--bfile", input_bfile,
      "--keep", keep_file,
      "--make-bed",
      "--out", output_bfile,
      "--threads 32",
      "--silent"
    )
    
    cat(paste(">>> 正在执行系统命令:", plink_cmd, "\n"))
    exit_code <- system(plink_cmd)
    
    if (exit_code == 0) {
      cat("\n=================================================================\n")
      cat(">>> [SUCCESS] 任务全部完成！\n")
      cat(">>> 1. 完整数据集 (500万): ", input_bfile, "\n")
      cat(">>> 2. 基因型子集 (160万): ", output_bfile, "\n")
      cat(">>> 3. 表型文件: 已更新 (前250万NA)\n")
      cat("=================================================================\n")
    } else {
      cat(">>> [ERROR] PLINK 执行失败，退出代码:", exit_code, "\n")
    }
  } else {
    cat(paste(">>> 错误: 总个体数 (", n_total, ") 少于需提取的子集数 (", n_subset, ")。\n"))
  }
} else {
  cat(">>> 错误: 未找到 .fam 文件。\n")
}
