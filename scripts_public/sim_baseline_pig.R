# ==============================================================================
# 使用 SIMER 模拟猪生长性状数据集 (纯净版: 仅 A + E)
#
# 目的:
# 1. 创建一个超大规模 (1.28M) 的数据集。
# 2. **完全去除** 固定效应、协变量和额外随机效应。
# 3. 表型构成: P = TBV + E
#
# 基本要求:
# - 群体规模: ~1,279,976 头 (G1=2625 ... G6=880804)
# - 世代数: 6 代
# - SNP位点: 50,000 个
# - 模型: 仅包含加性效应和残差
# ==============================================================================

# 步骤 0: 环境准备
if (!require("simer", quietly = TRUE)) {
  if (!require("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("xiaolei-lab/SIMER")
}
library(simer)

# ==============================================================================
# 步骤 1: 定义公共参数和输出路径
# ==============================================================================

# 设置输出路径
base_output_directory <- "./output/128w" 
if (!dir.exists(base_output_directory)) {
  dir.create(base_output_directory, recursive = TRUE)
}
output_dir_A <- file.path(base_output_directory)
if (!dir.exists(output_dir_A)) dir.create(output_dir_A)

# --- 群体结构参数 (目标 ~1280k, 6世代) ---
# 保持之前的结构以精确控制数量
p_ind <- 2625  # 创始群体大小
p_gen <- 6     # 世代数
p_prog <- 8    # 每个母本的后代数

# ==============================================================================
# 步骤 2: 模拟数据集 - 纯净加性模型 (P = A + E)
# ==============================================================================
cat(">>> 开始模拟数据集: 128万, 50k SNP, 无固定效应...\n")

# 定义遗传参数
h2A_multi <- list(tr1 = 0.5, tr2 = 0.4, tr3 = 0.3)
var_multi <- list(tr1 = 10, tr2 = 8, tr3 = 12) 

# 加性遗传相关
corA_multi <- matrix(c(1.0, 0.4, 0.2,
                       0.4, 1.0, 0.5,
                       0.2, 0.5, 1.0), 3, 3) 

# 残差相关
corE_multi <- matrix(c(1.0, 0.3, 0.1,
                       0.3, 1.0, 0.2,
                       0.1, 0.2, 1.0), 3, 3)

SP <- param.simer(
  # --- 输出设置 ---
  out = "data_A_1280k_simple",
  outpath = output_dir_A,
  out.format = "plink",
  
  # --- 群体和基因组参数 ---
  pop.ind = p_ind,
  pop.gen = p_gen,
  prog = p_prog,
  ps = c(0.5, 0.8),
  reprod.way = "randmate",
  species = "pig",
  pop.marker = 50000, # 5万 SNP
  
  # --- 选择权重 ---
  index.wt = c(0.4, 0.3, 0.3), 
  
  # --- 关键修改: 去除 pop.env ---
  pop.env = NULL, 
  
  qtn.model = "A",
  qtn.num = list(tr1 = 100, tr2 = 120, tr3 = 110),
  
  # --- 关键修改: 模型只包含 A 和 E ---
  phe.model = list(
    tr1 = "T1 = A + E",
    tr2 = "T2 = A + E",
    tr3 = "T3 = A + E"
  ),
  
  phe.h2A = h2A_multi,
  phe.var = var_multi,
  phe.corA = corA_multi,
  phe.corE = corE_multi,
  
  # --- 其他设置 ---
  seed.sim = 1280088, 
  ncpus = 8 
)

# 执行模拟
simer_results <- simer(SP)

# 打印信息
total_inds <- sum(sapply(simer_results$pheno$pop, nrow))
cat(">>> 完成: 总个体数:", total_inds, "\n")
cat(">>> 文件已输出至:", output_dir_A, "\n")
