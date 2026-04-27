#### 全部特征####
#### 两种综合评价指标的计算####
library(data.table)
library(stringr)

### ========== 1. 读取文件 ==========
external <- fread("E:/Protocol_dataset/TCGAbiolinks/TCGA_ALL/SNF/BRCA_SNF/SNFresults0817/Omics_Combination_Evaluation_SNF.csv")
internal <- fread("E:/Protocol_dataset/TCGAbiolinks/TCGA_ALL/SNF/BRCA_SNF/SNF_internal_metrics.csv")
purity   <- fread("E:/Protocol_dataset/TCGAbiolinks/TCGA_ALL/SNF/BRCA_SNF/SNFresults0817/Original_Purity_Gini_Entropy/Purity_Gini_Entropy_SNF_Summary.csv")
clinical <- fread("E:/Protocol_dataset/TCGAbiolinks/TCGA_ALL/SNF/BRCA_SNF/SNFresults0817/clinical_ECP/clinical_enrichment_SNF_summary.csv")
km       <- fread("E:/Protocol_dataset/TCGAbiolinks/TCGA_ALL/SNF/BRCA_SNF/SNFresults0817/KM_Plots/SNF_KM_summary_results.csv")

### ========== 2. 组合名称标准化函数 ==========
normalize_name <- function(x) {
   x <- gsub("SNF_", "", x)           # purity 文件前缀
   x <- gsub("\\s+", "", x)           # 删除空格
   x <- toupper(x) ## 全变大写，方便处理
   
   # 先把所有分隔符统一成 "_"
   x <- gsub("[+\\-\\.]", "_", x)
   
   # ——把 METHYLATION 统一成 MET——,因为前面已经将其都改为大写了
   x <- gsub("METHYLATION", "MET", x)
   
   # 最终统一为 "+"
   x <- gsub("_", "+", x)
   
   return(x)
}

# 找到组合列名
find_comb_col <- function(df) {#定义了一个名为 find_comb_col 的函数，它接受一个参数 df，通常是一个数据框（data frame
   cand <- c("Combination","Omics_Combination","Method","File","name") #字符向量 cand（candidate 的缩写），里面包含我们“候选”的列名
   existing <- intersect(cand, names(df))#intersect(cand, names(df)) 找出 cand 和 df 的列名之间的交集
   if (length(existing) > 0) return(existing[1]) #结果保存在变量 existing 中，它是一个字符向量，包含同时出现在 cand 和 df 列名中的名字，并且保留 cand 中的顺序（因为 intersect 在 R 中对第一个参数的顺序敏感）
   return(names(df)[1])#检查 existing 是否非空（即至少有一个候选列名存在于 df 中）。如果存在，则返回 existing 的第一个元素（也就是 cand 中优先级最高且确实存在的那个列名）
}

ext_col <- find_comb_col(external)
int_col <- find_comb_col(internal)
pur_col <- find_comb_col(purity)
cli_col <- find_comb_col(clinical)
km_col  <- find_comb_col(km)

# 应用组合标准化
external[, Combination_norm := normalize_name(get(ext_col))]#:= 是 data.table 特有的赋值操作符，用于按引用修改或添加列
internal[, Combination_norm := normalize_name(get(int_col))]
purity[,   Combination_norm := normalize_name(get(pur_col))]
clinical[, Combination_norm := normalize_name(get(cli_col))]
km[,       Combination_norm := normalize_name(get(km_col))]


### ========== 3. 构建 master 表（仅用于标准化） ==========
master <- data.table(Combination_norm = external$Combination_norm)


### ========== 4. 按顺序合并指标 ==========
# 外部指标
cols_ext <- intersect(c("NMI","ARI","Accuracy","MacroF1"),
                      names(external))
master <- merge(master, external[, c("Combination_norm", cols_ext), with = FALSE], 
                by="Combination_norm", all.x=TRUE,sort=FALSE)
#with = FALSE 表示：把 c(...) 当作列名（字符串）来处理，而不是当作变量或表达式。
#动态选列（cols_cli 可以是任意列名组合）
#all.x = TRUE：左连接，保留 master 所有行
#sort = FALSE：关键！ 告诉 merge() 不要对结果按键（Combination_norm）排序，而是保持 x（即当前 master）的原始行顺序。

# 内部指标（无 Dunn）
cols_int <- intersect(c("Silhouette","CH","DB"), names(internal))
master <- merge(master, internal[, c("Combination_norm", cols_int), with = FALSE], 
                by="Combination_norm", all.x=TRUE,sort=FALSE)

# purity（注意：已去 SNF_）
cols_purity <- intersect(c("Avg_Purity","Avg_Gini","Avg_Entropy"), names(purity))
master <- merge(master, purity[, c("Combination_norm", cols_purity), with = FALSE],
                by="Combination_norm", all.x=TRUE,sort=FALSE)

# 临床 ECP
if ("ECP" %in% names(clinical)) {
   clinical[, "ECP_num" := ECP]
} else if ("Total_Significant" %in% names(clinical)) {
   clinical[, "ECP_num" := Total_Significant]
}
cols_cli <- intersect(c("ECP_num"), names(clinical))
master <- merge(master, clinical[, c("Combination_norm",cols_cli), with=FALSE], 
                by="Combination_norm", all.x=TRUE,sort=FALSE)

# KM log-rank
cols_km <- intersect(c("log_rank_p","p_value","p.value","logrank_p"), names(km))
if (length(cols_km) > 0) {
   km[, LRT_p := get(cols_km[1])]
}
master <- merge(master, km[, c("Combination_norm",cols_km), with=FALSE],
                by="Combination_norm", all.x=TRUE,sort=FALSE)

#关于方向的问题：
#外部指标：
# 外部指标（NMI, ARI, Accuracy, MacroF1）
# 都是 越大越好 → 
#纯度类指标
# ✔ Purity / Entropy / Gini
# 方向：
# 指标	           越大越好?
# Avg_Purity	     ✔ 越大越好
# Avg_Gini	        ✘ 越小越好
# Avg_Entropy	     ✘ 越小越好

#log-rank p 越小代表生存差异越显著 → 更好
#但 p 值是 “越小越好”，其他指标（NMI、Purity、Silhouette 等）是 “越大越好”，直接标准化 p 值会导致方向相反，最终影响综合分数
#故log-rank p 是否需要先取 -log10

# 内部指标
# 注意方向：
# 指标	          越大越好？
# Silhouette	✔   越大越好
# CH	✔            越大越好
# DB	✘            越小越好


### ----- 1. 先处理方向性 -----
#对于那些越小越好的指标（例如 Avg_Gini、Avg_Entropy、Davies–Bouldin），
#不要直接把原始数值取负（-x），而应该先做 min–max 归一化（得到 0..1），
#然后用 1 - normalized（等价把最小值映射为 1，最大值映射为 0）
#为什么**不要直接把原始值取负（-x）**来反转方向？

#取负把原始分布的尺度保留了（如果不同指标数值量级差别大，负值尺度仍然不一致），之后如果直接做 min–max 还是能归一化，但在某些实现里直接把负数带入之后的解释容易混淆。
#更稳健、可解释的方法是：先 min–max 再 invert（1 - normalized）。这保证了所有指标按相同规则被压到 [0,1]，且方向一致。
#用 1 - normalized 更清晰地说明“原来越小越好，现在越大越好”。

#去除之前不对的预处理
master <- master[,-c("DB_rev","Avg_Gini_rev","Avg_Entropy_rev")]
master <- master[,-"ECP_score"]
# 1. 生存 p 值 → 转为越大越好，且避免 log(0)
master[, log_rank_score := -log10(pmax(log_rank_p, 1e-300))]  # 更安全 # 转为 -log10(p)，越大越好
#防止 p=0 → pmax(...,1e-300) 合理

# 2. 内部指标方向修正（全部转为 [0,1] 或有界正数）
master[, DB_rev_score := 1 / (1 + DB)]  # DB ∈ [0,∞) → score ∈ (0,1]

# 3. 外部指标中的纯度类指标方向修正（转为“纯度”风格，越大越好）
master[, Avg_Gini_score := 1 - Avg_Gini]               # Gini ∈ [0,1] → score ∈ [0,1]
master[, Avg_Entropy_score := 1 - Avg_Entropy / log(5)]  # Entropy ∈ [0, log(5)] → score ∈ [0,1]#
#保存
# master包含原始数据以及修改方向的数据
save(master,file = "./SNF_pre_combi.Rdata")
# 保存为CSV文件
write.csv(master, file = "SNF_orig.csv", row.names = FALSE)  # row.names = FALSE 表示不保存行名称
#由于 Silhouette 包含负值，必须确保 Min-Max 能处理全范围：
normalize_safe <- function(x) {
   x <- as.numeric(x)
   rng <- range(x, na.rm = TRUE)
   if (diff(rng) < .Machine$double.eps) return(rep(1, length(x)))#
   (x - rng[1]) / diff(rng) #经典 Min-Max 标准化,经典 Min-Max 标准化,输出区间为 0～1。
}
#diff(c(min, max)) = max - min  也就是向量的 数据范围。
#.Machine$double.eps = 2.22e−16是 R 语言中 双精度浮点能够区分的最小差值。相当于“计算机认为的零”。

# 应用于 Silhouette
# master[, Silhouette_norm := normalize_safe(Silhouette)]

#### ---------2.进行最大最小归一化--------------------------

tmp <- colnames(master)
# 4. 定义所有用于综合评分的指标（全部是“越大越好”的原始值）
metric_cols <- c(
   "NMI", "ARI", "Accuracy", "MacroF1", 
   "Avg_Purity","Avg_Gini_score", "Avg_Entropy_score",
   "Silhouette", "CH", "DB_rev_score",
   "log_rank_score", "ECP_num"   # ← 注意：这里是 ECP_num，不是 ECP_score
)
#提取修改方向一致的数据
snf_combi <- master[,c("Combination_norm",..metric_cols),with = FALSE]

# 5. 统一 Min-Max 标准化（只做这一次！）
normalize_safe <- function(x) {
   x <- as.numeric(x)
   rng <- range(x, na.rm = TRUE)
   if (diff(rng) < .Machine$double.eps) return(rep(1, length(x)))#如果所有值都一样，那么它们“理论上已经处于同一水平”，归一化后让它们都等于 1 是最保险的做法
   (x - rng[1]) / diff(rng)#经典 Min-Max 标准化,经典 Min-Max 标准化,输出区间为 0～1。
}

# 对每列标准化
snf_combi_norm <- copy(snf_combi)
for (col in metric_cols) {  # ← 注意：这里只循环 metric_cols
   snf_combi_norm[[col]] <- normalize_safe(snf_combi_norm[[col]])
}

#保存
save(master,snf_combi,snf_combi_norm,file = "./SNF_pre_combi.Rdata")
#  (min-max)
#master[, ECP_score := (ECP_num - min(ECP_num)) / (max(ECP_num) - min(ECP_num))]


####------根据参考文献的AWA原理计算综合评价指标------------------
#由于我们更偏向于乳腺癌因此去除了内部指标
compute_AWA_sum <- function(df,
                            external_norm_cols,   # 已经归一化好的外部指标列名
                            clinical_norm_cols,   # 已经归一化好的临床指标列名
                            w1 = 0.5,             # 外部组权重
                            w2 = 0.5) {           # 临床组权重
   
   dt <- copy(as.data.table(df))
   
   # ---- 1) 外部指标求和 ----
   if (length(external_norm_cols) == 0) {
      dt[, S_sum := 0]
      m1 <- 0
   } else {
      dt[, S_sum := rowSums(.SD, na.rm = TRUE), .SDcols = external_norm_cols]
      m1 <- length(external_norm_cols)
   }
   
   # ---- 2) 临床指标求和 ----
   if (length(clinical_norm_cols) == 0) {
      dt[, L_sum := 0]
      m2 <- 0
   } else {
      dt[, L_sum := rowSums(.SD, na.rm = TRUE), .SDcols = clinical_norm_cols]
      m2 <- length(clinical_norm_cols)
   }
   
   # ---- 3) AWA（使用 sum 而不是 mean）----
   denom <- m1 * w1 + m2 * w2
   if (denom <= 0) stop("Denominator cannot be zero.")
   
   dt[, AWA_sum := (S_sum * w1 + L_sum * w2) / denom]
   
   return(dt)
}

external_cols <- c("NMI","ARI","Accuracy","MacroF1",
                   "Avg_Purity","Avg_Gini_score","Avg_Entropy_score")

clinical_cols <- c("log_rank_score","ECP_num")

snf_res <- compute_AWA_sum(df = snf_combi_norm,
                       external_norm_cols = external_cols,
                       clinical_norm_cols = clinical_cols,
                       w1 = 0.5,
                       w2 = 0.5)

#保存基于AWA的综合得分
save(snf_combi_norm,snf_res,file="./SNF_AWA_combi.Rdata")

#综合评分1的画图
library(ggplot2)
library(data.table)

# 提取用于绘图的数据
plot_data <- snf_res[
   , .(Combination_norm, AWA_sum)#.() 更高效，且是 data.table 推荐写法，取出这两列
][order(-AWA_sum)]  # 降序：高分在上

# 可选：只显示 top N（例如 top 15），避免标签拥挤
# plot_data <- head(plot_data, 15)

# 创建组合名称的因子（确保排序正确）
plot_data[, Combination_norm := factor(Combination_norm, levels = rev(Combination_norm))]

# 绘图
p1 <- ggplot(plot_data, aes(x = Combination_norm, y = AWA_sum)) +
   geom_col(fill = "steelblue", width = 0.7) +
   coord_flip() +  # 横向显示，便于阅读长名称
   labs(
      title = "Comprehensive Performance of Multi-omics Combinations (SNF)",
      #subtitle = "Higher AWA_sum indicates better performance",
      x = "Multi-omics Combination",
      y = "AWA Score (Range: ~0–1)"
   ) +
   theme_minimal(base_size = 14) +#白底 + 无边框,只有水平浅灰线（Y 轴方向）,突出数据本身
   theme(
      plot.title = element_text(size =13.5,hjust = 0.5, face = "bold"),#hjust = 0.5：水平居中对齐,face = "bold"：加粗字体
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),   # 统一设置 x 和 y 轴标题字体, 坐标轴标题
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
      
      # 关键：实现只有左+下边框
      # panel.border = element_blank(),      # 移除默认边框
      # panel.grid = element_blank(),        # 移除所有网格线（可选）
      # 
      # axis.line = element_line(color = "black", linewidth = 0.5),  # 添加 L 型轴线
      # axis.ticks = element_line(color = "black"),                  # 显示刻度线
      # axis.ticks.length = unit(0.15, "cm")                         # 刻度长度
   )
# 显示图形
print(p1)

# 保存
ggsave(plot = p1,
       file = "SNF_AWA_Score_BarPlot.png",  
       width = 9, 
       height = 7, 
       units = "in",
       dpi = 300)

#
#保存为PDF或svg方便缩小放大而不失真。导出为 PDF（矢量图，无限缩放不失真）
#为了观感一致，将宽度调整为一半6，高度保持不变，老师推荐还是一样就好，少放几个算法，26个组合的一个代表，51的放一个。
ggsave("SNF_AWA_Score_BarPlot.pdf", plot = p1, width = 9, height = 7, units = "in")



external_cols <- c("NMI","ARI","Accuracy","MacroF1",
                   "Avg_Purity","Avg_Gini_score","Avg_Entropy_score")
internal_cols <- c("Silhouette", "CH", "DB_rev_score")
clinical_cols <- c("log_rank_score","ECP_num")

#S_sum := rowSums(.SD, na.rm = TRUE), .SDcols = external_norm_cols] 外部
# ---- 计算 AWA_external ----
snf_res[, AWA_external := rowMeans(.SD, na.rm = TRUE), .SDcols = external_cols]

# ---- 计算 AWA_internal ----不建议加入画图，会干扰主线,但可以自己私下看看
#snf_res[, AWA_internal := rowMeans(.SD, na.rm = TRUE), .SDcols = internal_cols]

# ---- 计算 AWA_clinical ----
snf_res[, AWA_clinical := rowMeans(.SD, na.rm = TRUE), .SDcols = clinical_cols]
## 2D图
library(ggplot2)
library(RColorBrewer)
distinct_colors <- c(
   "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
   "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#A6761D", "#666666",
   "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69",
   "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
)
p <- ggplot(snf_res, aes(x = S_sum, y = L_sum,#这里用S_sum和L_sum试试，自己计算的好像不大行
                    color = Combination_norm)) +
   geom_point(size = 5, alpha = 0.9) +
   theme_classic(base_size = 14) +# # base_size 控制默认字体，但会被下面覆盖#theme_bw,X 和 Y 方向都有（灰色实线）,有助于界定绘图区域，在多子图中更清晰
   labs(
      x = "AWA External Score",
      y = "AWA Clinical Score",
      color = "Omics Combination",
      title = "External vs Clinical AWA Scores(SNF)"
   ) +
   scale_color_manual(values = distinct_colors) +   # <<< 更强色差
   theme(
      legend.position = "right",
      legend.text = element_text(size = 9), # <<< legend 字体更小
      legend.title = element_text(size = 10),
      
      plot.title = element_text(size =13.5, hjust = 0.5),#调整标题的大小
      
      axis.title = element_text(size = 12),   # 统一设置 x 和 y 轴标题字体, 坐标轴标题
      axis.text = element_text(size = 12),    # 统一设置 x 和 y 轴刻度字体 , 坐标轴刻度标签（数字或因子标签）
      # 移除网格线
      panel.grid.major = element_blank()     # 去掉主网格线
      #panel.grid.minor = element_blank() # 去掉次网格线
      
      # 关键：实现只有左+下边框
      # panel.border = element_blank(),      # 移除默认边框
      # panel.grid = element_blank(),        # 移除所有网格线（可选）
      # 
      # axis.line = element_line(color = "black", linewidth = 0.5),  # 添加 L 型轴线
      # axis.ticks = element_line(color = "black"),                  # 显示刻度线
      # axis.ticks.length = unit(0.15, "cm")                         # 刻度长度
   ) 
#+
   #coord_equal(xlim = c(0.2, 0.9),  # 设置 x 轴范围
               #ylim = c(0.2, 0.9) )  # 设置 y 轴范围，与 x 相同)  # 添加此行以使x轴和y轴的标尺一致
p      

ggsave(p,
       file="./External vs Clinical AWA Scores(SNF).png",
       width = 8,
       height = 5,
       units = "in",
       dpi = 300)

#有internal指标的3D图，可做为补充图，但它不是主线
library(scatterplot3d)
library(dplyr)
library(ggplot2)
# 生成高区分度颜色（和你 2D 图一致）
distinct_colors <- scales::hue_pal()(length(unique(snf_res$Combination_norm)))

distinct_colors <- c(
   "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
   "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#A6761D", "#666666",
   "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69",
   "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
)
names(distinct_colors) <- unique(snf_res$Combination_norm)

# 打开 PNG 文件保存
png("SNF_AWA_3D.png", width = 13, height = 7, units = "in", res = 300)

# ---- 分区布局：左为 3D 图，右为 legend ----
layout(matrix(c(1,2), nrow = 1), widths = c(4,1))

# ====== 左侧：3D 绘图 ======
par(mar = c(5,5,4,2))  # 正常边距

s3d <- scatterplot3d(
   x = snf_res$AWA_external,
   y = snf_res$AWA_clinical,
   z = snf_res$AWA_internal,
   color = distinct_colors[snf_res$Combination_norm],
   pch = 19,
   cex.symbols = 2,
   main = "External vs Clinical vs Internal",
   xlab = "AWA External",
   ylab = "AWA Clinical",
   zlab = "AWA Internal",
   angle = 55,
   box = TRUE
)

# ====== 右侧：独立的 legend 面板 ======
par(mar = c(5,0,4,2)) # 给右侧留位置

plot.new()   # 新建空白图
legend(
   "left",
   legend = names(distinct_colors),
   col = distinct_colors,
   pch = 19,
   cex = 0.9,
   bty = "n"
)

dev.off()

####-- 根据Avg_Rank原理的综合评价指标-----------------------
#因为 Avg_Rank 基于排名，不需要提前标准化，所以直接用你合并后的原始指标即可
##### ---去掉内部指标计算-- ####
library(dplyr)
# === 1. 定义指标及其方向 ===
# 聚类指标
clustering_higher_better <- c("NMI", "ARI", "Accuracy", "MacroF1", 
                                       "Avg_Purity")
clustering_lower_better  <- c( "Avg_Gini", "Avg_Entropy")

# 临床指标
clinical_higher_better <- c("ECP_num")        # 越大越好 
clinical_lower_better  <- c("log_rank_p")    # 越小越好

# === 2. 计算每个指标的排名 ===
df_ranked <- master %>%
   select(Combination_norm, 
          all_of(c(clustering_higher_better, clustering_lower_better,
                   clinical_higher_better, clinical_lower_better)))
#all_of()： 这是一个 dplyr 的选择帮手函数（selection helper）。
#用途： 当你需要通过一个存储在变量中的字符向量来选择列时，就必须使用它
# 初始化一个列表存储排名
rank_list <- list()
#rank() 函数默认是：数值越小，排名越小（rank 越小）
# 因此越大越好时，导致排名是反的

# 越大越好 → rank(-x)
for (col in c(clustering_higher_better, clinical_higher_better)) {
   rank_list[[paste0(col, "_rank")]] <- rank(-df_ranked[[col]], ties.method = "average")
}
# 越小越好 → rank(x)
for (col in c(clustering_lower_better, clinical_lower_better)) {
   rank_list[[paste0(col, "_rank")]] <- rank(df_ranked[[col]], ties.method = "average")
}

# 合并排名到数据框
df_with_ranks <- df_ranked %>%
   bind_cols(rank_list)

# === 3. 提取排名列用于 Avg_Rank 计算 ===
clustering_rank_cols <- paste0(c(clustering_higher_better, clustering_lower_better), "_rank")
clinical_rank_cols     <- paste0(c(clinical_higher_better, clinical_lower_better), "_rank")

# === 4. 计算 Avg_Rank（按原论文公式结构）===论文中 Item_* 均为数据集个数，我们只有一个数据集 → 都是 1
n_clustering <- length(clustering_rank_cols)  # 应为 7
n_clinical   <- length(clinical_rank_cols)    # 应为 2
df_with_ranks <- as.data.frame(df_with_ranks)
total_clustering_rank <- rowSums(df_with_ranks[clustering_rank_cols])
total_clinical_rank   <- rowSums(df_with_ranks[clinical_rank_cols])

avg_rank <- (total_clustering_rank / (2 * n_clustering)) +
   (total_clinical_rank   / (2 * n_clinical))

# === 5. 输出最终结果 ===
result <- master %>%
   select(Combination_norm) %>%
   bind_cols(avg_rank = avg_rank) %>%
   arrange(avg_rank)

#保存数据
save(master,df_ranked,df_with_ranks,result,file = "SNF_Avg_Rank.Rdata")
snf_rank <- result
save(snf_rank,file = "snf_rank.Rdata")

#综合评分avg_rank的画图
#排名越小越好 
library(ggplot2)
library(data.table)

# 提取用于绘图的数据
plot_data2 <- result[
   , .(Combination_norm, avg_rank)#.() 更高效，且是 data.table 推荐写法，取出这两列
][order(avg_rank)]  # 升序：低分在上，因为排名越小越好

# 可选：只显示 top N（例如 top 15），避免标签拥挤
# plot_data <- head(plot_data, 15)

# 创建组合名称的因子（确保排序正确）
# 设置因子：将当前顺序反转，使“最好”成为最后一个 level → 出现在顶部
#因为画图时分数大的在上，导致最优结果在最底部，所以我们需要反转一下
plot_data2[, Combination_norm := factor(Combination_norm, levels = rev(Combination_norm))]

# 绘图
p2 <- ggplot(plot_data2, aes(x = Combination_norm, y = avg_rank)) +
   geom_col(fill = "lightseagreen", width = 0.7) + 
   coord_flip() +  # 横向显示，便于阅读长名称
   labs(
      title = "Overall Performance Ranking of Multi-omics Combinations (SNF)",
      x = "Multi-omics Combination",
      y = "Rank Score (lower is better)"
   ) +
   theme_minimal(base_size = 14) +#白底 + 无边框,只有水平浅灰线（Y 轴方向）,突出数据本身
   theme(
      plot.title = element_text(size =13.5,hjust = 0.5, face = "bold"),#hjust = 0.5：水平居中对齐,face = "bold"：加粗字体
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),   # 统一设置 x 和 y 轴标题字体, 坐标轴标题
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
      
      # 关键：实现只有左+下边框
      # panel.border = element_blank(),      # 移除默认边框
      # panel.grid = element_blank(),        # 移除所有网格线（可选）
      # 
      # axis.line = element_line(color = "black", linewidth = 0.5),  # 添加 L 型轴线
      # axis.ticks = element_line(color = "black"),                  # 显示刻度线
      # axis.ticks.length = unit(0.15, "cm")                         # 刻度长度
   )
# 显示图形
print(p2)

# 保存
ggsave(plot = p2,
       file = "SNF_Avg_rank_BarPlot.png",  
       width = 9, 
       height = 7, 
       units = "in",
       dpi = 300)

#保存为PDF或svg方便缩小放大而不失真。导出为 PDF（矢量图，无限缩放不失真）
#为了观感一致，将宽度调整为一半6，高度保持不变，老师推荐还是一样就好，少放几个算法，26个组合的一个代表，51的放一个。
ggsave("SNF_Avg_rank_BarPlot.pdf", plot = p2, width = 9, height = 7, units = "in")



#### 基于avg_rank的结果展示不同组学的贡献度-------------
library(data.table)
# 创建每种组学的包含标志
result[, RNA_in := grepl("RNA", Combination_norm)]
result[, CNV_in := grepl("CNV", Combination_norm)]
result[, MET_in := grepl("MET", Combination_norm)]
result[, PRO_in := grepl("PRO", Combination_norm)]
result[, MIR_in := grepl("MIR", Combination_norm)]

# 使用 avg_rank2 计算（它是不包含内部评价指标的）
contribution <- result %>%
   summarise(
      RNA_effect = mean(avg_rank[RNA_in]) - mean(avg_rank[!RNA_in]),
      CNV_effect = mean(avg_rank[CNV_in]) - mean(avg_rank[!CNV_in]),
      MET_effect = mean(avg_rank[MET_in]) - mean(avg_rank[!MET_in]),
      PRO_effect = mean(avg_rank[PRO_in]) - mean(avg_rank[!PRO_in]),
      MIR_effect = mean(avg_rank[MIR_in]) - mean(avg_rank[!MIR_in])
   )
# 当在 data.table 的 j 参数（即 DT[i, j, by] 中的 j）中使用时，:= 的作用是创建（或修改）一列，并且是 原位（By Reference） 进行的。
#负值表示“加入该组学能显著降低 avg_rank2 → 性能提升”
# 值越大，贡献越
# contribution
# RNA_effect CNV_effect MET_effect PRO_effect MIR_effect
#    -2.054113   2.076623  -3.970346   3.047403  -2.630952
library(ggplot2)

contrib_df <- data.frame(
   Omics = c("RNA", "CNV", "MET", "PRO", "MIR"),
   Effect = c( -2.054113,2.076623,-3.970346,3.047403,-2.630952)  # 示例值，替换为你的真实计算结果
)

ggplot(contrib_df, aes(x = reorder(Omics, Effect), y = Effect)) +
   geom_col(fill = "#80B1D3") +
   scale_y_continuous(limits = c(-4, 4), expand = c(0, 0)) +  # 调整 y 轴范围
   coord_flip() +
   labs(
      x = "Omics Type",
      y = "Avg_Rank Reduction (negative = more helpful)",
      title = "Contribution of Each Omics Layer to Clustering Performance (SNF)"
   ) +
   theme_minimal(base_size = 12) +
   theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.y = element_blank()
   )

#####

























