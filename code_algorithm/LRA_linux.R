
# R环境为3.6--------------
library(LRAcluster)
#install.packages("aricode")
library(aricode)
#install.packages("cluster")
library(cluster)
#options(repos = c(RSPM = "https://packagemanager.posit.co/cran/2020-02-29"))
#install.packages("caret")
library(caret)
#install.packages("dplyr")
library(dplyr)
#install.packages("clue")
library(clue)

#导入数据
#LRA的数据本来有
#突变数据（二进制数据，伯努利分布）
#CNV数据、DNA methy(高斯分布)
#RNA数据（泊松分布）

load("/data1/zxr/BRCA/all_filter_zscore.Rdata")
#输入数据需要是list,且均为行为特征，列为样本
met_exp_gene <- as.matrix(met_exp_gene)
snp_exp <- as.matrix(snp_exp)
data_list <- list(
   RNA = rna_exp,          # 基因表达矩阵
   Methylation = met_exp_gene, # DNA 甲基化矩阵
   CNV = cnv_log2ratio,             # 拷贝数变异矩阵,这里用的是连续值
   Somatic = snp_exp,           # 二项数据
   PRO = prot_exp, #   高斯分布
   MIR = mir_exp  #  泊松分布
   )
# 定义每个组学的数据类型（命名列表）
types <- list(
   RNA = "gaussian",#标准化后为高斯分布
   Methylation = "gaussian",
   CNV = "gaussian",
   Somatic = "binary",
   PRO = "gaussian",
   MIR = "gaussian"
)
#确保样本名称一致的代码
# common_samples <- Reduce(intersect, lapply(data_list, colnames))
# data <- lapply(data_list, function(x) x[, common_samples])
# identical(common_samples,Truelabel$Sample_ID)

# 定义显示名称
display_names <- c(
   RNA = "RNA-seq",
   Methylation = "Methylation",
   CNV = "CNV",
   Somatic = "Somatic Mutation",
   PRO = "Proteomics",
   MIR = "MicroRNA"
)

# 处理真实标签
Truelabel <- cnv_cli[, c(1:3, 16)]
Truelabel <- Truelabel %>%
   mutate(
      label = case_when(
         subtype == "Normal" ~ 1,
         subtype == "LumA" ~ 2,
         subtype == "LumB" ~ 3,
         subtype == "Her2" ~ 4,
         subtype == "Basal" ~ 5,
         TRUE ~ NA_real_
      )
   )
Tlabel <- Truelabel$label  # 确保列名正确

# 定义参数 ------------------------------------------------------------------
output_dir <- "/data1/zxr/BRCA/LRAresults"  # 与SNF结果目录同级不同名
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

C <- length(unique(Tlabel))  # 真实标签类别数

# 生成所有非空组合--------------------
omics_names <- names(data_list)
all_combinations <- unlist(
   lapply(2:length(omics_names), 
          function(k) combn(omics_names, k, simplify = FALSE)),
   recursive = FALSE
)
# 初始化结果存储 ------------------------------------------------------------
results <- data.frame(
   Combination = character(),
   NMI = numeric(),
   ARI = numeric(),
   Accuracy = numeric(),
   MacroF1 = numeric(),
   LabelFile = character(),
   NumClusters = integer(),
   stringsAsFactors = FALSE
)

# 主分析流程 ----------------------------------------------------------------
for (comb in all_combinations) {
   comb_name <- paste(comb, collapse = "+")
   
   tryCatch({
      # ---- 数据准备 ----
      data_subset <- data_list[comb]
      types_subset <- sapply(comb, function(x) types[[x]])
      names_subset <- display_names[comb]
      
      # ---- 降维分析 ----
      rlist <- LRAcluster(
         data = data_subset,
         types = types_subset,
         names = names_subset,
         dimension = 2
      )
      
      # ---- 聚类分析 ----
      coord <- t(rlist$coordinate)  # 转置为样本x维度
      set.seed(123)  # 保证可重复性
      rclust <- kmeans(coord, centers = C, nstart = 20) # C为聚类的个数
      
      # ---- 标签对齐 ----
      # 标签对齐（使用SNF相同的solve_LSAP方法）
      confusion <- table(Tlabel, rclust$cluster)
      map <- solve_LSAP(confusion, maximum = TRUE)
      aligned_labels <- map[rclust$cluster]
      
      # --- 保存标签文件（格式与SNF完全相同）
      safe_name <- gsub("[^A-Za-z0-9]", "_", comb_name)
      label_file <- file.path(output_dir, paste0("Labels_LRA_", safe_name, ".csv"))
      write.csv(
         data.frame(Sample = colnames(data_list[[1]]),  # 假设样本名为矩阵列名
                    Cluster = aligned_labels
         ), file = label_file, row.names = FALSE)
      
      # --- 保存未对齐的原始聚类标签（非常关键！）
      original_label_file <- file.path(output_dir, paste0("OriginalLabels_LRA_", safe_name, ".csv"))
      write.csv(
         data.frame(Sample = colnames(data_list[[1]]),
                    Cluster = rclust$cluster
         ), file = original_label_file, row.names = FALSE)
      
      # ---- 评估指标 ----
      #   评价指标
      nmi_value <- NMI(Tlabel, aligned_labels)
      ari_value <- ARI(Tlabel, aligned_labels)
      
      confusion_matrix <- confusionMatrix(
         factor(aligned_labels, levels = sort(unique(Tlabel))),
         factor(Tlabel, levels = sort(unique(Tlabel)))
      )
      
      # 记录结果
      results <- rbind(results, data.frame(
         Combination = comb_name,
         NMI = nmi_value,
         ARI = ari_value,
         Accuracy = confusion_matrix$overall["Accuracy"],
         MacroF1 = ifelse(C > 2,
                          mean(confusion_matrix$byClass[, "F1"], na.rm = TRUE),
                          confusion_matrix$byClass["F1"]),
         LabelFile = label_file,
         NumClusters = length(unique(aligned_labels))  # 实际聚类数
      ))
      
      print(paste("Success:", comb_name))
      
   }, error = function(e) {
      print(paste("Error in", comb_name, ":", e$message))
   })
}

# 结果保存 ------------------------------------------------------------------
write.csv(results, "/data1/zxr/BRCA/Omics_Combination_Evaluation_LRA.csv", row.names = FALSE)
save(results, file = "/data1/zxr/BRCA/output_LRA.Rdata")








