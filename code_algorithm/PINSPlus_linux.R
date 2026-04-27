##导入我们自己的数据
#需要把数据整理为列表的形式（包含生存数据patientID,survival,death）

#需要看文献中对连续和离散变量的处理（没有说离散的）
#对离散数据进行改造,不太行，不做离散数据
#
library(PINSPlus)
library(aricode)  # NMI, ARI
library(caret)    # Accuracy, F1
library(dplyr)
library(clue)     # solve_LSAP

#多类型的数据
#输入的是列表，同样列表中的每个矩阵为：行是样本，列是特征；且每个矩阵必须拥有相同的样本
load("/data1/zxr/BRCA/all_filter_zscore.Rdata")
#输入数据需要是list,且均为行是样本，列是特征

# 数据预处理 ----------------------------------------------------------------
data_list <- list(
   RNA = t(as.matrix(rna_exp)),          # 基因表达矩阵（连续）
   Methylation = t(as.matrix(met_exp_gene)),  # 甲基化矩阵（连续）
   CNV = t(as.matrix(cnv_log2ratio)),    # CNV矩阵（连续）
   #Somatic = t(as.matrix(snp_exp)),      # 体细胞突变（二分类0/1）
   PRO = t(as.matrix(prot_exp)),         # 蛋白质组（连续）
   MIR = t(as.matrix(mir_exp))           # miRNA（连续）
)
# 检查样本数和名称一致性（修正检查逻辑）----------------------------------
# 获取所有数据集的样本名（假设转置后行名为样本名）
sample_names <- lapply(data_list, rownames)

# 验证样本名完全一致
stopifnot(
   "样本名称不一致" = all(sapply(sample_names, function(x) identical(x, sample_names[[1]]))
   ))
   
   # 验证样本数一致
   num_samples <- sapply(data_list, nrow)
   stopifnot(
      "样本数量不一致" = all(num_samples == num_samples[1])
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
table(Tlabel)
# Tlabel
# 1   2   3   4   5 
# 12 191 174  44  86 

# 定义全局参数 --------------------------------------------------------------
output_dir <- "/data1/zxr/BRCA/PINSresults"  # 结果保存目录（与SNF分开）
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
C <- length(unique(Tlabel))  # 根据真实标签确定聚类数

# 生成所有组合（2~6个组学）------------------------------------------------
omics_names <- names(data_list)
all_combinations <- unlist(
   lapply(2:length(omics_names), function(k) combn(omics_names, k, simplify = FALSE)),
   recursive = FALSE
)

# 初始化结果存储结构 ------------------------------------------------------
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

# 主分析流程 --------------------------------------------------------------
for (combo in all_combinations) {
   combo_name <- paste(combo, collapse = "+")
   
   tryCatch({
      # --- 步骤1：提取当前组合数据 ---
      current_data <- data_list[combo]
      # 检查维度和NA情况
      dims <- sapply(current_data, dim)
      nas <- sapply(current_data, function(x) any(is.na(x)))
      if (any(nas)) stop(paste("NA in", combo_name))
      
      print(paste("Processing:", combo_name))
      print(dims)
      
      # --- 步骤2：运行PINS分型 ---
      set.seed(123)  # 固定随机种子
      results_subtype <- SubtypingOmicsData(
         current_data,
         k = C, #聚类数，如果k确定，上面两项将被忽略
         clusteringMethod = "kmeans",
         clusteringOptions = list(nstart = 50),
         agreementCutoff = 0.5,#一致阈值。默认值为 0.5。
         ncore = 4,
         verbose = T,
         sampledSetSize = 2000,
         knn.k = NULL,
      )
      
      # --- 步骤3：标签对齐 ---
      raw_labels <- results_subtype$cluster2
      # 转换成因子后再变成整数标签（且保留raw_label的顺序）
      predicted_label <- as.integer(factor(raw_labels, levels = unique(raw_labels)))
      #标签对齐
      
      confusion <- table(Tlabel, predicted_label)
      map <- solve_LSAP(confusion, maximum = TRUE)
      aligned_labels <- as.integer(map)[predicted_label]
      
      # --- 步骤5：保存标签文件 ---
      safe_name <- gsub("[^A-Za-z0-9]", "_", combo_name)
      label_file <- file.path(output_dir, paste0("Labels_PINS_", safe_name, ".csv"))
      write.csv(
         data.frame(
            Sample = rownames(current_data[[1]]),  # 样本ID
            Cluster = aligned_labels                # 对齐后的标签
         ), 
         file = label_file, 
         row.names = FALSE
      )
      # --- 保存未对齐的原始聚类标签（非常关键！）
      original_label_file <- file.path(output_dir, paste0("OriginalLabels_PINS_", safe_name, ".csv"))
      write.csv(
         data.frame(Sample = rownames(current_data[[1]]),
                    Cluster = predicted_label
         ), file = original_label_file, row.names = FALSE)
      
      
      # --- 步骤6：计算评估指标 ---
      confusion_matrix <- confusionMatrix(
         factor(aligned_labels, levels = sort(unique(Tlabel))),
         factor(Tlabel, levels = sort(unique(Tlabel)))
      )
      new_row <- data.frame(
         Combination = combo_name,
         NMI = NMI(Tlabel, aligned_labels),
         ARI = ARI(Tlabel, aligned_labels),
         Accuracy = confusion_matrix$overall["Accuracy"],
         MacroF1 = if (C > 2) {
            mean(confusion_matrix$byClass[, "F1"], na.rm = TRUE)
         } else {
            confusion_matrix$byClass["F1"]
         },
         LabelFile = label_file,
         NumClusters = length(unique(aligned_labels))
      )         
         # --- 步骤7：更新结果 ---
         results <- rbind(results, new_row)
         print(paste("Success:", combo_name))
         
   }, error = function(e) {
      print(paste("Error in", combo_name, ":", e$message))
   })
}

# 保存结果文件 ------------------------------------------------------------
write.csv(results, file.path(output_dir, "Omics_Combination_Evaluation_PINS.csv"), 
          row.names = FALSE)
save(results, file = file.path(output_dir, "output_PINS.Rdata"))





