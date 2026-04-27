# 加载必要的包
library(SNFtool)
library(aricode)
library(cluster)
library(caret)
library(gtools)
library(clue)
library(dplyr)

# 设定参数
K = 20   # 邻居数
alpha = 0.5  # 超参数
T = 20   # SNF算法的迭代次数
C = 5    # 聚类数目，对应PAM50
t = 20 #迭代次数

####-- SNF的算法输入数据的行为样本，列为特征------####
#SNF的算法中dist2函数适合计算连续型数据
load("/data1/zxr/BRCA/input_SNF.Rdata")
#数据在data_list
omics_names <- names(data_list)

# 标准化函数
Standard_Normalization <- function(x) {
   x <- as.matrix(x)
   mean <- apply(x, 2, mean)
   sd <- apply(x, 2, sd)
   sd[sd == 0] <- 1
   xNorm <- t((t(x) - mean) / sd)
   return(xNorm)
}
# 提前标准化数据
data_list <- lapply(data_list, Standard_Normalization)


# 生成所有组合（排除空组合）---------------------------------------------
all_combinations <- unlist(
   lapply(2:length(omics_names), function(k) combn(omics_names, k, simplify = FALSE)),
   recursive = FALSE
)

# 创建结果存储结构 ------------------------------------------------------
results <- data.frame(
   Combination = character(),
   NMI = numeric(),
   ARI = numeric(),
   Accuracy = numeric(),
   MacroF1 = numeric(),
   LabelFile = character(),
   stringsAsFactors = FALSE
)
# 动态确定聚类数（根据真实标签）-----------------------------------------
C <- length(unique(Tlabel))  # 示例中Tlabel有1~5，故C=5

# 循环处理不同的实验（比较算法性能和最佳组学组合筛选）
# **实验二：找到乳腺癌的最佳组学组合**
#多视图数据相似度网络融合
for (combo in all_combinations) {
   combo_name <- paste(combo, collapse = "+")
   
   # 获取当前组合对应的数据 ----------------------------------------------
   current_data <- lapply(combo, function(name) as.matrix(data_list[[name]]))
   
   # SNF融合流程 ---------------------------------------------------------
   tryCatch({
      # 计算距离矩阵
      dist_list <- lapply(current_data, function(x) dist2(x, x))
      
      # 构建相似性网络
      W_list <- lapply(dist_list, function(d) affinityMatrix(d, K, alpha))
      
      # 融合网络
      W_fused <- SNF(W_list, K, t)
      
      # 谱聚类
      group <- spectralClustering(W_fused, C)
      
      # 对齐标签（解决聚类编号与真实标签不一致问题）-----------------------
      confusion <- table(Tlabel, group)
      map <- solve_LSAP(confusion, maximum = TRUE)
      aligned_labels <- map[group]
      
      # 保存聚类标签（安全文件名）------------------------------------------
      output_dir <- "/data1/zxr/BRCA/SNFresults"
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)  # 确保递归创建目录
      # 生成文件名并保存，保存的是对齐后的标签
      safe_name <- gsub("[^A-Za-z0-9]", "_", combo_name)
      label_file <- file.path(output_dir, paste0("Labels_", safe_name, ".csv"))
      # 写入文件
      write.csv(data.frame(Sample = rownames(current_data[[1]]),Cluster = aligned_labels), label_file, row.names = FALSE)
     
       # --- 保存未对齐的原始聚类标签（非常关键！）
      original_label_file <- file.path(output_dir, paste0("OriginalLabels_SNF_", safe_name, ".csv"))
      write.csv(
         data.frame(Sample = rownames(current_data[[1]]),
                    Cluster = group
         ), file = original_label_file, row.names = FALSE)
      
      # 计算评价指标 ------------------------------------------------------
      nmi_value <- NMI(Tlabel, aligned_labels)
      nmi_value2 <- calNMI(aligned_labels, Tlabel)
      ari_value <- ARI(Tlabel, aligned_labels)
      
      confusion_matrix <- confusionMatrix(
         factor(aligned_labels, levels = sort(unique(Tlabel))),
         factor(Tlabel, levels = sort(unique(Tlabel)))
      )
      
      accuracy_value <- confusion_matrix$overall['Accuracy']
      macro_f1 <- ifelse(C > 2,
                         if (all(is.na(confusion_matrix$byClass[, "F1"]))) NA else mean(confusion_matrix$byClass[, "F1"], na.rm = TRUE),
                         confusion_matrix$byClass["F1"])
      #增加一个记录聚类数量的列到 results
      NumClusters <- length(unique(group))
      # 记录结果 ---------------------------------------------------------
      results <- rbind(results, data.frame(
         Combination = combo_name,
         NMI = nmi_value,
         NMI2 = nmi_value2,
         ARI = ari_value,
         Accuracy = accuracy_value,
         MacroF1 = macro_f1,
         LabelFile = label_file,
         NumClusters = NumClusters
      ))
      
      print(paste("Success:", combo_name))
   }, error = function(e) {
      print(paste("Error in", combo_name, ":", e$message))
   })
}
write.csv(results, "/data1/zxr/BRCA/Omics_Combination_Evaluation_SNF.csv", row.names = FALSE)
#save(results,file = "snf_output.Rdata")
save(results,file = "/data1/zxr/BRCA/output_SNF.Rdata")




