# 计算每个类的 purity、gini impurity
# 全局平均 purity、gini impurity
# -----这里用的标签是原始的标签，而非对齐后的标签--------
# 类内一致性评价模块
#用于评价
#每个聚类方法输出的 cluster 内部是否“纯”，是否对应同一种 PAM50 生物学亚型


# 加载必要R包
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
# install.packages("DescTools")
library(DescTools)  # For Gini

# 设置路径
label_dir <- "./SNFresults0817/"
output_dir <- file.path(label_dir, "Original_Purity_Gini_Entropy")
dir.create(output_dir, showWarnings = FALSE)
# 
# # PAM50 映射函数
# label_to_PAM50 <- function(cluster_labels) {
#    pam50_labels <- c("Normal", "LumA", "LumB", "Her2", "Basal")
#    factor(pam50_labels[cluster_labels], levels = pam50_labels)
# }

# 加载生存数据，含真实PAM50标签
load("./surv_data.Rdata")  # 应包含 Sample, PAM50_label 或 subtype 列
surv_data1 <- cnv_cli %>%
   select(Sample = Sample_ID,
          subtype = subtype)
# surv_data <- na.omit(surv_data)
# surv_data$PAM50_Label <- factor(cnv_cli$subtype, levels = c("Normal", "LumA", "LumB", "Her2", "Basal"))

# 定义计算 Entropy 函数
shannon_entropy <- function(p) {
   -sum(p * log2(p + 1e-10))  # 防止log(0)
}

# 获取所有原始聚类标签文件
label_files <- list.files(label_dir, pattern = "^OriginalLabels_.*\\.csv$", full.names = TRUE)

# 初始化结果收集
all_summary <- data.frame()

# 主循环
for (file in label_files) {
   method_name <- str_remove_all(basename(file), "^OriginalLabels_|\\.csv$")
   
   # 读取聚类结果
   label_df <- read.csv(file)
   colnames(label_df) <- c("Sample", "Cluster")
   
   # 合并真实标签
   #df <- merge(label_df, surv_data1, by = "Sample")样本顺序会变
   df <- left_join(label_df, surv_data1, by = "Sample")
   df$Cluster <- factor(df$Cluster)
   
   # 每个聚类内的指标
   stat <- df %>%
      group_by(Cluster) %>%
      summarise(n = n(),
                purity = max(table(subtype) / n),#该类中最大占比的 PAM50 亚型（纯度指标）
                gini = Gini(table(subtype) / n), #类内分布的 Gini impurity（越低表示越“纯”）
                entropy = shannon_entropy(table(subtype) / n)) %>%
      mutate(Method = method_name)
   
   # 加权平均（全局指标）
   summary_row <- data.frame(
      Method = method_name,
      Avg_Purity = weighted.mean(stat$purity, stat$n),
      Avg_Gini = weighted.mean(stat$gini, stat$n),
      Avg_Entropy = weighted.mean(stat$entropy, stat$n)
   )
   
   all_summary <- rbind(all_summary, summary_row)
   
   # 保存每方法详细信息
   write.csv(stat, file.path(output_dir, paste0("Purity_Gini_Entropy_", method_name, ".csv")), row.names = FALSE)
}

# 保存汇总
write.csv(all_summary, file.path(output_dir, "Purity_Gini_Entropy_SNF_Summary.csv"), row.names = FALSE)
   
 
 
   
  


