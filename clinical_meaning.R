#临床标签富集分析
#---------这里使用的是对齐后的聚类标签----------

#离散参数的富集：比如乳腺癌的ER/PR/HER2状态、组织学分级等分类变量，
#使用卡方检验（Chi-square test）来评估聚类结果与这些临床标签之间的关联性

#数值参数的富集：比如年龄等连续变量。使用Kruskal-Wallis检验，来评估聚类结果与这些数值变量之间的关联性

library(dplyr)
library(stringr)
library(readr)
library(tidyr)

# 设置路径
label_dir <- "./SNFresults0817/"
output_dir <- file.path(label_dir, "clinical_ECP")
dir.create(output_dir, showWarnings = FALSE)


# === 输入部分 ===
label_dir <- "./SNFresults0817"  # 替换为你的聚类标签文件夹路径
label_files <- list.files(label_dir, pattern = "^Labels_.*\\.csv$", full.names = TRUE)

load("./surv_data.Rdata") 
clinical <- cli_sur %>%
   select(Sample = Sample_ID,
          ER =ER_Status_nature2012,
          PR =PR_Status_nature2012,
          HER2= HER2_Final_Status_nature2012,
          Stage=Stage,
          T_stage=T_stage,
          N_stage=N_stage,
          M_stage=M_stage,
          Age=Age)# 包含 Sample 列 + 多个临床变量

discrete_vars <- c("ER", "PR", "HER2","Stage", "T_stage", "N_stage", "M_stage") # 离散变量
continuous_vars <- c("Age")  # 连续变量



# 用于存放详细结果
detailed_results <- list()
# 用于存放汇总结果
summary_results <- data.frame()

alpha <- 0.05  # 显著性阈值

for (label_file in label_files) {
   tryCatch({
      # 组合名去掉 Labels_ 和 .csv
      method_name <- str_remove_all(basename(label_file), "^Labels_|\\.csv$")
      
      # 读标签
      df_labels <- read_csv(label_file, show_col_types = FALSE) %>% drop_na()
      
      # 确保有 Sample 列
      if (!"Sample" %in% colnames(df_labels)) stop("缺少 Sample 列: ", label_file)
      
      # 合并临床信息
      merged_df <- inner_join(df_labels, clinical, by = "Sample") %>% drop_na()
      cluster_labels <- merged_df$Cluster
      
      # 离散型变量检验：卡方检验
      discrete_pvals <- sapply(discrete_vars, function(var) {
         tbl <- table(cluster_labels, merged_df[[var]])
         if (all(dim(tbl) > 1)) {
            suppressWarnings(chisq.test(tbl)$p.value)
         } else {
            NA
         }
      })
      
      # 连续型变量检验：Kruskal-Wallis
      continuous_pvals <- sapply(continuous_vars, function(var) {
         suppressWarnings(kruskal.test(merged_df[[var]] ~ cluster_labels)$p.value)
      })
      
      # 整理详细结果（长格式）
      all_pvals <- c(discrete_pvals, continuous_pvals)
      detailed_df <- data.frame(
         Omics_Combination = method_name,
         Variable = names(all_pvals),
         P_value = all_pvals,
         stringsAsFactors = FALSE
      )
      
      # 添加到详细结果列表
      detailed_results[[method_name]] <- detailed_df
      
      # 汇总统计显著变量个数
      sig_discrete <- sum(discrete_pvals < alpha, na.rm = TRUE)
      sig_continuous <- sum(continuous_pvals < alpha, na.rm = TRUE)
      total_sig <- sig_discrete + sig_continuous
      
      summary_results <- bind_rows(summary_results, data.frame(
         Omics_Combination = method_name,
         Significant_Discrete = sig_discrete,
         Significant_Continuous = sig_continuous,
         Total_Significant = total_sig
      ))
      
   }, error = function(e) {
      message("处理文件时出错: ", label_file, " -> ", e$message)
   })
}

# 生成详细表（宽格式：CNV+Meth  ER  p-value  PR  p-value ...）
detailed_wide <- bind_rows(detailed_results) %>%
   pivot_wider(
      names_from = Variable,
      values_from = P_value
   )

# 保存结果
write_csv(detailed_wide,file.path(output_dir, "clinical_enrichment_SNF_detailed.csv"))
write_csv(summary_results, file.path(output_dir,"clinical_enrichment_SNF_summary.csv"))

message("分析完成！已生成 clinical_enrichment_detailed.csv 和 clinical_enrichment_summary.csv")











