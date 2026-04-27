#### 临床评价指标 ####
# 1.生存分析   这里需要使用的对齐后的聚类标签
# 2.通路富集
#----------------------1.生存分析---------------- #
#---------这里需要使用的对齐后的聚类标签------------

# 1)加载R包
library(survival)
library(survminer)
library(dplyr)
library(stringr)
library(ggplot2)
library(readr)
library(survcomp)
# 2)载入数据

# 映射函数：数字标签 → PAM50 标签
label_to_PAM50 <- function(cluster_labels) {
   pam50_labels <- c("Normal", "LumA", "LumB", "Her2", "Basal")
   factor(pam50_labels[cluster_labels], levels = pam50_labels)
}

#读取标签数据
clinic_file <- load("./surv_data.Rdata")
survival_data<- surv_data
#设置聚类标签文件夹路径
label_dir <- "./SNFresults0817"
label_files <- list.files(label_dir, pattern = "^Labels_.*\\.csv$", full.names = TRUE)

#创建输出文件夹
output_dir <- file.path(label_dir, "KM_Plots")
dir.create(output_dir, showWarnings = FALSE)

# 6. 初始化结果记录表
summary_df <- data.frame(
   Method = character(),
   log_rank_p = numeric(),
   C_index = numeric(),
   stringsAsFactors = FALSE
)
   
# 7. 批量处理
for (label_file in label_files) {
   tryCatch({
   method_name <- str_remove_all(basename(label_file), "^Labels_|\\.csv$")

   # 读取聚类标签
   label_df <- read.csv(label_file)
   colnames(label_df) <- c("Sample", "Cluster")
   
   # 合并生存信息
   df <- merge(survival_data, label_df, by = "Sample")
   
   # 映射为 PAM50 标签
   df$PAM50_Label <- label_to_PAM50(df$Cluster)
   
   # 创建生存对象
   surv_obj <- Surv(time = df$time, event = df$event)
   
   # 生存差异 Log-rank 检验
   surv_diff <- survdiff(surv_obj ~ PAM50_Label, data = df)
   p_val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
   
   # C-index 计算（需将标签数字化，使用一致顺序）
   concordance <- concordance.index(
      x = df$Cluster, surv.time = df$time, surv.event = df$event, method = "noether"
   )
   c_index <- concordance$c.index
   
   # 绘图
   fit <- survfit(surv_obj ~ PAM50_Label, data = df)
   plot_path <- file.path(output_dir, paste0("KM_", method_name, ".png"))
   km_plot <- ggsurvplot(
      fit, data = df,
      pval = TRUE, pval.method = TRUE,
      legend.title = "PAM50 subtype",
      legend.labs = levels(df$PAM50_Label),
      risk.table = TRUE,
      title = paste("KM Curve -", method_name),
      ggtheme = theme_minimal(),
      
      # 1. 设置图主体（主图、标题、轴标签、刻度标签）的字体大小
      font.main = 14, #主标题
      font.x = c(12, "bold", "black"),#X 轴标签 的大小、样式和颜色
      font.y = c(12, "bold", "black"),#Y 轴标签 的大小、样式和颜色
      font.tickslab = c(10, "plain", "black"),#轴刻度上的数字/文字（如时间点和生存概率）的大小。
      font.legend = 8,#设置图例的字体大小
      
      # 2. 设置风险表 (Risk Table) 的字体大小
      risk.table.fontsize = 4 # ggsurvplot risk table 接受 0-10 的数字，4.5 左右接近 14pt 字体
   
   )
   ggsave(plot_path, km_plot$plot, width = 5, height = 4,units = "in",dpi = 300)
   
   # 记录结果
   summary_df <- rbind(summary_df, data.frame(
      Method = method_name,
      log_rank_p = round(p_val, 5),
      C_index = round(c_index, 5)
   ))
   
   message("Finished: ", method_name)
   
   }, error = function(e) {
      message("Error in ", label_file, ": ", e$message)
   })
}


   
# 8. 保存分析结果
write.csv(summary_df, file.path(output_dir, "SNF_KM_summary_results.csv"), row.names = FALSE)   
   
   
   
   
   
   
   



