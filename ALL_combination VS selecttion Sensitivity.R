### 用于计算所有特征最终推荐的稳定性，
##  利用筛选特征后的得到的排名与全部特征下的排名
# 在相同 ranking 规则下，Full vs Selected feature 会不会改变结论
# 这两个 rank 必须是“已经做完跨算法归一化后的最终 rank，不是某一个算法内的 Avg_rank
# 
#导入数据
load("E:/Protocol_dataset/TCGAbiolinks/TCGA_ALL/pre_Final_rank.Rdata")
load("E:/Protocol_dataset/TCGAbiolinks/TCGA_ALL/sel_pre_Final_rank.Rdata")
#由于排名导致的组学名称顺序不一致
#我们是为了让特征筛选验证全部特征，所以将特征筛选的顺序调整和全部特征的一致再进行合并
#用的是cross_algo_rank和sel_cross_algo_rank
#调整组学顺序一致
sel_cross_algo_rank1 <- sel_cross_algo_rank[match(cross_algo_rank$Combination_canon,sel_cross_algo_rank$Combination_canon),]
#整合到同一个数据框中
library(dplyr)
df_rank1 <- left_join(cross_algo_rank,sel_cross_algo_rank1,by="Combination_canon")

df_rank <-data.frame(Combination = df_rank1$Combination_canon ,
                     norm_rank_full = df_rank1$mean_norm_rank ,
                     norm_rank_selected = df_rank1$sel_mean_norm_rank) 

#标注 Top10 / Top20 / Others（基于最终推荐）
library(ggplot2)
library(dplyr)
df_plot <- df_rank %>%
   mutate(
      rank_group = case_when(
         norm_rank_full <= sort(norm_rank_full)[15] ~ "Top 15",
         norm_rank_full <= sort(norm_rank_full)[25] ~ "Top 25",
         TRUE ~ "Others"
      )
   )
# 绘制 Rank–Rank scatter
p_rankrank <- ggplot(
   df_plot,
   aes(
      x = norm_rank_full,
      y = norm_rank_selected,
      color = rank_group
   )
) +
   geom_point(size = 2.5, alpha = 0.8) +
   geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey50"
   ) +
   scale_x_continuous(limits = c(0, 1)) +
   scale_y_continuous(limits = c(0, 1)) +
   scale_color_manual(
      values = c(
         "Top 15"  = "#D73027",
         "Top 25"  = "#FC8D59",
         "Others"  = "grey70"
      )
   ) +
   labs(
      title = "Sensitivity of multi-omics combination rankings to feature selection",
      x = "Normalized Avg_rank (Full features)",
      y = "Normalized Avg_rank (Selected features)",
      color = "Combination group"
   ) +
   theme_bw() +
   # theme(
   #    plot.title = element_text(hjust = 0.5, face = "bold"),
   #    legend.position = "bottom",
   #    panel.grid.minor = element_blank()
   # )
  theme(
   plot.title = element_text(size = 13.5, hjust = 0.5, face = "bold"),
   axis.text.x = element_text(face = "bold", size = 10),
   axis.text.y = element_text(size = 10),
   axis.title = element_text(size = 12),
   legend.position = "bottom",
   panel.grid.minor = element_blank()
  )
    
print(p_rankrank)

ggsave(plot = p_rankrank,
       file = "../../rankrank scatter Plot.png",  
       width = 8, 
       height = 6, 
       units = "in",
       dpi = 300)
#保存为PDF或svg方便缩小放大而不失真。导出为 PDF（矢量图，无限缩放不失真）
#为了观感一致，将宽度调整为一半6，高度保持不变
ggsave("../../rankrank scatter Plot.pdf", plot = p_rankrank, width = 8, height = 6, units = "in")
#保存数据
save(df_rank1,df_rank,df_plot,file ="../../Final_rankrank_scatter.Rdata" )


# 强烈推荐：算一个 Spearman 相关系数
rho <- cor(
   df_plot$norm_rank_full,
   df_plot$norm_rank_selected,
   method = "spearman"
)
rho # 0.1721366

df_top25 <- df_plot %>%
   filter(rank_group %in% c("Top 15", "Top 25"))
rho_top25 <- cor(
   df_top25$norm_rank_full,
   df_top25$norm_rank_selected,
   method = "spearman"
)
rho_top25 # 0.02462959

df_top15 <- df_plot %>%
   filter(rank_group == "Top 15")
rho_top15 <- cor(
   df_top15$norm_rank_full,
   df_top15$norm_rank_selected,
   method = "spearman"
)
rho_top15 #0.2076993

#Top K overlap（不是 Spearman）
#现在必须算这个：
top25_full <- df_plot %>%
   arrange(norm_rank_full) %>%
   slice(1:25) %>%
   pull(Combination)

top25_sel <- df_plot %>%
   arrange(norm_rank_selected) %>%
   slice(1:25) %>%
   pull(Combination)

jaccard_top25 <- length(intersect(top25_full, top25_sel)) /
   length(union(top25_full, top25_sel))
jaccard_top25 #0.3888889











