#算所有算法所有组合的总排名以及贡献度,以及标准化后的可比性
#### 1.计算总排名推荐-------------
load("E:/Protocol_dataset/TCGAbiolinks/TCGA_ALL/SNF/BRCA_SNF/snf_rank.Rdata")
load("E:/Protocol_dataset/TCGAbiolinks/TCGA_ALL/BCC/BCC/bcc_rank.Rdata")
load("E:/Protocol_dataset/TCGAbiolinks/TCGA_ALL/iClusterPlus/iClusterPlus/iCluster_rank.Rdata")
load("E:/Protocol_dataset/TCGAbiolinks/TCGA_ALL/LRA/lra_rank.Rdata")
load("E:/Protocol_dataset/TCGAbiolinks/TCGA_ALL/PINSPlus/PINSPlus/pins_rank.Rdata")

library(dplyr)
library(tidyr)

rank_list <- list(
   BCC      = bcc_rank,
   iCluster = icluster_rank,
   LRA      = lra_rank,
   PINS     = pins_rank,
   SNF      = snf_rank
)
#1.单算法内归一化名次，解决不同算法组合数量不同，不同算法rank尺度不可比的问题
rank_long <- bind_rows( #将 lapply 处理完的多个小表格，像叠罗汉一样垂直拼接成一个名为 rank_long 的大表
   lapply(names(rank_list), function(a) { #代码遍历 rank_list 中的每一个算法名称
      df <- rank_list[[a]]# 取出数据，拿到该算法对应的原始排名表
      
      df %>%
         arrange(avg_rank) %>%  #重新排序，确保在该算法内部，表现最好的组合排在第一行 ，avg_rank 越小越好
         mutate(
            Algorithm = a, # 标记这个数据来自哪个算法
            Rank = row_number(), #根据排序结果分配一个整数排名
            N_a = n(),  #算法a支持的组合总数，计算该算法一共测试了多少个组合
            NormRank = (Rank - 1) / (N_a - 1)#核心的数学处理
         ) %>% #消除不同算法因组合数量不同导致的排名尺度差异
         select(Combination_norm, Algorithm, Rank, N_a, NormRank)
   })
)
#用的是 Rank-based normalization，不是数值归一化

#2.跨算法聚合，计算该组合在所有可支持它的算法中的平均相对排名
#由于类似像RNA+MET;MET+RNA会被识别为两种组合，因此我们需要先解决它的组合名排序导致的问题
# 定义组合规范化函数
normalize_combo <- function(x) {
   parts <- unlist(strsplit(x, "\\+")) #根据加号 + 将字                                                                                                       符串拆分成多个部分
   parts <- trimws(parts) #去除字符串首尾的空格
   paste(sort(parts), collapse = "+") #将拆分出来的部分按字母顺序重新排序
}
#在 long 表中生成“规范组合名”
#这一步一定要在 group_by 之前
rank_long <- rank_long %>%
   mutate(
      Combination_canon = sapply(Combination_norm, normalize_combo)
   )
#sapply 的作用： 它像一个“复印机”或者“循环器”，把 Combination_norm 这一列里的每一个组合名称，
#挨个儿丢进 normalize_combo 函数里处理，然后再把结果收集起来，拼成一个新的向量。

##---- 中间加的results需要给出5个算法的归一化的结果-----

# 跨算法聚合
library(dplyr)
cross_algo_rank <- rank_long %>%
   group_by(Combination_canon) %>%
   summarise(
      mean_norm_rank = mean(NormRank),#该组合在所有算法中的平均归一化排名（越小越好）
      sd_norm_rank   = sd(NormRank),#排名波动性（标准差
      n_algorithm    = n_distinct(Algorithm),#参与评估的算法数量（有些组合可能某些算法不支持）
      .groups = "drop"#避免保留分组信息
   ) %>%
   arrange(mean_norm_rank)#按平均排名升序排列 → 排名越靠前（值越小）的组合排在前面
# 标记 Top 10 / Top 20（画图用）
cross_algo_rank <- cross_algo_rank %>%
   mutate(
      rank_global = row_number(),#给排序后的组合编号（第1名到第N名）
      rank_group = case_when( 
         rank_global <= 10 ~ "Top 10",#Top10” 是指“全局综合表现最好的10个组合”，不是“某个算法下的Top10”
         rank_global <= 20 ~ "Top 20",
         TRUE ~ "Others"
      )
   )

cross_algo_rank$rank_group
#把这个信息 merge 回 rank_long：
rank_long2 <- rank_long %>%
   left_join(
      cross_algo_rank %>%
         select(Combination_canon, rank_group),
      by = "Combination_canon"
   )

#高亮TOP10/20
p_5pre <- ggplot(rank_long2,
       aes(x = Algorithm,
           y = NormRank,
           color = rank_group)) +
   
   geom_boxplot(
      outlier.shape = NA,
      fill = "grey95",
      color = "grey50",
      width = 0.5
   ) +
   
   geom_jitter(
      width = 0.25,
      alpha = 0.8,
      size = 2.5
   ) +
   
   scale_color_manual(
      values = c(
         "Top 10" = "#D73027",
         "Top 20" = "#FC8D59",
         "Others" = "grey70"
      )
   ) +
   
   labs(
      title = "Cross-algorithm comparison using normalized composite ranking",
      x = "Algorithm",
      y = "Normalized Avg_rank",
      color = "Global ranking tier"
   ) +
   
   theme_bw() +
   # theme(
   #    plot.title = element_text(hjust = 0.5, face = "bold"),
   #    axis.text.x = element_text(face = "bold")
   # )
  theme(
   plot.title = element_text(size = 13.5, hjust = 0.5, face = "bold"),
   axis.text.x = element_text(face = "bold", size = 10),
   axis.text.y = element_text(size = 10),
   axis.title = element_text(size = 12)
   #panel.grid.major.y = element_blank()
   #panel.grid.minor = element_blank()
  )

p_5pre

ggsave(plot = p_5pre,
       file = "../../pre_Final_Recommendation_jitter_Plot.png",  
       width = 10, 
       height = 7, 
       units = "in",
       dpi = 300)
#保存为PDF或svg方便缩小放大而不失真。导出为 PDF（矢量图，无限缩放不失真）
#为了观感一致，将宽度调整为一半6，高度保持不变
ggsave("../../pre_Final_Recommendation_jitter_Plot.pdf", plot = p_5pre, width = 10, height = 7, units = "in")
#保存数据
save(cross_algo_rank,rank_long2,rank_long,file ="../../pre_Final_rank.Rdata" )


#3.用规范的组合名做跨算法聚合
final_rank <- rank_long %>%
   group_by(Combination_canon) %>% #将数据按照 Combination_canon 这一列进行分组
   summarise(
      Final_AvgRank = mean(NormRank),#计算归一化排名（NormRank）的平均值。
      Support_Algorithms = n_distinct(Algorithm),#计算该组合涉及到了多少个不重复的算法
      Algorithms = paste(sort(unique(Algorithm)), collapse = ", "),#提取该组合对应的所有算法名称，按字母顺序排序，用逗号把这些名称串联成一个字符串
      .groups = "drop"#计算完成后，取消分组状态。这是一个良好的编程习惯，防止后续操作意外受到分组影响。
   ) %>%
   arrange(Final_AvgRank)#按照平均排名从低到高（升序）排列。
#保存最终排名
save(final_rank,file = "../../Final_rank.Rdata")
write.csv(final_rank, file = "../../Final_rank_Summary.csv", row.names = FALSE)
# final_rank <- rank_long %>%
#    group_by(Combination_norm) %>%
#    summarise(
#       Final_AvgRank = mean(NormRank),
#       Support_Algorithms = n(),
#       Algorithms = paste(sort(unique(Algorithm)), collapse = ", "),
#       .groups = "drop"
#    ) %>%
#    arrange(Final_AvgRank)
# Final_AvgRank：最终推荐分数（越小越好）
# Support_Algorithms：该组合被多少算法支持（非常重要，建议报告）
# Algorithms：支持该组合的算法列表（可用于补充材料）

library(ggplot2)
topN <- 20
TOP20 <- ggplot(final_rank[1:topN, ],
       aes(x = reorder(Combination_canon, -Final_AvgRank),#reorder(Combination_canon, Final_AvgRank) 是按数值从小到大排列
           y =  Final_AvgRank,
           fill = Support_Algorithms)) +
   geom_col(width = 0.7) +
   coord_flip() + #会翻转视觉方向，导致「最小值在最下面」 y 轴（原本的 x 轴）是从下往上绘制的,所以最小值（排名第 1 的组合）会出现在图表的最底部。
   scale_fill_gradient(low = "#2c7bb6", high = "#d7191c") +
   labs(
      x = "Omics Combination",
      y = "Final Normalized Rank Score(lower is better)",
      fill = "#Algorithms",
      title = "Final Recommendation of Multi-omics Combinations for BRCA(Top20)"
   ) +
   theme_classic(base_size = 14)
TOP20
ggsave(plot = TOP20,
       file = "../../Final_Recommendation_Top20_BarPlot.png",  
       width = 10, 
       height = 7, 
       units = "in",
       dpi = 300)
#保存为PDF或svg方便缩小放大而不失真。导出为 PDF（矢量图，无限缩放不失真）
#为了观感一致，将宽度调整为一半6，高度保持不变
ggsave("../../Final_Recommendation_Top20_BarPlot.pdf", plot = TOP20, width = 10, height = 7, units = "in")


p1 <- ggplot(final_rank,
             aes(x = reorder(Combination_canon, -Final_AvgRank),#reorder(Combination_canon, Final_AvgRank) 是按数值从小到大排列
                 y =  Final_AvgRank,
                 fill = Support_Algorithms)) +
   geom_col(width = 0.7) +
   coord_flip() + # y 轴（原本的 x 轴）是从下往上绘制的,所以最小值（排名第 1 的组合）会出现在图表的最底部。
   scale_fill_gradient(low = "#2c7bb6", high = "#d7191c") +
   labs(
      x = "Omics Combination",
      y = "Final Normalized Rank Score(lower is better)",
      fill = "#Algorithms",
      title = "Final Recommendation of Multi-omics Combinations for BRCA"
   ) +
   theme_classic(base_size = 14)
p1

ggsave(plot = p1,
       file = "../../Final_Recommendation_BarPlot.png",  
       width = 10, 
       height = 9, 
       units = "in",
       dpi = 300)

#保存为PDF或svg方便缩小放大而不失真。导出为 PDF（矢量图，无限缩放不失真）
#为了观感一致，将宽度调整为一半6，高度保持不变
ggsave("../../Final_Recommendation_BarPlot.pdf", plot = p1, width = 10, height = 9, units = "in")


####2.根据排名推荐统计组学贡献度，用于回答对于乳腺癌分型来说，哪些组学是核心的，哪些是可选，哪些贡献度有限-----
#根据Final_rank计算组学贡献度
#1.把组合拆成组学
library(dplyr)
library(tidyr)
ocs_long <- final_rank %>%
   mutate(Contribution = 1 - Final_AvgRank) %>%   # Rank → Contribution
   separate_rows(Combination_canon, sep = "\\+") %>%
   rename(Omics = Combination_canon)
#2.计算OCS
ocs_result <- ocs_long %>%
   group_by(Omics) %>%
   summarise(
      OCS = mean(Contribution),
      Frequency = n(),
      .groups = "drop"
   ) %>%
   arrange(desc(OCS))

#3.组学贡献度柱状图
library(ggplot2)

p <- ggplot(ocs_result,
            aes(x = reorder(Omics, OCS),
                y = OCS,
                fill = OCS)) +
   geom_col(width = 0.6) +
   coord_flip() +
   scale_fill_gradient(low = "#91bfdb", high = "#d73027") +
   labs(
      x = "Omics Type",
      y = "Omics Contribution Score (OCS)",
      title = "Omics Contribution to BRCA Subtyping",
      subtitle = "Based on cross-algorithm normalized Rank_Score"
   ) +
   theme_classic(base_size = 14)
p

ggsave(plot = p,
       file = "../../Omics_Contribution_BRCA_Subtype_BarPlot.png",  
       width =6 , 
       height = 5, 
       units = "in",
       dpi = 300)

#保存为PDF或svg方便缩小放大而不失真。导出为 PDF（矢量图，无限缩放不失真）
#为了观感一致，将宽度调整为一半6，高度保持不变
ggsave("../../Omics_Contribution_BRCA_Subtype_BarPlot.pdf", plot = p, width = 6, height = 5, units = "in")


####-2-计算TOP20中的组学贡献度
#1.选取top20最优组合
topN <- 20
top_combos <- final_rank %>%
   arrange(Final_AvgRank) %>%
   slice(1:topN)
#2.拆分组合 → 单个组学
library(tidyr)
library(dplyr)

omics_freq <- top_combos %>%
   separate_rows(Combination_canon, sep = "\\+") %>%
   group_by(Combination_canon) %>%
   summarise(Frequency = n(), .groups = "drop")

#3.计算相对贡献度
omics_freq <- omics_freq %>%
   mutate(
      Contribution = Frequency / sum(Frequency)
   ) %>%
   arrange(desc(Contribution))

#4.柱状图
library(ggplot2)

p <- ggplot(omics_freq,
       aes(x = reorder(Combination_canon, Contribution),
           y = Contribution,fill = Contribution)) +
   geom_col(width = 0.7) +   
   scale_fill_gradient(low = "#2c7bb6", high = "#d73027")+
   coord_flip() +
   labs(
      x = "Omics Type",
      y = "Contribution Frequency",
      title = paste0(
         "Omics Contribution Frequency in Top20")
   ) +
   theme_classic(base_size = 14)
p

ggsave(plot = p,
       file = "../../Top20_frequency.png",  
       width = 6, 
       height = 5, 
       units = "in",
       dpi = 300)

#保存为PDF或svg方便缩小放大而不失真。导出为 PDF（矢量图，无限缩放不失真）
#为了观感一致，将宽度调整为一半6，高度保持不变
ggsave("../../Top20_frequency.pdf", plot = p, width = 6, height = 5, units = "in")

save(ocs_result,ocs_long,top_combos,omics_freq,file = "../../Final_rank_contribution_frequency.Rdata")









