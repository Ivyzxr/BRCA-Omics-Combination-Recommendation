#不同组学组合的评价指标比较
#同一算法的不同组学组合整合的评价指标：

library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)

# ---1.读取数据
df <- read.csv("./SNF_orig.csv",stringsAsFactors = FALSE)
df1 <- df[,c(1:5)]
# 转换为长格式
library(tidyr)
df_long <- pivot_longer(df1, cols = -Combination_norm, 
                        names_to = "Metric", values_to = "Value")

#2. 分组柱状图
library(ggplot2)

p1 <- ggplot(df_long, aes(x = Combination_norm, y = Value, fill = Metric)) +
   geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
   geom_text(aes(label = ifelse(Value > 0.4, round(Value, 2), "")), 
             position = position_dodge(0.8), vjust = -0.5, size = 3) +
   scale_fill_brewer(palette = "Set2") +
   labs(title = "Performance Metrics Across Combinations(SNF)",
        x = "Combinations",
        y = "Value of the evaluation",
        fill = "Evaluation") +
   theme_minimal() +
   theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "right",
      panel.grid.major = element_line(colour = "grey90"),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
   )

ggsave(p1,
       filename = "./SNF_Combinations.png",
       width    = 12,#为了使观感一致
       height   = 7,#
       dpi      = 300,
       units    = "in"
)
#保存为PDF或svg方便缩小放大而不失真。导出为 PDF（矢量图，无限缩放不失真）
#为了观感一致，将宽度调整为一半6，高度保持不变,老师推荐还是一样就好，少放几个算法，26个组合的一个代表，51的放一个。
ggsave("SNF_Combinations.pdf", plot = p1, width = 12, height = 7, units = "in")

#3. 折线图（趋势分析）

# p2 <- ggplot(df_long, aes(x = Combination_norm, y = Value, group = Metric, color = Metric)) +
#    geom_line(linewidth = 1.2) +
#    geom_point(size = 3) +
#    geom_text(aes(label = round(Value, 2)), 
#              vjust = -1, size = 3, show.legend = FALSE,  # 禁止文本生成图例
#              check_overlap = TRUE) +
#    scale_color_brewer(palette = "Dark2") +
#    labs(title = "Scaling Trends in Combinatorial Omics Analysis",#组学组合扩展趋势分析
#         x = "Combinations",
#         y = "Value of the evaluation",
#         color = "Evaluation") +
#    theme_bw() +
#    theme(
#       plot.title = element_text(hjust = 0.5, face = "bold"),
#       axis.text.x = element_text(angle = 90, hjust = 1),#横坐标的标注为垂直于坐标轴
#       panel.grid.minor = element_blank(),
#       legend.position = "bottom",
#       plot.background  = element_rect(fill = "white", colour = NA),
#       panel.background = element_rect(fill = "white", colour = NA)
#    )
# 
# ggsave(p2,
#        filename = "./LRA_Combinations_tendency.png",
#        width    = 12,
#        height   = 7,
#        dpi      = 300,
#        units    = "in"
# )

##为了让趋势图按照2-3-4的顺序排列
library(dplyr)
library(stringr)
# df_long1 <- df_long %>%
#    mutate(
#       combo_size = str_count(Combination_norm, "\\+") + 1#需要对+做转义
#    )
df_long1 <- df_long %>%
   mutate(                                 #对每个拆分后的列表元素，计算其长度（即有多少个成分)
      combo_size = sapply(str_split(Combination_norm, "\\+"), length)#str_count(Combination_norm, "\\+") + 1
   ) %>%
   arrange(combo_size, Combination_norm) %>%#按组合大小排序，再按组合名称字母顺序排序
   mutate(
      Combination_norm = factor(#将 Combination_norm 列转换为因子（factor），并固定其水平（levels）为当前排序后的唯一值顺序
         Combination_norm,
         levels = unique(Combination_norm)
      )
   )
p2 <- ggplot(df_long1, aes(x = Combination_norm, y = Value, group = Metric, color = Metric)) +
   geom_line(linewidth = 1.2) +
   geom_point(size = 3) +
   geom_text(aes(label = round(Value, 2)), 
             vjust = -1, size = 3, show.legend = FALSE,  # 禁止文本生成图例
             check_overlap = TRUE) +
   scale_color_brewer(palette = "Dark2") +
   labs(title = "Scaling Trends in Combinatorial Omics Analysis(SNF)",#组学组合扩展趋势分析
        x = "Combinations",
        y = "Value of the evaluation",
        color = "Evaluation") +
   theme_bw() +
   theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1),#横坐标的标注为垂直于坐标轴
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
   )
ggsave(p2,
       filename = "./SNF_Combinations_tendency.png",
       width    = 12,
       height   = 7,
       dpi      = 300,
       units    = "in"
)

#保存为PDF或svg方便缩小放大而不失真。导出为 PDF（矢量图，无限缩放不失真）
#为了观感一致，将宽度调整为一半6，高度保持不变，老师推荐还是一样就好，少放几个算法，26个组合的一个代表，51的放一个。
ggsave("SNF_Combinations_tendency.pdf", plot = p2, width = 12, height = 7, units = "in")


#4. 雷达图（标准化处理）
library(ggradar)
library(scales)
library(fmsb)
library(RColorBrewer)
# 标准化数据预处理
df1 <- df[,c(1:5)]
df_radar <- df1
df_radar[, -1] <- apply(df_radar[, -1], 2, rescale)  # Min-Max标准化
# = 2. 生成带透明度的颜色（关键！）===
n_combos <- nrow(df_radar)
# 使用高区分度调色板，并添加透明度（alpha = 0.6）
base_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_combos)
transparent_colors <- paste0(base_colors, "99")  # "99" ≈ 60% 不透明度 (十六进制 alpha)
p3 <- ggradar(
   df_radar,
   values.radar = c("0", "0.5", "1"),
   grid.min = 0,
   grid.mid = 0.5,
   grid.max = 1,
   group.colours = transparent_colors,  # 颜色 + 透明度
   group.line.width = 0.3,               # ✅ 等价于 geom_path(size=0.3)
   group.point.size = 1.0,               # ✅ 等价于 geom_point(size=1.0)
   plot.title = "Standardized radar chart for metric balance analysis(SNF)",
   legend.title = "Combinations",
   legend.position = "right"
) +
   theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.text = element_text(size = 8,face = "bold"),# 字体缩小到 6pt
      legend.title = element_text(size = 10, face = "bold"),#设置图例标题的字体大小为8点，并且将其加粗显示，以突出图例标题的重要性。
      legend.box.margin = margin(t = -10, b = -10),# 紧凑图例
      legend.spacing.x = unit(0.2, "cm"), #设置图例项之间的水平间距为0.2厘米
      legend.key.size = unit(0.8, "lines"),# 控制图例中代表不同数据系列的颜色块的大小。这里设定为0.8行高
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
   )   

ggsave(p3,
       filename = "./SNF_Combinations_radar.png",
       width    = 18,
       height   = 9,
       dpi      = 300,
       units    = "in"
)
#保存为PDF或svg方便缩小放大而不失真。导出为 PDF（矢量图，无限缩放不失真）
#为了观感一致，将宽度调整为一半6，高度保持不变，老师推荐还是一样就好，少放几个算法，26个组合的一个代表，51的放一个。
ggsave("SNF_Combinations_radar.pdf", plot = p3, width = 18, height = 9, units = "in")

#### 不同组合的单独雷达图---##
df2 <- df %>%
   mutate(
      combo_size = sapply(str_split(.[[1]], "\\+"), length)
   )
plot_radar_by_group <- function(df_sub, title_suffix, file_suffix) {
   
   df_radar <- df_sub[, 1:5]
   df_radar[, -1] <- apply(df_radar[, -1], 2, rescale)
   
   n_combos <- nrow(df_radar)
   base_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_combos)
   transparent_colors <- paste0(base_colors, "99")
   
   p <- ggradar(
      df_radar, 
      values.radar = c("0", "0.5", "1"),
      grid.min = 0,
      grid.mid = 0.5,
      grid.max = 1,
      group.colours = transparent_colors,
      group.line.width = 0.3,
      group.point.size = 1.0,
      plot.title = paste(
         "Standardized radar chart for metric balance analysis(SNF)",
         title_suffix
      ),
      legend.title = "Combinations",
      legend.position = "right"
   ) +
      theme(
         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
         legend.text = element_text(size = 8, face = "bold"),
         legend.title = element_text(size = 10, face = "bold"),
         legend.key.size = unit(0.8, "lines"),
         plot.background  = element_rect(fill = "white", colour = NA),
         panel.background = element_rect(fill = "white", colour = NA)
      )
   
   ggsave(
      p,
      filename = paste0("./(SNF)Combinations_radar_", file_suffix, ".png"),
      width = 14,
      height = 10,
      dpi = 300,
      units = "in"
   )
}
# 2-omics
plot_radar_by_group(
   df2 %>% filter(combo_size == 2),
   "(2-omics combinations)",
   "2omics"
)

# 3-omics
plot_radar_by_group(
   df2 %>% filter(combo_size == 3),
   "(3-omics combinations)",
   "3omics"
)

# 4-omics
plot_radar_by_group(
   df2 %>% filter(combo_size == 4),
   "(4-omics combinations)",
   "4omics"
)

# 5–6 omics
plot_radar_by_group(
   df2 %>% filter(combo_size >= 5),
   "(5–6 omics combinations)",
   "5to6omics"
)

#### 对其他指标，放在补充图中，
# Purity，Entropy，Gini  画柱状图
#其中 Purity越大越好，而Entropy，Gini越小越好
#2. 分组柱状图
df_purity <- df[,c(1,9:11)]
# 转换为长格式
library(tidyr)
df_long2 <- pivot_longer(df_purity, cols = -Combination_norm, 
                         names_to = "Metric", values_to = "Value")
library(ggplot2)

p4 <- ggplot(df_long2, aes(x = Combination_norm, y = Value, fill = Metric)) +
   geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
   #geom_text(aes(label = ifelse(Value > 0.4, round(Value, 2), "")), 
   # position = position_dodge(0.8), vjust = -0.5, size = 3) +  #将数字标注去掉
   scale_fill_brewer(palette = "Set2") +
   labs(title = "Classification Purity Metrics Across Combinations(SNF)",
        x = "Combinations",
        y = "Value of the evaluation ",
        fill = "Evaluation") +
   theme_minimal() +
   theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "right",
      panel.grid.major = element_line(colour = "grey90"),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
   )

ggsave(p4,
       filename = "./SNF_class_combinations.png",
       width    = 12,
       height   = 7,
       dpi      = 300,
       units    = "in"
)

#保存为PDF或svg方便缩小放大而不失真。导出为 PDF（矢量图，无限缩放不失真）
#为了观感一致，将宽度调整为一半6，高度保持不变，老师推荐还是一样就好，少放几个算法，26个组合的一个代表，51的放一个。
ggsave("SNF_class_combinations.pdf", plot = p4, width = 12, height = 7, units = "in")


#### 对其他指标，放在补充图中，
# ECP  画柱状图

#2. 分组柱状图
df_ECP <- df[,c(1,12)]
# 转换为长格式
library(tidyr)
df_long3 <- pivot_longer(df_ECP, cols = -Combination_norm, 
                         names_to = "Metric", values_to = "Value")
library(ggplot2)

p5 <- ggplot(df_long3, aes(x = Combination_norm, y = Value, fill = Metric)) +
   geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
   #geom_text(aes(label = ifelse(Value > 0.4, round(Value, 2), "")), 
   # position = position_dodge(0.8), vjust = -0.5, size = 3) +  #将数字标注去掉
   scale_fill_brewer(palette = "Set2") +
   labs(title = "Number of Enriched Clinical Parameters(SNF)",
        x = "Combinations",
        y = "Number of ECPs",
        fill = "Evaluation") +
   theme_minimal() +
   theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "right",
      panel.grid.major = element_line(colour = "grey90"),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
   )

ggsave(p5,
       filename = "./SNF_ECP_combinations.png",
       width    = 12,
       height   = 7,
       dpi      = 300,
       units    = "in"
)
#保存为PDF或svg方便缩小放大而不失真。导出为 PDF（矢量图，无限缩放不失真）
#为了观感一致，将宽度调整为一半6，高度保持不变，老师推荐还是一样就好，少放几个算法，26个组合的一个代表，51的放一个。
ggsave("SNF_ECP_combinations.pdf", plot = p5, width = 12, height = 7, units = "in")


#ggplot2 的默认字体大小
# 主标题 (plot.title)	14 pt
# X/Y 轴标题 (axis.title.x, axis.title.y)	11 pt
# X/Y 轴刻度标签 (axis.text.x, axis.text.y)	11 pt
# 图例标题 (legend.title)	11 pt
# 图例项标签 (legend.text)	11 pt
#
#### 按组的数目排序
df_long31 <- df_long3  %>%
   mutate(                                 #对每个拆分后的列表元素，计算其长度（即有多少个成分)
      combo_size = sapply(str_split(Combination_norm, "\\+"), length)#str_count(Combination_norm, "\\+") + 1
   ) %>%
   arrange(combo_size, Combination_norm) %>%#按组合大小排序，再按组合名称字母顺序排序
   mutate(
      Combination_norm = factor(#将 Combination_norm 列转换为因子（factor），并固定其水平（levels）为当前排序后的唯一值顺序
         Combination_norm,
         levels = unique(Combination_norm)
      )
   )

p6 <- ggplot(df_long31, aes(x = Combination_norm, y = Value, group = Metric, color = Metric)) +
   geom_line(linewidth = 1.2) +
   geom_point(size = 3) +
   geom_text(aes(label = round(Value, 2)), 
             vjust = -1, size = 3, show.legend = FALSE,  # 禁止文本生成图例
             check_overlap = TRUE) +
   scale_color_brewer(palette = "Dark2") +
   labs(title = "Scaling Trends of Number of ECPs(SNF)",#组学组合扩展趋势分析
        x = "Combinations",
        y = "Number of ECPs",
        color = "Evaluation") +
   theme_bw() +
   theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1),#横坐标的标注为垂直于坐标轴
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
   )
ggsave(p6,
       filename = "./SNF_ECP_tendency.png",
       width    = 12,
       height   = 7,
       dpi      = 300,
       units    = "in"
)
