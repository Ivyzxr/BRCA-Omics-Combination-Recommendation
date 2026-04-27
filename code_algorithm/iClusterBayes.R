#!/usr/bin/env Rscript

##### 0. 并行库 #####
library(parallel)

##### 1. 加载必要的包 #####
library(iClusterPlus)   # 提供 iClusterBayes
library(aricode)        # NMI, ARI
library(cluster)        # clustering utilities
library(caret)          # confusionMatrix（Accuracy & F1）
library(clue)           # solve_LSAP（匈牙利算法）
library(dplyr)          # 数据操作

##### 2. 加载数据 & Somatic 清洗 #####
load("/data1/zxr/BRCA/all_filter_zscore.Rdata")                        # rna_exp, met_exp_gene, snp_exp, prot_exp, mir_exp, cnv_cli
load("/data1/zxr/BRCA/TCGA_BRCA_CNV_Processed_icluster.Rdata")  # cnv_chrseg_final

mat <- as.matrix(snp_exp)
good_rows <- apply(mat, 1, function(x) any(x==0) && any(x==1))
somatic_clean <- t(mat[good_rows, , drop=FALSE])
#输入数据为行为样本，列为特征
##### 3. 构建 data_list & 真实标签 #####
data_list <- list(
   RNA         = t(rna_exp),
   Methylation = t(met_exp_gene),
   CNV         = cnv_chrseg_final_z,
   Somatic     = somatic_clean,
   PRO         = t(prot_exp),
   MIR         = t(mir_exp)
)

sample_ids <- rownames(data_list[[1]])
Truelabel   <- cnv_cli %>%
   filter(cnv_cli[,1] %in% sample_ids) %>%
   arrange(match(cnv_cli[,1], sample_ids)) %>%
   mutate(label = case_when(
      subtype == "Normal" ~ 1,
      subtype == "LumA"   ~ 2,
      subtype == "LumB"   ~ 3,
      subtype == "Her2"   ~ 4,
      subtype == "Basal"  ~ 5,
      TRUE                ~ NA_real_
   ))
Tlabel <- Truelabel$label
C      <- length(unique(na.omit(Tlabel)))  # 实际类别数

##### 4. 各组学分布映射 #####
type_map <- list(
   RNA         = "gaussian",
   Methylation = "gaussian",
   CNV         = "gaussian",
   Somatic     = "binomial",
   PRO         = "gaussian",
   MIR         = "gaussian"   # log2RPM 当作高斯
)

##### 5. 枚举所有组合 2~6 视图 #####
omics_names     <- names(data_list)
all_combinations <- unlist(
   lapply(2:length(omics_names), function(k)
      combn(omics_names, k, simplify = FALSE)
   ),
   recursive = FALSE
)

##### 6. 输出目录 & 结果容器 #####
output_dir <- "/data1/zxr/BRCA/iClusterBayesResults"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

##### 7. 并行 mclapply #####
n_cores <- 4
cat("使用", n_cores, "核并行运行 iClusterBayes ...\n")

results_list <- mclapply(all_combinations, function(combo) {
   combo_name <- paste(combo, collapse = "+")
   cat("→ 处理：", combo_name, "\n")
   on.exit(gc(), add=TRUE)
   
   tryCatch({
      # — 构造 iClusterBayes 参数 —
      sub_data <- lapply(combo, function(nm) data_list[[nm]])
      types    <- sapply(combo, function(nm) type_map[[nm]])
      L        <- length(sub_data)
      
      # 基础参数（可根据需要调整）
      args <- list(
         type        = types,
         K           = C - 1,
         n.burnin    = 1000,
         n.draw      = 1200,
         prior.gamma = rep(0.1, L),
         sdev        = 0.5,
         beta.var.scale = 1,
         thin        = 1,
         pp.cutoff   = 0.5
      )
      # 挂载 dt1...dt6
      for (i in seq_len(L)) {
         args[[paste0("dt", i)]] <- sub_data[[i]]
      }
      
      # 运行 iClusterBayes
      fit <- do.call(iClusterBayes, args)
      raw <- fit$clusters  # 长度 = 样本数
      
      # — 标签对齐 —
      map     <- solve_LSAP(table(Tlabel, raw), maximum = TRUE)
      aligned <- as.integer(map[raw])
      
      # — 保存标签文件 —对齐后的标签文件
      label_file <- file.path(output_dir,
                              paste0("Labels_iClusterBayes_", gsub("[^A-Za-z0-9]", "_", combo_name), ".csv"))
      write.csv(data.frame(Sample = sample_ids, Cluster = aligned),
                file = label_file, row.names = FALSE)
      
      # --- 保存未对齐的原始聚类标签（非常关键！）
      original_label_file <- file.path(output_dir, paste0("OriginalLabels_iClusterBayes_", combo_name, ".csv"))
      write.csv(
         data.frame(Sample = sample_ids,
                    Cluster = raw
         ), file = original_label_file, row.names = FALSE)
      
      # — 计算评价指标 —
      nmi_val <- NMI(Tlabel, aligned)
      ari_val <- ARI(Tlabel, aligned)
      cm      <- confusionMatrix(
         factor(aligned, levels = sort(unique(Tlabel))),
         factor(Tlabel,   levels = sort(unique(Tlabel)))
      )
      acc_val <- cm$overall["Accuracy"]
      f1s     <- suppressWarnings(cm$byClass[,"F1"])
      mac_f1  <- if (C > 2) mean(f1s, na.rm=TRUE) else f1s
      
      # — 返回结果行 —
      data.frame(
         Combination = combo_name,
         NMI         = as.numeric(nmi_val),
         ARI         = as.numeric(ari_val),
         Accuracy    = as.numeric(acc_val),
         MacroF1     = as.numeric(mac_f1),
         LabelFile   = label_file,
         NumClusters = length(unique(aligned)),
         stringsAsFactors = FALSE
      )
   }, error = function(e) {
      warning("组合 ", combo_name, " 出错：", conditionMessage(e))
      NULL
   })
}, mc.cores = n_cores, mc.preschedule = FALSE)

##### 8. 合并 & 保存 #####
results <- do.call(rbind, Filter(Negate(is.null), results_list))
write.csv(results,
          file.path(output_dir, "Omics_Combination_Evaluation_iClusterBayes.csv"),
          row.names = FALSE)
save(results,
     file = file.path(output_dir, "output_iClusterBayes.Rdata"))

cat("全部完成！结果已保存到：\n",
    file.path(output_dir, "Omics_Combination_Evaluation_iClusterBayes.csv"), "\n",
    file.path(output_dir, "output_iClusterBayes.Rdata"), "\n")


