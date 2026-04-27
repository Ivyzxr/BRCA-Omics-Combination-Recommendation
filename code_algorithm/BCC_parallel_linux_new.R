##### 1. 加载必要的包 #####
library(bayesCC)    # BCC 算法
library(aricode)    # NMI, ARI
library(cluster)    # clustering utilities（如 pam/kmeans）
library(caret)      # confusionMatrix（Accuracy & F1）
library(dplyr)      # 数据操作
library(clue)       # solve_LSAP（匈牙利算法）
library(foreach)    # 并行循环
library(doParallel) # 并行后端

##### 2. 数据加载 & 初步预处理 #####
load("/data1/zxr/BRCA/all_filter_zscore.Rdata")  
# 假设 all_filter.Rdata 中包含：rna_exp, met_exp_gene, cnv_log2ratio, prot_exp, mir_exp, cnv_cli

# 构建数据列表（行=特征，列=样本）
data_list <- list(
   RNA          = as.matrix(rna_exp),
   Methylation  = as.matrix(met_exp_gene),
   CNV          = as.matrix(cnv_log2ratio),
   PRO          = as.matrix(prot_exp),
   MIR          = as.matrix(mir_exp)
)

# 对齐“公共”样本（列名）并剔除任何意外缺失
common_samples <- Reduce(intersect, lapply(data_list, colnames))
data_list <- lapply(data_list, function(mat) mat[, common_samples, drop = FALSE])

# 准备真实标签向量
# 假设 cnv_cli 的第一列是 Sample_ID，与上面 colnames() 对应
sample_id_col <- colnames(cnv_cli)[1]
Truelabel <- cnv_cli %>%
   filter(.data[[sample_id_col]] %in% common_samples) %>%
   mutate(label = case_when(
      subtype == "Normal" ~ 1,
      subtype == "LumA"   ~ 2,
      subtype == "LumB"   ~ 3,
      subtype == "Her2"   ~ 4,
      subtype == "Basal"  ~ 5,
      TRUE                ~ NA_real_
   ))

# 按 common_samples 的顺序排列
Truelabel <- Truelabel[match(common_samples, Truelabel[[sample_id_col]]), ]
Tlabel <- Truelabel$label  # 长度 = length(common_samples)
C <- length(unique(na.omit(Tlabel)))  # 聚类数

##### 3. 并行环境设置 #####
n_cores <- 8
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Registered", n_cores, "cores for parallel execution.\n")

##### 4. 生成所有 2~5 视图组合 #####
omics_names    <- names(data_list)
all_combinations <- unlist(
   lapply(2:length(omics_names), function(k) combn(omics_names, k, simplify = FALSE)),
   recursive = FALSE
)

##### 5. 主并行循环：对每个组合运行 bayesCC #####
output_dir <- "/data1/zxr/BRCA/BCCresults"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

results <- foreach(combo = all_combinations,
                   .combine = rbind,
                   .packages = c("bayesCC","aricode","clue","caret","dplyr")) %dopar% {
                      combo_name      <- paste(combo, collapse = "+")
                      safe_combo_name <- paste(combo, collapse = "_")
                      
                      # 提取当前组合的数据子集
                      current_data <- setNames(
                         lapply(combo, function(nm) data_list[[nm]]),
                         combo
                      )
                      
                     
                      # 执行 BCC
                      out <- tryCatch({
                         cat("Running BCC for:", combo_name, "...\n")
                         bcc_res <- bayesCC(X = current_data,
                                            K = C,
                                            IndivAlpha = TRUE,
                                            maxiter = 1000)
                         
                         # 提取硬聚类标签
                         bestC  <- bcc_res$Cbest              # Nsamples × K 矩阵
                         labels <- apply(bestC, 1, which.max) # 长度 = Nsamples
                         
                         # 对齐标签（匈牙利算法）
                         tab   <- table(Tlabel, labels)
                         map   <- solve_LSAP(tab, maximum = TRUE)
                         aligned_labels <- as.integer(map[labels])
                         
                         # 保存标签文件----对齐后的标签
                         label_file <- file.path(output_dir, paste0("Labels_BCC_", safe_combo_name, ".csv"))
                         write.csv(data.frame(Sample = common_samples,
                                              Cluster = aligned_labels),
                                   file = label_file, row.names = FALSE)
                         # --- 保存未对齐的原始聚类标签（非常关键！）
                         original_label_file <- file.path(output_dir, paste0("OriginalLabels_BCC_", safe_combo_name, ".csv"))
                         write.csv(
                            data.frame(Sample = common_samples,
                                       Cluster = labels
                            ), file = original_label_file, row.names = FALSE)
                         
                         # 评估指标
                         nmi_val  <- NMI(Tlabel, aligned_labels)
                         ari_val  <- ARI(Tlabel, aligned_labels)
                         cm       <- confusionMatrix(
                            factor(aligned_labels, levels = sort(unique(Tlabel))),
                            factor(Tlabel,       levels = sort(unique(Tlabel))))
                         acc_val  <- cm$overall["Accuracy"]
                         f1_bycls <- suppressWarnings(cm$byClass[ , "F1"])
                         mac_f1   <- if (C > 2) mean(f1_bycls, na.rm = TRUE) else f1_bycls["F1"]
                         nc       <- length(unique(aligned_labels))
                         
                         # 返回一行结果
                         data.frame(
                            Combination = combo_name,
                            NMI         = unname(nmi_val),
                            ARI         = unname(ari_val),
                            Accuracy    = unname(acc_val),
                            MacroF1     = unname(mac_f1),
                            LabelFile   = label_file,
                            NumClusters = nc,
                            stringsAsFactors = FALSE
                         )
                         
                      }, error = function(e) {
                         warning("Error in ", combo_name, " : ", e$message)
                         return(NULL)
                      })
                      
                      out
                   }

##### 6. 关闭并行 #####
stopCluster(cl)
cat("Parallel cluster stopped.\n")

##### 7. 保存汇总结果 #####
final_csv <- "/data1/zxr/BRCA/Omics_Combination_Evaluation_BCC.csv"
final_rds <- "/data1/zxr/BRCA/output_BCC.Rdata"
write.csv(results, final_csv, row.names = FALSE)
save(results, file = final_rds)
cat("All done! Results saved to:\n", final_csv, "\n", final_rds, "\n")









