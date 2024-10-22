


library(MultiBaC)

# 首先,检查每个批次的样本数
n_samples_A <- ncol(A.rna)
n_samples_B <- ncol(B.rna)
n_samples_C <- ncol(C.rna)

# 创建新的实验设计
experimentalDesign <- list(
  "A" = rep(c("X1", "X2"), length.out = n_samples_A),
  "B" = rep(c("X1", "X2"), length.out = n_samples_B),
  "C" = rep(c("X1", "X2"), length.out = n_samples_C)
)

# 重新创建mbac对象
data_RNA <- createMbac(
  inputOmics = list(A.rna, B.rna, C.rna), 
  batchFactor = c("A", "B", "C"),
  experimentalDesign = experimentalDesign,
  omicNames = "RNA"
)

custom_col <- c("brown2", "dodgerblue", "forestgreen")
custom_pch <- c(19,19,19,1,1,1,
                19,19,1,1,
                19,19,1,1)

# 现在尝试运行ARSyNbac
arsyn_result <- ARSyNbac(
  data_RNA, 
  modelName = "RNA", 
  Variability = 0.95, 
  batchEstimation = TRUE, 
  filterNoise = FALSE, 
  Interaction = FALSE
)

## ----arsyn2, message = FALSE, fig.cap = "Batch correction when Interaction = TRUE. Left: Explained variance plot. Default plot when showplot = TRUE. It represents the number of components (x axis) needed to explain a certain variability (y axis) of the effect of interest (batch effect). The number of components selected in the model is indicated with a triangle symbol. Gray dashed line indicates the threshold given by the Variability argument (in percentage). Right: PCA plot of corrected gene expression data considering the interaction batch-condition."----
par(mfrow = c(1,2))
arsyn_2 <- ARSyNbac(data_RNA, modelName = "RNA", Variability = 0.95, 
                    batchEstimation = TRUE, filterNoise = FALSE, Interaction = TRUE)
plot(arsyn_2, typeP="pca.cor", bty = "L",
     pch = custom_pch, cex = 3, col.per.group = custom_col,
     legend.text = c("Color: Batch", names(data_RNA$ListOfBatches),
                     "Fill: Cond.", unique(colData(data_RNA$ListOfBatches$A)$tfactor)),
     args.legend = list("x" = "topright",
                        "pch" = c(NA, 15, 15, 15, 
                                  NA, 15, 0),
                        "col" = c(NA, "brown2", "dodgerblue", "forestgreen",
                                  NA, 1, 1),
                        "bty" = "n",
                        "cex" = 1.2))


## ----arsyn3, message = FALSE, fig.cap = "Noise reduction mode. Left: Explained variance plot. Default plot when showplot = TRUE. It represents the percentage of variability in the residuals (y axis) explained by a model with a given number of principal components (x axis). The number of selected components in the final model is indicated with a triangle symbol, and computed to explain beta times the average variability of the residuals. Right: PCA plot of corrected gene expression data."----
par(mfrow = c(1,2))
arsyn_3 <- ARSyNbac(data_RNA, modelName = "RNA", beta = 0.5, 
                    batchEstimation = FALSE, filterNoise = TRUE)
plot(arsyn_3, typeP="pca.cor", bty = "L",
     pch = custom_pch, cex = 3, col.per.group = custom_col,
     legend.text = c("Color: Batch", names(data_RNA$ListOfBatches),
                     "Fill: Cond.", unique(colData(data_RNA$ListOfBatches$A)$tfactor)),
     args.legend = list("x" = "topright",
                        "pch" = c(NA, 15, 15, 15, 
                                  NA, 15, 0),
                        "col" = c(NA, "brown2", "dodgerblue", "forestgreen",
                                  NA, 1, 1),
                        "bty" = "n",
                        "cex" = 1.2))


## ----arsyn4, message = FALSE, fig.cap = "Both modes together. Left: Explained variance plot. Default plot when showplot = TRUE. It represents the percentage of variability in the residuals (y axis) explained by a model with a given number of principal components (x axis). The number of selected components in the final model is indicated with a triangle symbol, and computed to explain beta times the average variability of the residuals. Right: PCA plot of corrected gene expression data."----
par(mfrow = c(1,2))
arsyn_4 <- ARSyNbac(data_RNA, modelName = "RNA", beta = 0.5, 
                    batchEstimation = TRUE, filterNoise = TRUE,
                    Interaction = TRUE,
                    Variability = 0.95)
plot(arsyn_4, typeP="pca.cor", bty = "L",
     pch = custom_pch, cex = 3, col.per.group = custom_col,
     legend.text = c("Color: Batch", names(data_RNA$ListOfBatches),
                     "Fill: Cond.", unique(colData(data_RNA$ListOfBatches$A)$tfactor)),
     args.legend = list("x" = "topright",
                        "pch" = c(NA, 15, 15, 15, 
                                  NA, 15, 0),
                        "col" = c(NA, "brown2", "dodgerblue", "forestgreen",
                                  NA, 1, 1),
                        "bty" = "n",
                        "cex" = 1.2))

## ----design, echo = FALSE, fig.cap = "Scheme of the yeast example data structure."----
knitr::include_graphics("designScheme.png", dpi = 30)
# 获取校正后的数据
rna_data <- assay(arsyn_4$CorrectedData$A, "RNA")

## -----------------------------整理结果------------------------------------------------
# 提取并合并所有批次的RNA数据
all_rna_df <- do.call(cbind, lapply(names(arsyn_4$CorrectedData), function(batch) {
  rna_df <- as.data.frame(assay(arsyn_4$CorrectedData[[batch]], "RNA"))
  colnames(rna_df) <- paste0(batch, "_", colnames(rna_df))
  return(rna_df)
}))

# 添加基因ID
all_rna_df$gene_id <- rownames(all_rna_df)

# 查看结果
head(all_rna_df)
colnames(all_rna_df)


# 查看各批次列名
lapply(names(arsyn_4$CorrectedData), function(batch) {
  grep(paste0("^", batch, "_"), colnames(all_rna_df), value = TRUE)
})

# 残差的平均变异度是指残差(实际观测值与预测值之间的差异)的方差或标准差的平均水平。具体来说:
#   
#   残差定义:
#   残差 = 实际观测值 - 预测值
#   残差方差:
#     残差方差衡量了残差围绕其均值(通常为0)的分散程度。
#   平均变异度计算:
#     计算每个数据点的残差
#   计算这些残差的方差或标准差
#   如果有多个组或批次,可以取这些方差或标准差的平均值

# ARSyNbac算法中beta参数的作用,其含义如下:
#   
# beta是一个数值参数。它用于识别残差中的系统噪声。
# 具体工作原理:
#   计算残差的平均变异度。
#   将平均变异度乘以beta值,得到一个阈值。
#   任何解释变异度超过这个阈值的主成分,都被识别为系统噪声。
#   这个参数主要在"噪声减少模式"下使用,即当批次效应未知时。
#   默认值为2,意味着默认情况下,变异度超过平均值2倍的主成分会被视为系统噪声。
#   调整beta值可以控制噪声移除的程度:
#     增大beta值会导致识别为噪声的主成分数量减少,保留更多的变异。
#     减小beta值会导致更多的主成分被识别为噪声,移除更多的变异。
# 
# 这个参数的设计允许用户根据数据特征和分析需求灵活调整噪声识别的敏感度。它在保留有意义的生物学变异和移除不必要的系统噪声之间取得平衡,是ARSyNbac方法灵活性的重要体现。


  
