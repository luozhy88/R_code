

##############################
#A.rna列为样本，行为特征
#此处方法需要考虑分组信息进行批次矫正
# 不同批次之间的特征数目是否影响运行需要测试
# <!-- 根据代码分析，我来总结ARSyNbac方法去除批次效应的主要流程： -->
# 
# <!-- ## 核心流程 -->
# 
# <!-- **数据预处理**： -->
# <!-- - 将所有批次的数据合并成一个大矩阵X -->
# <!-- - 构建实验设计矩阵designA和批次设计矩阵designB -->
# <!-- - 对数据进行中心化处理，去除整体平均值的影响 -->
# 
# <!-- **模型分解**： -->
# <!-- 主要通过ASCA (ANOVA Simultaneous Component Analysis)方法将数据分解为以下几个部分： -->
# <!-- - Model.a：实验设计效应 -->
# <!-- - Model.b：批次效应（当Interaction=FALSE时） -->
# <!-- - Model.bab：批次与实验设计的交互效应（当Interaction=TRUE时） -->
# <!-- - Model.res：残差/噪声部分 -->
# 
# <!-- ## 校正原理 -->
# 
# <!-- **批次效应去除**： -->
# <!-- - 从原始数据中减去Model.b（或Model.bab）的影响 -->
# <!-- - 如果filterNoise=TRUE，还会额外去除噪声部分（Model.res） -->
# <!-- - 保留实验设计效应（Model.a）和生物学相关变异 -->
# 
# <!-- **参数控制**： -->
# <!-- - beta参数：控制主成分数量的选择 -->
# <!-- - Variability参数：控制要解释的总体变异比例 -->
# <!-- - batchEstimation：控制是否进行批次效应估计 -->
# <!-- - Interaction：控制是否考虑交互效应 -->
# 
# <!-- ## 技术细节 -->
# 
# <!-- 该方法使用ASCA算法进行数据分解： -->
# <!-- - 利用主成分分析（PCA）降维 -->
# <!-- - 通过方差分解识别不同来源的变异 -->
# <!-- - 分别建模并去除非期望的变异源[1] -->


##############################
rm(list=ls())
library(MultiBaC)
data('multiyeast')

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
# custom_pch <- c(19,19,19,1,1,1,
#                 19,19,1,1,
#                 19,19,1,1)

## ----arsyn4, message = FALSE, fig.cap = "Both modes together. Left: Explained variance plot. Default plot when showplot = TRUE. It represents the percentage of variability in the residuals (y axis) explained by a model with a given number of principal components (x axis). The number of selected components in the final model is indicated with a triangle symbol, and computed to explain beta times the average variability of the residuals. Right: PCA plot of corrected gene expression data."----
par(mfrow = c(1,2))
arsyn_4 <- ARSyNbac(data_RNA, modelName = "RNA", beta = 0.5, 
                    batchEstimation = TRUE, filterNoise = TRUE,
                    Interaction = FALSE,
                    showplot=T,
                    Variability = 0.95)
# trace(ARSyNbac,edit=T)
plot(arsyn_4, typeP="pca.cor", bty = "L",
     # pch = custom_pch,
     cex = 3, col.per.group = custom_col,
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
# knitr::include_graphics("designScheme.png", dpi = 30)
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
