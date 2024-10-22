


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
