# ComBat采用经验贝叶斯方法进行批次效应校正，主要包含以下步骤：
# 
# 标准化处理
# 
# 首先对数据进行标准化，使基因表达值在不同批次间具有可比性
# 计算每个基因的均值和方差，进行标准化转换1
# 
# 批次效应估计
# 
# 使用线性模型估计每个基因在不同批次中的位置参数(γ)和尺度参数(δ)1
# 
# 位置参数反映批次间的系统性偏差
# 尺度参数反映批次内的变异程度
# 
# 优势特点
# 
# 适用性广
# 
# 适用于小样本量数据
# 可以处理多个批次3
# 
# 稳健性
# 
# 采用经验贝叶斯框架避免过度校正
# 保持生物学差异的同时去除技术性批次效应3
# 
# 使用建议
# 
# 推荐使用随机分布的样本设计
# 不建议使用非随机的样本分布方式

# 1. 创建测试数据
set.seed(123)

# 创建表达矩阵
n_genes <- 1000  
n_samples <- 30  
expression_matrix <- matrix(rnorm(n_genes * n_samples, mean=10, sd=2), 
                            nrow=n_genes, 
                            ncol=n_samples)

# 添加批次效应
batch <- rep(c(1,2,3), each=10)  
batch_effect <- c(rep(2, 10), rep(-1, 10), rep(0.5, 10))
for(i in 1:n_samples) {
  expression_matrix[,i] <- expression_matrix[,i] + batch_effect[i]
}

# 设置行名和列名
rownames(expression_matrix) <- paste0("gene", 1:n_genes)
colnames(expression_matrix) <- paste0("sample", 1:n_samples)

# 2. 加载必要的包
library(sva)
library(ggplot2)

# 3. 运行ComBat校正
combat_data <- ComBat(dat = expression_matrix, 
                      batch = batch,
                      par.prior = TRUE,
                      mean.only = FALSE)

# 4. 可视化比较
# 校正前的PCA
pca_before <- prcomp(t(expression_matrix))
pca_df_before <- data.frame(PC1 = pca_before$x[,1],
                            PC2 = pca_before$x[,2],
                            Batch = factor(batch))

# 校正后的PCA
pca_after <- prcomp(t(combat_data))
pca_df_after <- data.frame(PC1 = pca_after$x[,1],
                           PC2 = pca_after$x[,2],
                           Batch = factor(batch))

# 绘制校正前的PCA图
p1 <- ggplot(pca_df_before, aes(x=PC1, y=PC2, color=Batch)) +
  geom_point(size=3) +
  theme_bw() +
  ggtitle("before PCA") +
  theme(plot.title = element_text(hjust = 0.5))

# 绘制校正后的PCA图
p2 <- ggplot(pca_df_after, aes(x=PC1, y=PC2, color=Batch)) +
  geom_point(size=3) +
  theme_bw() +
  ggtitle("after PCA") +
  theme(plot.title = element_text(hjust = 0.5))

# 并排显示两个图
library(gridExtra)
grid.arrange(p1, p2, ncol=2)
