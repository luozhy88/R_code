# 通过距离判断PCA是否存在差异显著

# 加载必要的包
library(vegan)
library(ggplot2)
library(factoextra)

# 创建示例数据
set.seed(123)
# 创建20个样本，10个特征的矩阵
data_matrix <- matrix(rnorm(200), nrow=20)
rownames(data_matrix) <- paste0("S", 1:20)
colnames(data_matrix) <- paste0("Feature", 1:10)

# 创建分组信息
groups <- factor(rep(c("A","B","C"), length.out=20))
metadata <- data.frame(Group = groups)
rownames(metadata) <- rownames(data_matrix)

#########################PCA##############################


# 加载必要的包
library(ggplot2)
library(factoextra)

# 对数据进行PCA分析
pca_result <- prcomp(data_matrix, scale. = TRUE)

# 使用ggplot2创建更美观的PCA图
pca_data <- as.data.frame(pca_result$x)
pca_data$Group <- groups

# 绘制PCA图
ggplot(pca_data, aes(x=PC1, y=PC2, color=Group)) +
  geom_point(size=3) +
  stat_ellipse(type="t", level=0.95) +
  theme_bw() +
  labs(title="PCA分析结果",
       x=paste0("PC1 (", 
                round(summary(pca_result)$importance[2,1]*100, 1), 
                "%)"),
       y=paste0("PC2 (",
                round(summary(pca_result)$importance[2,2]*100, 1),
                "%)")) +
  scale_color_manual(values=c("A"="#E41A1C", "B"="#377EB8", "C"="#4DAF4A")) +
  theme(legend.position="right",
        plot.title = element_text(hjust = 0.5))

# 添加置信椭圆
# 输出主成分贡献率
summary(pca_result)

########################判断显著#############################
# 计算Bray-Curtis距离矩阵
dist_matrix <- vegdist(data_matrix, method="bray")

# 执行Adonis检验
adonis_result <- adonis2(dist_matrix ~ Group, 
                         data = metadata, 
                         permutations = 999)
print(adonis_result)
adonis_result$`Pr(>F)`

