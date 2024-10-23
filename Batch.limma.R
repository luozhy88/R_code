
#lima可以不放入group信息做批次矫正，批次中如果第一批次需要做ref，需要factor按照levels顺序进行
# 
# # 使用limma包
# library(limma)
# batch_corrected_data <- removeBatchEffect(data_matrix, batch=batch_info)



ngenes <- 10
nsamples <- 8
y <- matrix(rnorm(ngenes*nsamples),ngenes,nsamples)
group <- factor(c("A","A","A","A","B","B","B","B"))
batch <- factor(c(1,1,2,2,1,1,2,2))
colnames(y) <- paste(group,batch,sep=".")
y[,batch==2] <- y[,batch==2] + 5
y[,group=="B"] <- y[,group=="B"] + 1

# trace(removeBatchEffect,edit=TRUE)
y.corrected <- removeBatchEffect(y, batch=batch)


oldpar <- par(mfrow=c(1,2))
plotMDS(y,main="Original")
plotMDS(y.corrected,main="Batch corrected")
par(oldpar)



