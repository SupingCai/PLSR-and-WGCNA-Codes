# 加载 R 包
#install.packages("WGCNA")
library(WGCNA)


f_exp <- "gene_regional_expression.txt"

#=====================================================================================
#
#  1.导入表达数据
#
#=====================================================================================
{
  # 不要将文件中的字符串转换为因子
  options(stringsAsFactors = FALSE)
  
  # 读取表达数据
  exp0 <- read.delim(f_exp, row.names=1)
  
  # 过滤掉在所有样品中表达之和小于10 的基因
  # exp <- exp0[rowSums(exp0) > 10, ]
  
  # 转置
  datExpr = t(exp0)
  dim(datExpr)
  rownames(datExpr)
  
  # 保存数据
  write.table(datExpr, file = "datExpr.txt",
              sep = "\t", quote = F, row.names = T)
}

#=====================================================================================
#
#  2.寻找最佳 β 值
#
#=====================================================================================
{
  # 开启多线程模式
  enableWGCNAThreads(nThreads = 6)
  
  # 通过对 power 的多次迭代，确定最佳 power
  powers <- c(1:20)
  sft = pickSoftThreshold(datExpr, 
                          powerVector = powers, 
                          verbose = 5,
                          networkType = "signed"
                          )
  # 计算出来的最佳 β 存放在
  sft$powerEstimate
  
  # 画图
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2))
  cex1 = 0.9
  #  R2 ~ soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.8,col="red")
  # Mean connectivity ~ soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  sft$powerEstimate = 10
}

#=====================================================================================
#
#  3. 构建网络
##  3.1. 计算相关系数
##  3.2. 计算邻接矩阵
##  3.3. 计算 TOM 矩阵
##  3.4. 聚类并划分模块
##  3.5. 合并相似模块
#=====================================================================================
{
  net = blockwiseModules(
    # 0.输入数据
    datExpr, 
    
    # 1. 计算相关系数
    corType = "pearson", # 相关系数算法，pearson|bicor
    
    # 2. 计算邻接矩阵
    power = sft$powerEstimate, # 前面得到的 soft power
    networkType = "unsigned", # unsigned | signed | signed hybrid
    
    # 3. 计算 TOM 矩阵
    TOMType = "unsigned", # none | unsigned | signed
    saveTOMs = TRUE, # 是否保存
    saveTOMFileBase = "blockwiseTOM", # 保存文件前缀
    
    # 4. 聚类并划分模块
    deepSplit = 2, # 0|1|2|3|4, 值越大得到的模块就越多越小
    minModuleSize = 50,
    
    # 5. 合并相似模块
    ## 5.1 计算模块特征向量（module eigengenes， MEs），即第一个主成分（1st PC）
    ## 5.2 计算 MEs 与 datTrait 之间的相关性
    ## 5.3 对距离小于 mergeCutHeight （1-cor）的模块进行合并
    mergeCutHeight = 0.15, 

    # 其他参数
    numericLabels = TRUE, # 以数字命名模块
    nThreads = 0 # 0 则使用所有可用线程
    )
  # 查看每个模块包含基因数目
  table(net$colors) 
}


#=====================================================================================
#
#  4. 可视化
#
#=====================================================================================
{
  # 模块名称修改为颜色
  moduleColors = labels2colors(net$colors)
  # 同时绘制聚类图和模块颜色
  pdf(file = "plotDendroAndColors.pdf")
  plotDendroAndColors(
    net$dendrograms[[1]], 
    moduleColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE, 
    addGuide = TRUE)
  dev.off()
}

save(datExpr, sft, net, moduleColors, file = "wgcna-network.Rdata")
