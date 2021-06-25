# 加载 R 包
library(WGCNA)

# 导入数据
load(file = "wgcna-network.Rdata")

#f_trait <- "MCI_S-MCI_AD_traits.txt"
#f_trait <- "NC-MCI_S_traits.txt"
f_trait <- "traits.txt"

#=====================================================================================
#
#  5. 导入性状数据
#
#=====================================================================================
{
  datTraits <- read.delim(f_trait, row.names=1, 
                          quote="")
  # 确保 datTraits 与 datExpr 的 rownames 顺序一致
  rownames(datTraits)
  traitRows = match(rownames(datExpr), rownames(datTraits))
  datTraits = datTraits[traitRows, ]
}

#=====================================================================================
#
#  6. 模块与性状关系
## 6.1 计算模块特征向量（moduleEigengenes， MEs）
## 6.2 计算模块（MEs）与性状之间的相关性
## 6.3 相关性 heatmap
#=====================================================================================
{
  # 6.0 获得样本数目和基因数目
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  
  # 6.1 计算模块特征向量 MEs
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  
  # 6.2 计算模块（MEs）与性状之间的相关性
  moduleTraitCor = cor(MEs, datTraits, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

  # 6.3 相关性 heatmap
  sizeGrWindow(10,6)
  ## 连接相关性和 pvalue
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  
  
  ## heatmap 画图
  pdf(file = "Module-trait_relationships.pdf",height = 4,width = 12)
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = t(moduleTraitCor),
                 yLabels = names(datTraits),
                 xLabels = names(MEs),
                 xSymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = t(textMatrix),
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  }

#=====================================================================================
#
#  6. 基因与性状关系（GS）& 基因与模块关系（MM）
## 6.1 计算 module membership (MM): 基因（TMP）与模块（MEs）的相关性
## 6.2 计算 Gene Significance (GS): 基因（TMP）与性状的相关性
#   
#=====================================================================================
{
  modNames = substring(names(MEs), 3)
  traitNames = names(datTraits)
  
  # 计算 MM
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                            nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  
  # 计算 GS
  geneTraitSignificance = as.data.frame(cor(datExpr, datTraits, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                            nSamples))
  names(geneTraitSignificance) = paste("GS.", traitNames, sep="");
  names(GSPvalue) = paste("p.GS.", traitNames, sep="");

  geneInfo<-cbind(geneModuleMembership, MMPvalue, geneTraitSignificance, GSPvalue)
  write.table(geneInfo, file = "geneInfo2.txt", 
              sep = "\t", 
              quote = F)
}

save(datExpr, datTraits, MEs, geneModuleMembership, geneTraitSignificance, 
     file = "trait_analysis.RData")
