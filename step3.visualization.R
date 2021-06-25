# 加载 R 包
library(WGCNA)

# 导入数据
load(file = "wgcna-network.Rdata")
load(file = "trait_analysis.RData")
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#=====================================================================================
#
#  7. 基因共表达网络可视化
#   
#=====================================================================================
{
  ## 7.1. 使用所有基因绘制
  # 计算 dissTOM, dissTOM = 1 - TOM
  dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate);
  
  # 取 7 次方，仅为展示更显著
  plotTOM = dissTOM^7
  
  # 绘图
  diag(plotTOM) = NA;
  sizeGrWindow(9,9)
  geneTree = net$dendrograms[[1]]
  TOMplot(plotTOM, geneTree, moduleColors,col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'), main = "Network heatmap plot, all genes")
  
  ## 7.2 随机选择 400 条绘图
  nSelect = 400
  set.seed(10)
  select = sample(nGenes, size = nSelect)
  selectTOM = dissTOM[select, select]
  # 重新聚类
  selectTree = hclust(as.dist(selectTOM), method = "average")
  selectColors = moduleColors[select]
  # Open a graphical window
  sizeGrWindow(9,9)
  plotDiss = selectTOM^7
  # 绘图
  diag(plotDiss) = NA;
  TOMplot(plotDiss, selectTree, selectColors,col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'), main = "Network heatmap plot, selected genes")
}

#=====================================================================================
#
#  8. 模块特征向量网络可视化
#   
#=====================================================================================
{
  # 重新计算MEs
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  # Plot the relationships among the eigengenes and the trait
  sizeGrWindow(5,7.5)
  par(cex = 0.9)
  plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                        = 90)
  
  # Plot the dendrogram
  sizeGrWindow(6,6)
  par(cex = 1.0)
  plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                        plotHeatmaps = FALSE)
  # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  par(cex = 1.0)
  plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                        plotDendrograms = FALSE, xLabelsAngle = 90)
}

#=====================================================================================
#
#  9. 生成 Cytoscape 的输入文件
#   
#=====================================================================================
{
  # 重新计算 TOM
  TOM = TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate);
  # 要可视化的模块
  modules = c("cyan")
  # 要可视化的基因
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule]
  # 候选基因的 TOM
  modTOM = TOM[inModule, inModule]
  
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 nodeAttr = moduleColors[inModule])
}

