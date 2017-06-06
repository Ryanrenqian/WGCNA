file='all_sample_gene.fpkm.xls'
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
data=read.table(file,header=T)
rownames(data)=data[,1]
datExpr=data[,-1]
datExpr=as.data.frame(t(datExpr))
rm(data)
#check the value
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

#
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}


powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

pdf(file='gene.pdf')
dev.off()

softPower=9
# now here is the solution
bwnet = blockwiseModules(datExpr, maxBlockSize = 10000,
                         power = softPower, TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "gene",
                         verbose = 3)

# Load the results of single-block analysis
load(file = "gene-block.1.RData");
load(file="gene-block.2.RData")
# Relabel blockwise modules
#but i can not find moduleLables
moduleLabels=bwnet$colors
bwLabels = matchLabels(bwnet$colors, moduleLabels);
# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)
#look your labels
table(bwLabels)


sizeGrWindow(6,6)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# congrulations!!! now you can export data to cytoscape


# notice that we have make several gene blocks, which means we can only convert them to several cytoscape files one by one!
#set the module you want to plot# first you need to choose two module to plot
modules = c("brown", "turquoise")
# Here I choose the first block to plot. 
includeColnames = is.finite(match(bwModuleColors[bwnet$blockGenes[[1]]], modules))
nodeNames = names(datExpr[bwnet$blockGenes[[1]]])
exportNetworkToCytoscape(TOM,threshold = 0.5,
                         edgeFile = paste("gene-FPKM-edges-",paste(modules,collapse="-"),".txt",sep=""),
                         nodeFile = paste("gene-FPKM-node-",paste(modules,collapse="-"),".txt",sep=""),
                         weighted = TRUE,nodeNames = nodeNames,includeColNames = includeColnames)




