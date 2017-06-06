# clean up your data at first

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


#find the outlier
sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut, set the h to difine the place you want to cut
abline(h = 8e+04, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 8e+04, minSize = 10)
table(clust)


# clust 1 contains the samples we want to keep.

pdf(file='sample clustering.pdf')
dev.off()

# cut the outlier
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)




