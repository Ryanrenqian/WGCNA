#!/usr/bin/env Rscript

# memory use caution: 20000 genes --> 16G 30000 genes-->32g 8000-10000-->4G
getdata<-function(file){
	data=read.table(file,header=T,row.names = 1)
	return (data)
}
cleandata<-function(datExpr){
	gsg = goodSamplesGenes(datExpr, verbose = 3)
	if (!gsg$allOK)
	{
		# Optionally, print the gene and sample names that were removed:
		if (sum(!gsg$goodGenes)>0)
			printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
		if (sum(!gsg$goodSamples)>0)
			printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
		# Remove the offending genes and samples from the data:
		datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
		rm(gsg)
	}
	return(datExpr)
}
ChoosePower<-function(datExpr,level=0.9){
	powers = c(c(1:15), seq(from = 16, to=20, by=2))
	sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
	pdf(file="PowerSelect.pdf")
	sizeGrWindow(9, 5)
	par(mfrow = c(1,2))
	cex1 = 0.9
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",
	     ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
	abline(h=level,col="red")
	plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
	dev.off()
	for(i in 1:length(sft$fitIndices[,1])){
		if (-sign(sft$fitIndices[i,3])*sft$fitIndices[i,2]>=level){
			softPower=sft$fitIndices[i,1]
			break
		}
	}
	return(softPower)
}
# lack of memory
solution2=function(datExpr,softPower,name){
	bwnet = blockwiseModules(datExpr, maxBlockSize = 10000,
				 power = softPower, TOMType = "unsigned", minModuleSize = 30,
				 reassignThreshold = 0, mergeCutHeight = 0.25,
				 numericLabels = TRUE,
				 saveTOMs = TRUE,
				 saveTOMFileBase = name,
				 verbose = 3)
	return(bwnet)
}
# enough memory
solution1=function(datExpr,softPower,name,minMod=30,MEDissThres=0.75){
	adjacency = adjacency(datExpr, power = softPower)
	TOM = TOMsimilarity(adjacency);
	dissTOM = 1-TOM
	geneTree = hclust(as.dist(dissTOM), method = "average");
	pdf(file='Genecluster.pdf')
	plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
	     labels = FALSE, hang = 0.04)
	dev.off()
	minModuleSize = minMod
	dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
				    deepSplit = 2, pamRespectsDendro = FALSE,
				    minClusterSize = minModuleSize);
	table(dynamicMods)
	dynamicColors = labels2colors(dynamicMods)
	table(dynamicColors)
	pdf(file='DynamicGenecluster.pdf')
	plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
			    dendroLabels = FALSE, hang = 0.03,
			    addGuide = TRUE, guideHang = 0.05,
			    main = "Gene dendrogram and module colors")
	dev.off()
	# ploting heat map

	plotTOM = dissTOM^softPower;
	# Set diagonal to NA for a nicer plot
	diag(plotTOM) = NA;
	# Call the plot function, actually if you want plot all gene the memory-use might spill out
	pdf(file='gene heatmap.pdf')
	sizeGrWindow(9,9)
	TOMplot(plotTOM, geneTree, dynamicColors,dynamicColors, main = "Network heatmap plot, all genes")
}
solution3=function(datExpr,softPower,name){
  net = blockwiseModules(datExpr, power = softPower,
                         TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = FALSE,
                         verbose = 3)
  sizeGrWindow(12, 9)
  mergedColors = labels2colors(net$colors)
  pdf(file="Genecluster.pdf")
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  TOM = TOMsimilarityFromExpr(datExpr, power = softPower)
  nodeNames = colnames(datExpr)
  for (modules in unique(mergedColors)){
    includeColnames = is.finite(match(mergedColors, modules))
    modNames=nodeNames[includeColnames]
    modTOM = TOM[includeColnames,includeColnames];
    dimnames(modTOM) = list(modNames, modNames)
    cyt=exportNetworkToCytoscape(modTOM,threshold = 0.5,
                                 edgeFile = paste(name,"-edges-",paste(modules,collapse="-"),".txt",sep=""),
                                 nodeFile = paste(name,"-node-",paste(modules,collapse="-"),".txt",sep=""),
                                 weighted = TRUE,nodeNames = modNames,nodeAttr = mergedColors[includeColnames])
  }
  
}



library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
library(getopt)
opt=getopt(matrix(c(
		    'help', 'h', 0, "logical",
		    'input', 'i', 1, "character",
		    'sample','s',1,"character",
		    'lilelihood','l',1,'float',
		    'threshod','t',1,'float'),
		  ncol=4,byrow=TRUE))
if( !(is.null(opt$help)) || is.null(opt$input)||is.null(opt$sample))
{
	cat( paste0(
		    "Usage: Rscript test_edgeR.R [arguments]\n",
		    "-h     --help          Show this help.\n",
		    "-i     --input          Input gene reads_counts.\n",
		    "-s             --sample       output name.\n",
		    "\n",
		    "\n"
		    ))
	quit()
}

file=opt$input
name=opt$sample
data=getdata(file)
datExpr=cleandata(t(data))
softPower=ChoosePower(datExpr)
solution3(datExpr,softPower,name)
