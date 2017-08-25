# WGCNA
# This is my own Rscript for WGCNA in bioinformatic analysis
WGCNA provides us based on topological correlation analysis for us and fclust for those modules. 
While there too many genes in the networks, our program might burning all memory.
Hence, it provides two kinds modules for this analysis.
First method is applied when there are not too much genes to analysis and can provides more accurate results, and its result can be easily export to cyberscape.
However, while genes are numerous, by separating them into several small group according to their adjancency, we can analysis them one by one and integrate them togather.
