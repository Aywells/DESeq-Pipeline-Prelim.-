setwd('C:/Users/aywel/Desktop/Research/Barabasi Lab/Nomix/Code/Data (New)')

library(xopen) # General Library
library(readxl)
library(readr)
library(plyr)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(dplyr)
library(vsn)
library(rgl)

library(GEOquery) # Gene Library
library(Biobase)
library(BiocGenerics)
library(BiocManager)
library(biomaRt)

library(sva) # Batch Effect Library
library(pamr)
library(limma)
library(DESeq2)

data.total <- read.csv("Raw_cl.csv")    # RAW DATA
data.pheno <- data.total[,3:15]    # Phenotypic logicals for diseases
data.counts <- data.total[,19:4653]    # Raw counts for every sample
data.batch <- data.total[,1:2]    # Just sample/batch IDs

# Trial Analysis ---------------------------------------------------------------

# (Need better and more raw_count patients)
raw_counts <- matrix(as.matrix(data.total[,2:602]), ncol = ncol(data.total[,2:602]), dimnames = NULL)
rownams <- (data.total[,1])
colnams <- colnames(data.total)
colnams <- (colnams[2:602])
colnames(raw_counts) <- colnams
rownames(raw_counts) <- rownams
raw_counts <- as.matrix(raw_counts)    # Raw_counts mtx for all samples

col.data <- data.frame(data.pheno[,2]) # (Need more phenotype data)
rownames(col.data) <- colnams    # Essentially phenotypes

DESeq.ds <- DESeqDataSetFromMatrix ( countData = raw_counts,
                                     colData = col.data,
                                     design = ~ Batch )    # DEseq format for raw_counts input

DESeq.ds <- DESeq.ds [ rowSums ( counts ( DESeq.ds ) ) > 0 , ]

DESeq.ds <- estimateSizeFactors ( DESeq.ds )
sizeFactors ( DESeq.ds ) # Size factor normalization

counts.sf_norm <- counts ( DESeq.ds , normalized = TRUE )    # normalized counts

counts.sf_normlog2 <- log2( counts.sf_norm + 1)    # normalized and log2 transformed counts

#DESeq.rlog <- rlog ( DESeq.ds , blind = TRUE )    # super slow varience shrinkage transform
#rlog.norm.counts <- assay ( DESeq.rlog )  

DESeq.vst <- vst(DESeq.ds)    # faster varience shrinkage transform (norm and log2)
vst.norm.counts <- assay(DESeq.vst)

distance.m_vst <- as.dist (1 - cor( vst.norm.counts , method = "pearson" ) )
plot ( hclust ( distance.m_vst, method = "average" ) , # Hierarchical clustering + visual
       labels = colnames ( vst.norm.counts ) ,
       main = " VST Read Counts : Pearson correlation " )

str ( colData ( DESeq.ds )$Batch )
colData ( DESeq.ds )$Batch <- relevel ( colData ( DESeq.ds )$Batch , "Asthma1" )

# True Analysis ---------------------------------------------------------------

DESeq.ds <- DESeq (DESeq.ds) # DESeq analysis

DGE.results <- results ( DESeq.ds , independentFiltering = TRUE , alpha = 0.05)
summary ( DGE.results )

hist ( DGE.results $ pvalue , # P-value histogram for all genes
       col = "grey" , border = "white" , xlab = " " , ylab = " " ,
       main = " Frequencies of P-values " ) 

plotMA( DGE.results , alpha = 0.05 , main = "MA Plot" , # Scatter plot for all genes (significant in red)
        ylim = c( -11 ,30) )

DGE.results.sorted <- DGE.results [ order ( DGE.results$padj ) , ]
DGEgenes <- rownames ( subset ( DGE.results.sorted , padj<0.05) )

upGenes <- subset(DGE.results, log2FoldChange > 0)
# upGenes <- upGenes[,2]

downGenes <- subset(DGE.results, log2FoldChange < 0)
# downGenes <- downGenes[,2]

L1000CDS2(upGenes,downGenes,signature = reverse)

# Heatmaps ---------------------------------------------------------------

hm.mat_DGEgenes <- counts.sf_normlog2 [DGEgenes , ]
aheatmap ( hm.mat_DGEgenes , Rowv = NA , Colv = NA )

aheatmap ( hm.mat_DGEgenes ,
           Rowv = TRUE , Colv = TRUE , # add dendrograms to rows and columns
           distfun = "euclidean" , hclustfun = "average" )

aheatmap ( hm.mat_DGEgenes ,
           Rowv = NA , Colv = NA ,
           distfun = "euclidean" , hclustfun = "average" ,
           scale = "row" ) # values are transformed into distances from the center of the row - specific average : ( actual value - mean of the group ) / standard deviation

# Mean-SD Plot ---------------------------------------------------------------

msd_plot <- meanSdPlot ( raw_counts ,
                         ranks = FALSE ,
                         plot = FALSE )
msd_plot $gg +
  ggtitle ( " Raw Read Counts " ) +
  ylab ( " Standard Deviation " ) +
  xlab ( " Mean " )

# PCA (2-D) ---------------------------------------------------------------

Lcombo.pca <- prcomp(t(vst.norm.counts))

PCAscores <- Lcombo.pca$x
PCAloadings <- Lcombo.pca$rotation
PCAcolors <- c("green","violet","blue", "pink", "yellow", "red", 'orange', "cyan", "magenta", "grey0")

xmin <- -100
xmax <- 100
ymin <- -200
ymax <- 200

Part1PC <- 2
Part2PC2 <- 3

plot(PCAscores[1:39,Part1PC:Part2PC2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(1)],   
     bg=PCAcolors[c(1)],
     cex=1.0,         
     main="Scores"     
)
par(new=T)
plot(PCAscores[40:93,Part1PC:Part2PC2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(2)],   
     bg=PCAcolors[c(2)],
     cex=1.0,         
     main="Scores"     
)
par(new=T)
plot(PCAscores[94:140,Part1PC:Part2PC2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(3)],   
     bg=PCAcolors[c(3)],
     cex=1.0,         
     main="Scores"     
)
par(new=T)
plot(PCAscores[141:198,Part1PC:Part2PC2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(4)],   
     bg=PCAcolors[c(4)],
     cex=1.0,         
     main="Scores"     
)
par(new=T)
plot(PCAscores[199:292,Part1PC:Part2PC2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(5)],   
     bg=PCAcolors[c(5)],
     cex=1.0,         
     main="Scores"
     
)
par(new=T)
plot(PCAscores[293:494,Part1PC:Part2PC2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(6)],   
     bg=PCAcolors[c(6)],
     cex=1.0,         
     main="Scores"
)
par(new=T)
plot(PCAscores[495:500,Part1PC:Part2PC2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(7)],   
     bg=PCAcolors[c(7)],
     cex=1.0,        
     main="Scores"
     
)
par(new=T)
plot(PCAscores[501:531,Part1PC:Part2PC2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(8)],   
     bg=PCAcolors[c(8)],
     cex=1.0,         
     main="Scores"
)
par(new=T)
plot(PCAscores[532:594,Part1PC:Part2PC2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(9)],   
     bg=PCAcolors[c(9)],
     cex=1.0,         
     main="Scores",
)
par(new=T)
plot(PCAscores[532:594,Part1PC:Part2PC2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(9)],   
     bg=PCAcolors[c(9)],
     cex=1.0,         
     main="Scores",
)
par(new=T)
plot(PCAscores[595:601,Part1PC:Part2PC2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(10)],   
     bg=PCAcolors[c(10)],
     cex=1.0,         
     main="Scores",
     legend("topright",
            c("L8","L3","L4","L14","L1","L6","L9","L7","L15","Asthma"),
            fill=c("green","violet","blue", "pink", "yellow", "red", 'orange', "cyan", "magenta", "grey0"),cex = 0.5)
)

plot(PCAloadings[,Part1PC:Part2PC2], 
     pch=21,             
     bg="black",        
     cex=1,             
     main="Loadings"  
)
text(PCAloadings[,1:2],           
     labels=rownames(PCAloadings))

# PCA (3-D) ---------------------------------------------------------------

Lcombo.pca <- prcomp(t(raw_counts))

PCAscores <- Lcombo.pca$x
PCAloadings <- Lcombo.pca$rotation
PCAcolors <- c("green","violet","blue", "pink", "yellow", "red", 'orange', "cyan", "magenta", "grey0")

xlab <- "PC1"
ylab <- "PC2"
zlab <- "PC3"

plot3d(PCAscores[1:39,1:3],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
       xlab = xlab, ylab = ylab, zlab = zlab,
       col=PCAcolors[c(1)],   
       size = 5,
       main = "", sub = "", ann = FALSE, axes = FALSE,
)
box3d()
par(new=T)
plot3d(PCAscores[40:93,1:3],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
       xlab = xlab, ylab = ylab, zlab = zlab,
       col=PCAcolors[c(2)],   
       size = 5,
       main = "", sub = "", ann = FALSE, axes = FALSE,
)
box3d()
par(new=T)
plot3d(PCAscores[94:140,1:3],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
       xlab = xlab, ylab = ylab, zlab = zlab,
       col=PCAcolors[c(3)],   
       size = 5,
       main = "", sub = "", ann = FALSE, axes = FALSE,
)
box3d()
par(new=T)
plot3d(PCAscores[141:198,1:3],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
       xlab = xlab, ylab = ylab, zlab = zlab,
       col=PCAcolors[c(4)],   
       size = 5,
       main = "", sub = "", ann = FALSE, axes = FALSE,
)
box3d()
par(new=T)
plot3d(PCAscores[199:292,1:3],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
       xlab = xlab, ylab = ylab, zlab = zlab,
       col=PCAcolors[c(5)],   
       size = 5,
       main = "", sub = "", ann = FALSE, axes = FALSE,
)
box3d()
par(new=T)
plot3d(PCAscores[293:494,1:3],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
       xlab = xlab, ylab = ylab, zlab = zlab,
       col=PCAcolors[c(6)],   
       size = 5,
       main = "", sub = "", ann = FALSE, axes = FALSE,
)
box3d()
par(new=T)
plot3d(PCAscores[495:500,1:3],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
       xlab = xlab, ylab = ylab, zlab = zlab,
       col=PCAcolors[c(7)],   
       size = 5,
       main = "", sub = "", ann = FALSE, axes = FALSE,
)
box3d()
par(new=T)
plot3d(PCAscores[501:531,1:3],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
       xlab = xlab, ylab = ylab, zlab = zlab,
       col=PCAcolors[c(8)],   
       size = 5,
       main = "", sub = "", ann = FALSE, axes = FALSE,
)
box3d()
par(new=T)
plot3d(PCAscores[532:594,1:3],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
       xlab = xlab, ylab = ylab, zlab = zlab,
       col=PCAcolors[c(9)],   
       size = 5,
       main = "", sub = "", ann = FALSE, axes = FALSE,
)
box3d()
par(new=T)
plot3d(PCAscores[595:601,1:3],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
       xlab = xlab, ylab = ylab, zlab = zlab,
       col=PCAcolors[c(10)],   
       size = 5,
       main = "", sub = "", ann = FALSE, axes = FALSE,
)
box3d()


# Function (L1000 Query) ---------------------------------------------------------------
L1000CDS2 = function(upGenes= NULL, downGenes = NULL, signature = NULL,
                     mimic = TRUE, combination = TRUE,dbVersion = 'latest'){
  baseLink = 'http://amp.pharm.mssm.edu/L1000CDS2/query'
  if(is.null(signature)){
    assertthat::assert_that(!is.null(upGenes) & !is.null(downGenes))
    payload = list(data = list(upGenes = upGenes,
                               dnGenes = downGenes),
                   config = list(aggravate = mimic,
                                 searchMethod = 'geneSet',
                                 share = FALSE,
                                 combination = combination,
                                 `db-version` = dbVersion))
    
    response = httr::POST(baseLink,body = payload,encode = 'json')
  } else{
    payload = list(data  = list(genes = names(signature),vals = unlist(unlist(signature))),
                   config = list(aggravate = mimic,
                                 searchMethod = 'CD',
                                 share = FALSE,
                                 combination = combination,
                                 `db-version` = dbVersion))
    response = httr::POST(baseLink,body = payload,encode = 'json')
    
  }
  
  response <- jsonlite::fromJSON(httr::content(response, 'text'))
  return(response)
}
