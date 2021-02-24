library(xopen) # General Library
library(readxl)
library(readr)
library(plyr)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(dplyr)

library(GEOquery) # Gene Library
library(Biobase)
library(BiocGenerics)
library(BiocManager)
library(biomaRt)

library(sva) #Batch Effect Library
library(pamr)
library(limma)
library(DESeq2)

data.combo <- read.table("Control_Group2.txt")    # RAW DATA
expData <- read_excel("Express_tot.xlsx")    # Raw data + sample number
phenoData <- read_excel("Pheno.xlsx")    # Batch numbers for each sample number
batchno <- as.data.frame(phenoData[,2])    # Just batch numbers


# --------------------------------------------------------------- Normalization
vstd <- vst(as.matrix(data.combo[,2:(595)])) # Quicker version of below (same variance)
vsd <- varianceStabilizingTransformation(as.matrix(data.combo[,2:595]))
rld <- rlog(as.matrix(data.combo[,2:595]))

# --------------------------------------------------------------- SVA
model <- model.matrix(~as.factor(X), data = phenoData)
null <- model.matrix(~1, data = phenoData)

svobj <- sva(expData,model,null)

# --------------------------------------------------------------- Limma
limobj <- removeBatchEffect(data.combo[,2:595])


limobj[,1] <- as.numeric(as.factor(limobj[,1]))
limobj[,2:919] <- as.numeric(limobj[,2:919])

# --------------------------------------------------------------- Scanorama

# --------------------------------------------------------------- SVA + Limma

# --------------------------------------------------------------- EdgeR

# --------------------------------------------------------------- ComBat

Cbtobj <- ComBat(data.combo[,2:(595)], batchno)

# --------------------------------------------------------------- PCA

Lcombo.pca <- prcomp(t(data.combo[,2:(595)]))

PCAscores <- Lcombo.pca$x
PCAloadings <- Lcombo.pca$rotation
PCAcolors <- c("green","violet","blue", "pink", "yellow", "red", 'orange', "cyan", "grey0")

xmin <- -300000
xmax <- 1750000
ymin <- -2.5e+5
ymax <- 6e+5


plot(PCAscores[1:39,1:2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(1)],   
     bg=PCAcolors[c(1)],
     cex=1.0,         
     main="Scores"     
)
par(new=T)
plot(PCAscores[40:93,1:2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(2)],   
     bg=PCAcolors[c(2)],
     cex=1.0,         
     main="Scores"     
)
par(new=T)
plot(PCAscores[94:140,1:2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(3)],   
     bg=PCAcolors[c(3)],
     cex=1.0,         
     main="Scores"     
)
par(new=T)
plot(PCAscores[141:198,1:2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(4)],   
     bg=PCAcolors[c(4)],
     cex=1.0,         
     main="Scores"     
)
par(new=T)
plot(PCAscores[199:292,1:2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(5)],   
     bg=PCAcolors[c(5)],
     cex=1.0,         
     main="Scores"
     
)
par(new=T)
plot(PCAscores[293:494,1:2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(6)],   
     bg=PCAcolors[c(6)],
     cex=1.0,         
     main="Scores"
)
par(new=T)
plot(PCAscores[495:500,1:2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(7)],   
     bg=PCAcolors[c(7)],
     cex=1.0,        
     main="Scores"
     
)
par(new=T)
plot(PCAscores[501:531,1:2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(8)],   
     bg=PCAcolors[c(8)],
     cex=1.0,         
     main="Scores"
)
par(new=T)
plot(PCAscores[532:594,1:2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(9)],   
     bg=PCAcolors[c(9)],
     cex=1.0,         
     main="Scores",
)
par(new=T)
plot(PCAscores[532:594,1:2],xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     pch=21,
     col=PCAcolors[c(9)],   
     bg=PCAcolors[c(9)],
     cex=1.0,         
     main="Scores",
     legend("topright",
            c("L8","L3","L4","L14","L1","L6","L9","L7","L15"),
            fill=c("green","violet","blue", "pink", "yellow", "red", 'orange', "cyan", "grey0"),cex = 0.5)
)


plot(PCAloadings[,1:2], 
     pch=21,             
     bg="black",        
     cex=1,             
     main="Loadings"  
)
text(PCAloadings[,1:2],           
     labels=rownames(PCAloadings) 
)
