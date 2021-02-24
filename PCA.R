library(xopen) # General Library
library(readxl)
library(readr)
library(plyr)
library(openxlsx)
library(limma)
library(ggplot2)
library(ggfortify)
library(dplyr)

library(GEOquery) # Gene Library
library(Biobase)
library(BiocGenerics)
library(BiocManager)
library(biomaRt)

# Read specific tables and text for sufficent subject mtx.
L1 <- read_excel('GSE64813_raw.xlsx')        #*
# L2 <- read.csv('GSE97356_raw.csv')
L3 <- read_excel('GSE107991_raw.xlsx')       #*
L4 <- read_excel('GSE107992_raw.xlsx')       #*
# L5 <- read_excel('GSE107994_raw.xlsx')
L6 <- read.table('GSE112057_raw.txt')        #*
L7 <- read.csv('GSE114037_raw.csv')          #*
L8 <- read.csv('GSE114407_raw.csv')          #*
L9 <- read_excel('GSE33701_raw.xlsx')        #*
# L10 <- read.table('GSE85531_raw.txt')
# L11 <- read_excel('GSE110487_raw.xlsx')
# L12 <- read.table('GSE112087_raw.txt')
# L13 <- read.table('GSE112594_raw.txt')
L14 <- read.csv('GSE120913_raw.csv')         #*
L15 <- read_excel('GSE124180_raw.xlsx')      #*
# L16 <- read.csv('GSE124284_raw.csv')
# L17 <- read.csv('GSE124400_raw.csv')

# Sets row 1 to header for specific files
header.true <- function(df) {
        names(df) <- as.character(unlist(df[1,]))
        df[-1,]
}
L6 <- header.true(L6)
# L10 <- header.true(L10)
# L12 <- header.true(L12)
# L13 <- header.true(L13)

# Transulator for ensemble ID's to gene symbols w/ biomaRt
# ensemblsIDS <- L11[,1]
# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# new <- getBM(attributes='hgnc_symbol', 
#       filters = 'ensembl_gene_id', 
#       values = ensemblsIDS, 
#       mart = ensembl)

# convert usable files to numeric
i1 <- c(2:95)
L1[ , i1] <- apply(L1[ , i1], 2, function(x) as.numeric(as.character(x)))

#i2 <- c(2:325)
#L2[ , i2] <- apply(L2[ , i2], 2, function(x) as.numeric(as.character(x)))

i3 <- c(2:55)
L3[ , i3] <- apply(L3[ , i3], 2, function(x) as.numeric(as.character(x)))

i4 <- c(2:48)
L4[ , i4] <- apply(L4[ , i4], 2, function(x) as.numeric(as.character(x)))

i6 <- c(2:203)
L6[ , i6] <- apply(L6[ , i6], 2, function(x) as.numeric(as.character(x)))

i7 <- c(2:32)
L7[ , i7] <- apply(L7[ , i7], 2, function(x) as.numeric(as.character(x)))

i8 <- c(2:40)
L8[ , i8] <- apply(L8[ , i8], 2, function(x) as.numeric(as.character(x)))

i9 <- c(2:7)
L9[ , i9] <- apply(L9[ , i9], 2, function(x) as.numeric(as.character(x)))

i14 <- c(2:59)
L14[ , i14] <- apply(L14[ , i14], 2, function(x) as.numeric(as.character(x)))

i15 <- c(2:64)
L15[ , i15] <- apply(L15[ , i15], 2, function(x) as.numeric(as.character(x)))


#combo <- inner_join(unique(L2),unique(L8), by = 'X')
combo2 <- inner_join(unique(L8),unique(L3), by = 'X')
combo3 <- inner_join(unique(combo2),unique(L4), by = 'X')
combo4 <- inner_join(unique(combo3),unique(L14), by = 'X')
combo5 <- inner_join(unique(combo4),unique(L1), by = 'X')
combo6 <- inner_join(unique(combo5),unique(L6), by = 'X')
combo7 <- inner_join(unique(combo6),unique(L9), by = 'X')
combo8 <- inner_join(unique(combo7),unique(L7), by = 'X')
combo9 <- inner_join(unique(combo8),unique(L15), by = 'X')


x <- combo9[1,]
y <- combo9[,1]
data <- expand.grid(X=x, Y=y)

# Heatmap 
ggplot(data, aes(X, Y)) + 
        geom_tile()


Lcombo.pca <- prcomp(t(combo9[, 2:919]))

PCAscores <- Lcombo.pca$x
PCAloadings <- Lcombo.pca$rotation
PCAcolors <- c("green","violet","blue", "pink", "yellow", "red", 'orange', "cyan", "grey0","bisque")

plot(PCAscores[0:324,1:2],xlim=c(-300000,1750000), ylim=c(-2.5e+5,6e+5),
     pch=21,
     col=PCAcolors[c(1)],   
     bg=PCAcolors[c(1)],
     cex=1.0,         
     main="Scores"     
)
par(new=T)
plot(PCAscores[325:363,1:2],xlim=c(-300000,1750000), ylim=c(-2.5e+5,6e+5),
     pch=21,
     col=PCAcolors[c(2)],   
     bg=PCAcolors[c(2)],
     cex=1.0,         
     main="Scores"     
)
par(new=T)
plot(PCAscores[364:417,1:2],xlim=c(-300000,1750000), ylim=c(-2.5e+5,6e+5),
     pch=21,
     col=PCAcolors[c(3)],   
     bg=PCAcolors[c(3)],
     cex=1.0,         
     main="Scores"     
)
par(new=T)
plot(PCAscores[418:464,1:2],xlim=c(-300000,1750000), ylim=c(-2.5e+5,6e+5),
     pch=21,
     col=PCAcolors[c(4)],   
     bg=PCAcolors[c(4)],
     cex=1.0,         
     main="Scores"     
)
par(new=T)
plot(PCAscores[465:522,1:2],xlim=c(-300000,1750000), ylim=c(-2.5e+5,6e+5),
     pch=21,
     col=PCAcolors[c(5)],   
     bg=PCAcolors[c(5)],
     cex=1.0,         
     main="Scores"
     
)
par(new=T)
plot(PCAscores[523:616,1:2],xlim=c(-300000,1750000), ylim=c(-2.5e+5,6e+5),
     pch=21,
     col=PCAcolors[c(6)],   
     bg=PCAcolors[c(6)],
     cex=1.0,         
     main="Scores"
)
par(new=T)
plot(PCAscores[617:818,1:2],xlim=c(-300000,1750000), ylim=c(-2.5e+5,6e+5),
     pch=21,
     col=PCAcolors[c(7)],   
     bg=PCAcolors[c(7)],
     cex=1.0,        
     main="Scores"
     
)
par(new=T)
plot(PCAscores[819:824,1:2],xlim=c(-300000,1750000), ylim=c(-2.5e+5,6e+5),
     pch=21,
     col=PCAcolors[c(8)],   
     bg=PCAcolors[c(8)],
     cex=1.0,         
     main="Scores"
)
par(new=T)
plot(PCAscores[825:855,1:2],xlim=c(-300000,1750000), ylim=c(-2.5e+5,6e+5),
     pch=21,
     col=PCAcolors[c(9)],   
     bg=PCAcolors[c(9)],
     cex=1.0,         
     main="Scores"
)
par(new=T)
plot(PCAscores[856:918,1:2],xlim=c(-300000,1750000), ylim=c(-2.5e+5,6e+5),
     pch=21,
     col=PCAcolors[c(10)],   
     bg=PCAcolors[c(10)],
     cex=1.0,         
     main="Scores",
     legend("topright",
            c("L2","L3","L4","L8","L14","L1","L6","L9","L7","L15"),
            fill=c("green","violet","blue", "pink", "yellow", "red", 'orange', "cyan", "grey0","bisque"),cex = 0.5)
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


# L1.pca <- prcomp(L1[2:27975, 2:189])
# autoplot((prcomp(combo4[, 2:652]))
