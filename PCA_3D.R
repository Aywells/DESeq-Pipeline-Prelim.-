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

library(PEPPER)

# Data set specific parameters
geo.id <- "GSE7621" # GEO id of the data set
probe.conversion <- "ENTREZ_GENE_ID" # column name for gene id mapping
conversion.map <- NULL # probe to gene mapping, if NULL uses the mapping in the data set
conversion.mapping.function <- NULL # modify probe names using this function 
sample.mapping.column <- "characteristics_ch1" # column to use for sample mapping
geo.id.sub <- NULL # the platform to use if there are multiple platform annotations
reprocess <- "affy" # reprocessing type for raw data
output.dir <- "./"

# Get the expression and sample mapping info from the reprocessed data set
d <- fetch.expression.data(geo.id, sample.mapping.column = sample.mapping.column, do.log2 = NULL, probe.conversion = probe.conversion, conversion.map = conversion.map, conversion.mapping.function = conversion.mapping.function, output.dir = output.dir, geo.id.sub = geo.id.sub, reprocess = reprocess)

expr <- d$expr
sample.mapping <- d$sample.mapping

# Get z scores
out.file <- "z.dat"
cutoff <- 2.5

z = get.z.matrix(expr, sample.mapping, method="mean", out.file=out.file)
indices <- apply(abs(z), 2, function(x) { which(x >= cutoff)})
geneids <- lapply(indices, function(x) { rownames(z)[x] })

# Alteratively you can use get.peeps.from.z.matrix
peeps <- get.peeps.from.z.matrix(z, cutoff=2, convert.to.pvalues=F) 

