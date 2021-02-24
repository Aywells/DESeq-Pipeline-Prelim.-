
library(NomiX)
ls("package:NomiX")

require(data.table)
require(magrittr)
require(dplyr)
require(pROC)
require(ggplot2)
require(OptimalCutpoints)

setwd('C:/Users/aywel/Desktop/Research/Barabasi Lab/Nomix/Code/Data (New)')
data.total <- read.csv("Raw_cl.csv")    # RAW DATA
data.pheno <- data.total[,3:15]    # Phenotypic logicals for diseases
data.counts <- data.total[,19:4653]    # Raw counts for every sample
data.batch <- data.total[,1:2]    # Just sample/batch IDs

data.norm.counts <- as.data.frame(scale(data.counts)) #temporary calulation

#DEG Function
deg <- DEG(Data = data.norm.counts,
          Disease_vector = as.vector(data.pheno$Asthma),
          k = 8,
          Disease_name = 'Asthma',
          method = "MaxDOR")

#DEG-Aux Function
deg2 <- DEG_aux(
  i=8,
  Disease_vector = as.vector(data.pheno[,5]),
  Disease_name = 'Asthma',
  Data = data.norm.counts,
  method = "MaxDOR")

#Score-Patient Function
score_P <- Score_Patient(
  Disease = "Asthma",
  data.counts,
  dir = deg,
  normalize = T
  )

#Score-Individual Function
score_I <- Score_Individual(
  Disease = "Asthma",
  RNAseq = data.norm.counts,
  dir = "C:/Users/aywel/Desktop/Research/Barabasi Lab/Nomix/Code/Data (New)/DEGs"
)

