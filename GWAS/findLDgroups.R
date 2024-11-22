#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: findLDgroups.R
##
## Purpose of script: Find groups of variants with r2 linkage value superior to X
##
## Author: Victor Loegler
##
## Date Created: 2022-05-06
##
## ---------------------------
##
## Notes: Takes as argument 
##        - the list of SNPs, one SNP per line in file
##        - the matrix of r2, same order than the SNP list, without headers
##        - The path to the output file
##        - minimum value of R2 to group snps
##
## ---------------------------
#library(reshape2)
library("reshape2", lib.loc="/ccc/work/cont007/fg0006/loeglerv/Soft/Rlibraries")
## ---------------------------

# Get paths of input files
args = commandArgs(trailingOnly=TRUE)
snpsPath <- args[1]
ldPath <- args[2]
outPath <- args[3]
r2threshold <- as.numeric(args[4])

signifSNPs <- read.table(file = snpsPath)[,1]
ld <- read.table(file = ldPath, header = F, sep = "\t")
colnames(ld) <- signifSNPs
rownames(ld) <- signifSNPs
ld$SNP <- rownames(ld)

data <- melt(ld, id.vars = "SNP")
colnames(data) <- c("SNP1", "SNP2", "R2")
data$SNP1 <- as.character(data$SNP1)
data$SNP2 <- as.character(data$SNP2)
data <- data[data$SNP1 != data$SNP2,]

# Clustering
# Hierarchical clustering, average method
distMatrix <- ld
distMatrix[,"SNP"] <- NULL
maxR2 <- max(data$R2)
while (maxR2 >= r2threshold){
  SNP1 <- data$SNP1[data$R2 == max(data$R2)][1]
  SNP2 <- data$SNP2[data$R2 == max(data$R2)][1]
  newSNP <- paste(SNP1, SNP2, sep = ",")
  distMatrix[newSNP,] <- ( distMatrix[SNP1,] + distMatrix[SNP2,] ) / 2
  distMatrix[,newSNP] <- ( distMatrix[,SNP1] + distMatrix[,SNP2] ) / 2
  distMatrix <- distMatrix[-match(c(SNP1, SNP2), rownames(distMatrix)),]
  distMatrix[,c(SNP1, SNP2)] <- NULL

  distMatrix$SNP <- rownames(distMatrix)
  data <- melt(distMatrix, id.vars = "SNP")
  colnames(data) <- c("SNP1", "SNP2", "R2")
  data$SNP1 <- as.character(data$SNP1)
  data$SNP2 <- as.character(data$SNP2)
  data <- data[data$SNP1 != data$SNP2,]
  distMatrix[,"SNP"] <- NULL

  maxR2 <- max(data$R2)
}


# Write to output
write(colnames(distMatrix), outPath)

