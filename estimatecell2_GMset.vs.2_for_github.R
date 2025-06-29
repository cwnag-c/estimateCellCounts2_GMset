library(minfi)
library(wateRmelon)
library(data.table)
library(magrittr)
library(FlowSorted.Blood.EPIC)
#library(FlowSorted.Blood.450k)
#-------------------------set pathway-------------#
source("./estimateCellCounts2_GMset")
#------------------GSE201287------
targets<-read.csv("./GSE201287_MDD_Blood_pheno.csv")
rgSet <- read.metharray.exp(targets = targets,extended = T)
targets$Sample_Name<-paste(targets$Sample_Name,targets$Slide,targets$Array,sep = "_")#与rgSet的samplenames对齐
targets<- targets[match(sampleNames(rgSet), targets$Sample_Name),]

referenceset<-load("./FlowSorted.Blood.EPIC.rda")
cellp <-estimateCellCounts2(rgSet,
                            compositeCellType = "Blood",
                            processMethod = "preprocessQuantile",
                            probeSelect = "IDOL",
                            referenceset=referenceset,
                            cellTypes = base::c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"))
cellp <-cellp$prop
M<-preprocessRaw(rgSet)
GM<-mapToGenome(M)
all(match(sampleNames(GM),targets$Sample_Name))
GM$Sex<-targets$Gender
cellp2 <-estimateCellCounts2_GMset(GM,referenceset=referenceset)
cellp2 <-cellp2$prop

df<-cellp-cellp2
