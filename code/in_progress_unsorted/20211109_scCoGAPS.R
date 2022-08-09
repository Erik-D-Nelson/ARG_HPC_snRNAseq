###Code for running scCoGAPS on all nuclei
###We're going to start with 80 patterns

library(SingleCellExperiment)
library(scRNAseq)
library(EnsDb.Mmusculus.v79)
library(scater)
library(scran)
library(CoGAPS)

##load sce.subset
load("/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50_clusteredPCs.rda")

##Pull out logcounts from sce.subset
data<-assay(sce.subset,'logcounts')
data<-as.matrix(data)

##set parameters
params <- new("CogapsParams",sparseOptimization=T)
params <- setDistributedParams(params, nSets=100)
params <- setParam(params, "nIterations", 1000)
params <- setParam(params, "nPatterns", 80)

##run scCoGAPS
results<-CoGAPS(data, params, messages=T, transposeData=F, distributed="single-cell")

##save everything
save(params,results,file='20211109_scCoGAPS.rda')
sessionInfo()