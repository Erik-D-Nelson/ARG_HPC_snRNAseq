###Code for running scCoGAPS on 500 randomly sampled nuclei
###We're going to start with 50 patterns

library(SingleCellExperiment)
library(scRNAseq)
library(EnsDb.Mmusculus.v79)
library(scater)
library(scran)
library(CoGAPS)

##load sce.subset
load("/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50_clusteredPCs.rda")

##sample sce.subset
sce.sample<-sce.subset[,sample(ncol(sce.subset),500)]

##set parameters
params <- new("CogapsParams")
params <- setDistributedParams(params, nSets=5)
params <- setParam(params, "nIterations", 1000)
params <- setParam(params, "nPatterns", 50)

##run scCoGAPS
results<-scCoGAPS(sce.sample, params, messages=T, transposeData=F)

##save everything
save(sce.sample,params,results,file='test_scCoGAPS.rda')
sessionInfo()