
library(SingleCellExperiment)
library(CoGAPS)

##load sce.subset
load("/users/enelson/newSceSubset_geneFiltered_under20_asMatrix.rda")

##set parameters
params <- new("CogapsParams",sparseOptimization=T)
params <- setDistributedParams(params, nSets=20)
params <- setParam(params, "nIterations", 12000)
params <- setParam(params, "nPatterns", 70)
params <- setParam(params, "singleCell", T)

##run scCoGAPS
results<-CoGAPS(data=data, params, messages=T, transposeData=F, distributed='single-cell', outputFrequency=100, 
checkpointInterval=500, checkpointOutFile="/users/enelson/scCoGAPS.out")

##save everything
save(params,results,file='processed_data/2022_scCoGAPS.rda')
sessionInfo()