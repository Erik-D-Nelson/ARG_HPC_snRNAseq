###Code for running scCoGAPS on 500 randomly subsetd nuclei
###We're going to start with 50 patterns

library(SingleCellExperiment)
library(scRNAseq)
library(EnsDb.Mmusculus.v79)
library(scater)
library(scran)
library(CoGAPS)

##load sce.subset
load("/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50_clusteredPCs.rda")

##subset sce.subset
sce.subset<-sce.subset[,subset(ncol(sce.subset),500)]

##set parameters
params <- new("CogapsParams",sparseOptimization=T)
params <- setDistributedParams(params, nSets=5)
params <- setParam(params, "nIterations", 1000)
params <- setParam(params, "nPatterns", 50)

##run scCoGAPS
results<-scCoGAPS(data, params, messages=T, transposeData=F)

##save everything
save(sce.subset,params,results,file='test_scCoGAPS.rda')
sessionInfo()

cellType<-list()

for (i in 1:length(levels(sce.subset$cellType))){
  cellType[[i]]<-sce.subset$cellType==levels(sce.subset$cellType)[i]
}
names(cellType)<-levels(sce.subset$cellType)
dimnames<-list(colnames(sce.subset),names(cellType))
cellType<-matrix(unlist(cellType),ncol=15,nrow=15769,dimnames=dimnames)
map_cellType<-cor(cellType,result)
map_cellType<-map_cellType[,c(order)]
pheatmap(map_cellType,cluster_cols=F,cluster_rows=F)

cellType<-list()

for (i in 1:length(levels(sce.subset$cellType))){
  cellType[[i]]<-sce.subset$cellType==levels(sce.subset$cellType)[i]
}
names(cellType)<-levels(sce.subset$cellType)
dimnames<-list(colnames(sce.subset),names(cellType))
cellType<-matrix(unlist(cellType),ncol=15,nrow=15769,dimnames=dimnames)
map_cellType<-cor(cellType,result)
map_cellType<-map_cellType[,c(order)]
pheatmap(map_cellType,cluster_cols=F,cluster_rows=F)

condition<-list()
for (i in 1:length(levels(sce.subset$condition))){
  condition[[i]]<-sce.subset$condition==levels(sce.subset$condition)[i]
}
names(condition)<-levels(sce.subset$condition)
dimnames<-list(colnames(sce.subset),names(condition))
condition<-matrix(unlist(condition),ncol=2,nrow=15769,dimnames=dimnames)

map_condition<-cor(condition,result)
map_condition<-map_condition[,order]
pheatmap(map_condition,cluster_cols=F,cluster_rows=F)


##clusterprofiler
enrich_go <- enrichGO(gene = x$ID[x$FDR  < 0.05],
                      OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "BP",
                      pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)

## Visualize enrichment results
pdf('barplot_GOanalysis_30cats.pdf')
barplot(enrich_go, font.size = 7,showCategory=30)
dev.off()


head(sort(features[,12], decreasing=TRUE), n=10)


fixTo <- rownames(marks[[1]])

for(s in names(marks)){
  marks[[s]]$rank <- 1:nrow(marks[[s]])
  marks[[s]] <- marks[[s]][fixTo, ]
}

ts <- sapply(marks, function(x){x$rank})
rownames(ts) <- fixTo
ts <- ts[rownames(ts) %in% rowData(sce.hsap.sub)$Symbol.uniq, ]
rownames(ts) <- rowData(sce.hsap.sub)$JAX.geneID[match(rownames(ts), rowData(sce.hsap.sub)$Symbol.uniq)]

features <- features[match(rownames(marks[[1]]), rownames(features)), ]

sce.subset$label_seizure<-as.character(sce.subset$annotation)
sce.subset$label_seizure<-factor(ifelse(sce.subset$condition=='seizure',
                                 paste0(sce.subset$annotation,'.seizure'),
                                 paste0(sce.subset$annotation,'.control')))
paste0(sce.fig3$cellType,'.seizure')

rowData(sce.subset)$EG = mapIds(org.Mm.eg.db,
                     keys=rowData(sce.subset)$ID,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

#get list of lists of nuclei by cell type
cellType<-list()
for (i in 1:length(levels(sce.subset$cellType))) {
  cellType[[i]]<-colnames(sce.subset[,sce.subset$cellType==levels(sce.subset$cellType)[i]])}

condition<-list()
for (i in 1:length(levels(sce.subset$condition))) {
  condition[[i]]<-colnames(sce.subset[,sce.subset$condition==levels(sce.subset$condition)[i]])}