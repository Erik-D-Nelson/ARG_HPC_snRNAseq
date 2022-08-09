library(SingleCellExperiment)
library(scRNAseq)
library(EnsDb.Mmusculus.v79)
library(scater)
library(scran)
library(scry)
library(DropletUtils)
library(jaffelab)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(bluster)
library(pheatmap)

load("/users/enelson/sceSubset.rda")


### Clustering ============================================================
g25 <- buildSNNGraph(sce.subset, k=25, use.dimred = 'PCA')
clust25 <- igraph::cluster_walktrap(g25)$membership
colData(sce.subset)$k_25_label <- factor(clust20)
table(colData(sce.subset)$k_25_label)

###Annotation=============================================================
annotationTab<-data.frame("cluster"=c(1:28),"annotation"=NA)
annotationTab$annotation[c(26,15,9)] <- paste0("GABA.", c("3","4","5"))
annotationTab$annotation[c(21,8)] <- paste0("GABA.", c("1","2"))
annotationTab$annotation[c(6,4)]<-paste0("DG.", c("1","2"))
annotationTab$annotation[c(22,2,17,18,5)]<-c('CA4','CA3.1','CA3.2','CA2','CA1')
annotationTab$annotation[c(3,11)]<-c("PS.1","PS.2")
annotationTab$annotation[13]<-("Sub")
annotationTab$annotation[c(1,27,24,10,14,23,16)]<-paste0("L2/3.", c("1","2","3","4","5","6","7"))
annotationTab$annotation[c(20,25)]<-paste0("L5/Po.", c("1","2"))
annotationTab$annotation[c(28,12,7)]<-paste0("L6.", c("1","2","3"))
annotationTab$annotation[c(19)]<-("L6b")

annotationTab$cluster<-factor(annotationTab$cluster)
annotationTab$annotation<-factor(annotationTab$annotation)
sce.subset$annotation <- annotationTab$annotation[match(sce.subset$k_25_label,
                                                           annotationTab$cluster)]

###Marker gene detection=====================================================
marks<-findMarkers(sce.subset,pval.type='all',
                  direction='up',group=sce.subset$annotation)

##make cellType factor for plots
sce.subset$cellType<-as.character(sce.subset$annotation)
sce.subset$cellType<-ifelse(sce.subset$k_25_label %in% c(1,27,24,10,14,23,16),'L2/3',
         ifelse(sce.subset$k_25_label %in% c(20,25),'L5/Po',
                ifelse(sce.subset$k_25_label %in% c(28,12,7,19),'L6/6b',
                                sce.subset$cellType)))
                                                                        
sce.subset$cellType<-factor(sce.subset$cellType,
                            levels=c('DG.1','DG.2','CA4','CA3.1','CA3.2','CA2',
                                     'CA1','PS.1','PS.2','Sub',
                                     'L2/3','L5/Po','L6/6b','GABA.1',
                                     'GABA.2','GABA.3','GABA.4','GABA.5')
)

plotUMAP(sce.subset,colour_by='cellType',text_by='label')
###Figure 1 plots===============================================================

##UMAP by annotation for figure 1
pdf('20210913_UMAP_byAnnotation.pdf',h=7,w=7)
plotUMAP(sce.subset,colour_by='k_20_label',
         text_by='annotation',text_size=3,add_legend=F) 
dev.off()

##dotplot for figure 1
features=c('Slc17a7','Gad2','Prox1',
           'Calb2','Nectin3','St18',
           'Ntsr2','Man1a','Dcn','Ndst4','Fn1','Satb2',
           'Pamr1','Foxp2','Ndnf','Htr3a',
           'Lhx6','Lamp5','Pvalb','Sst')

pdf('dotPlot_clusterMarkers.pdf',h=6.5,w=12)
plotDots(sce,group='cellType',features=features,color=c('white','red')) + 
  scale_y_discrete(limits=rev(features)) +
  theme(axis.text.x = element_text(angle = 45,vjust=0.75),text = element_text(size = 20))
dev.off()

CA3.1_cor<-list()
for(i in 1:length(enrich_CA3.1$geneID)){
CA3.1_cor[[i]]<-as.vector(strsplit(enrich_CA3.1$geneID, "/")[[i]])}


##CA field and DG-specific violin plots by cellType for figure S1
markers.custom<-list(
  'CA fields + dentate gyrus' = 'Zbtb20',
  'CA1' = 'Fibcd1',
  'DG' = 'Prox1',
  'CA3 + CA4' = 'Galnt3'
)
pdf('20210923_violinPlots_bycellType.pdf',h=6,w=4)
for(i in 1:length(markers.custom)){
  print(
    plotExpression(sce.subset, exprs_values = "logcounts", features=c('Zbtb20'),
                   x="k_25_label", colour_by="k_20_label", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
                   width = 0.3) + 
      ggtitle(label=(names(markers.custom)[i])) +
      theme(axis.text.x = element_text(angle = 90))
  )}
dev.off()

save(sce.subset, hdg, pc.choice.hpc,marks,file="/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&20_clusteredPCs.rda")
##Proceed to figure 2 plots

sce$annotation<-factor(sce$annotation,
    levels=c("DG.1","DG.2","CA4","CA3.1","CA3.2","CA2",
              "CA1","PS.1","PS.2","Sub","L2/3.1","L2/3.2",
              "L2/3.3","L2/3.4","L2/3.5","L2/3.6","L2/3.7",
              "L5/Po.1","L5/Po.2","L6.1","L6.2","L6.3","L6b",
              "GABA.1","GABA.2","GABA.3","GABA.4","GABA.5"))

sce$cellType<-factor(sce$cellType,
                       levels=c("DG.1","DG.2","CA4","CA3.1","CA3.2","CA2",
                                "CA1","PS.1","PS.2","Sub","L2/3",
                                "L5/Po","L6/6b",
                                "GABA.1","GABA.2","GABA.3","GABA.4","GABA.5"))
