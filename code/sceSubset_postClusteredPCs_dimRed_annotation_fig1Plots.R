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

load("/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&20_clusteredPCs.rda")
##getClusteredPCs() suggested an optimal # of PCs to be 45 for 10000 HDG

### Dimensionality Reduction ============================================================
set.seed(1000) 
sce.subset <- runPCA(sce.subset,exprs_values="poisson_pearson_residuals", 
                     subset_row=hdg) 
reducedDimNames(sce.subset)

set.seed(100000)
sce.subset <- runUMAP(sce.subset, dimred="PCA")

#Verify length
ncol(reducedDim(sce.subset, "PCA"))

g20 <- buildSNNGraph(sce.subset, k=20, use.dimred = 'PCA')
clust20 <- igraph::cluster_walktrap(g20)$membership
colData(sce.subset)$k_20_label <- factor(clust20)
table(colData(sce.subset)$k_20_label)

g10 <- buildSNNGraph(sce.subset, k=10, use.dimred = 'PCA')
clust10 <- igraph::cluster_walktrap(g10)$membership
colData(sce.subset)$k_10_label <- factor(clust10)
table(colData(sce.subset)$k_10_label)

g5 <- buildSNNGraph(sce.subset, k=5, use.dimred = 'PCA')
clust5 <- igraph::cluster_walktrap(g5)$membership
colData(sce.subset)$k_5_label <- factor(clust5)
table(colData(sce.subset)$k_5_label)

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
#findMarkers by annotated cluster
#marks<-findMarkers(sce.subset,pval.type='all',
#                  direction='up',group=sce.subset$cellType)

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

features=c('Slc17a7','Prox1',
           'Calb2','Nectin3','Man1a','Dcn','Fn1','Satb2',
            'Foxp2','Gad2')

pdf('dotPlot_clusterMarkers.pdf',h=6.5,w=12)
plotDots(sce,group='cellType',features=features,color=c('white','red')) + 
  scale_y_discrete(limits=rev(features)) +
  theme(axis.text.x = element_text(angle = 45,vjust=0.75),text = element_text(size = 20))
dev.off()

CA3.1_cor<-list()
for(i in 1:length(enrich_CA3.1$geneID)){
CA3.1_cor[[i]]<-as.vector(strsplit(enrich_CA3.1$geneID, "/")[[i]])}


##CA field and DG-specific violin plots by cellType for figure 1
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

#get list of


annotationTab<-data.frame("cluster"=c(1:28),"annotation"=NA)
annotationTab$annotation[c(8,12,16,23)] <- paste0("GABA.", c("1","2","3","4","5"))
annotationTab$annotation[c(22,17,23,6,19,18)]<-("RHP.L2/3.")# , c("1","2","3","4","5","6"))
annotationTab$annotation[c(9,16,15,1)]<-paste0("RHP.L5/6/6b.", c("1","2","3","4"))
annotationTab$annotation[c(2,14)]<-paste0("CA3.", c("1","2"))
annotationTab$annotation[c(5,3)]<-paste0("DG.", c("1","2"))
annotationTab$annotation[c(7,4,10,12,20)]<-c("CA1.1","CA1.2","Sub","CA2","CA4")
annotationTab$cluster<-factor(annotationTab$cluster)
annotationTab$annotation<-factor(annotationTab$annotation)
sce.subset$annotation <- annotationTab$annotation[match(sce.subset$k_20_label,
                                                        annotationTab$cluster)]


#fde725, #d8e219, #b0dd2f
#89d548
#65cb5e
#46c06f
#2eb37c
#21a585
#1f978b
#23898e
#297b8e
#2e6d8e
#355e8d
#3d4e8a
#433d84
#472a7a
#481769
#440154

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
 
"GABA.1"  "GABA.2"  "GABA.3"  "GABA.4"  "GABA.5" "L5/Po.1" "L5/Po.2"
"L6.1"    "L6.2"    "L6.3"    "L6b"    
