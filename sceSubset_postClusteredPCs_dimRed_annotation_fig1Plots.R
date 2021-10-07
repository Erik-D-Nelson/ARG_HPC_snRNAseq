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

load("/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50_clusteredPCs.rda")
##getClusteredPCs() suggested an optimal # of PCs to be 45 for 10000 HDG

### Dimensionality Reduction ============================================================
set.seed(1000) 
sce.subset <- runPCA(sce.subset,exprs_values="poisson_pearson_residuals", 
                     subset_row=hdg, BSPARAM=BiocSingular::RandomParam(),ncomponents=45) 
reducedDimNames(sce.subset)

set.seed(100000)
sce.subset <- runUMAP(sce.subset, dimred="PCA")

#Verify length
ncol(reducedDim(sce.subset, "PCA"))

g50 <- buildSNNGraph(sce.subset, k=50, use.dimred = 'PCA')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(sce.subset)$k_50_label <- factor(clust50)
table(colData(sce.subset)$k_50_label)

annotationTab<-data.frame("cluster"=c(1:23),"annotation"=NA)
annotationTab$annotation[c(11,21,8,13)] <- paste0("GABA.", c("1","2","3","4"))
annotationTab$annotation[c(22,17,23,6,19,18)]<-("RHP.L2/3.")# , c("1","2","3","4","5","6"))
annotationTab$annotation[c(9,16,15,1)]<-paste0("RHP.L5/6/6b.", c("1","2","3","4"))
annotationTab$annotation[c(2,14)]<-paste0("CA3.", c("1","2"))
annotationTab$annotation[c(5,3)]<-paste0("DG.", c("1","2"))
annotationTab$annotation[c(7,4,10,12,20)]<-c("CA1.1","CA1.2","Sub","CA2","CA4")
annotationTab$cluster<-factor(annotationTab$cluster)
annotationTab$annotation<-factor(annotationTab$annotation)
sce.subset$annotation.23 <- annotationTab$annotation[match(sce.subset$k_50_label,
                                                           annotationTab$cluster)]
##findMarkers by annotated cluster
marks<-findMarkers(sce.subset,pval.type='all',
                   direction='up',group=sce.subset$annotation)

##make cellType factor for plots
sce.subset$cellType<-as.character(sce.subset$annotation)
sce.subset$cellType<-ifelse(sce.subset$k_50_label %in% c(6,17,18,19,22,23),'RHP.L2/3',
         ifelse(sce.subset$k_50_label %in% c(1,9,15,16),'RHP.L5/6/6b',
                ifelse(sce.subset$annotation=='CA1.1','CA1.D',
                       ifelse(sce.subset$annotation=='CA1.2','CA1.V',
                              ifelse(sce.subset$annotation=='CA3.1','CA3.D',
                                     ifelse(sce.subset$annotation=='CA3.2','CA3.V',
                                            ifelse(sce.subset$annotation=='DG.1','DG.seizure',
                                                   ifelse(sce.subset$annotation=='DG.2','DG.control',
                                                          ifelse(sce.subset$annotation=='GABA.1','GABA.Sst',
                                                                 ifelse(sce.subset$annotation=='GABA.2','GABA.PV',
                                                                        ifelse(sce.subset$annotation=='GABA.3','GABA.NGF',
                                                                               ifelse(sce.subset$annotation=='GABA.4','GABA.CGE',
                                                                                      sce.subset$cellType)
                                                                        ))))))))))))
sce.subset$cellType<-factor(sce.subset$cellType,
                            levels=c('CA1.D','CA1.V','CA2',
                                     'CA3.D','CA3.V','CA4',
                                     'Sub','DG.control','DG.seizure',
                                     'RHP.L2/3','RHP.L5/6/6b','GABA.Sst',
                                     'GABA.PV','GABA.NGF','GABA.CGE')
)

###Figure 1 plots===============================================================

##UMAP by annotation for figure 1
pdf('20210913_UMAP_byAnnotation.pdf',h=7,w=7)
plotUMAP(sce.subset,colour_by='k_50_label',
         text_by='annotation',text_size=3,add_legend=F) 
dev.off()

##dotplot for figure 1
features=c('Slc17a7','Gad2','Fibcd1',
           'Strip2','Nectin3','Col6a1',
           'Glis1','Fn1','Prox1',
           'Satb2','Foxp2','Sst',
           'Pvalb','Lamp5','Htr3a')

pdf('20210924_dotPlot_clusterMarkers.pdf',h=6,w=6)
plotDots(sce.subset,group='cellType',features=features) + 
  scale_y_discrete(limits=rev(features)) +
  theme(axis.text.x = element_text(angle = 45))
dev.off()


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
    plotExpression(sce.subset, exprs_values = "logcounts", features=markers.custom[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
                   width = 0.3) + 
      ggtitle(label=(names(markers.custom)[i])) +
      theme(axis.text.x = element_text(angle = 90))
  )}
dev.off()

save(sce.subset, hdg, pc.choice.hpc,marks,file="/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50_clusteredPCs.rda")
##Proceed to figure 2 plots
