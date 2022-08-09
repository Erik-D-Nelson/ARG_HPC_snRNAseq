################################################################################
### LIBD pilot 10x snRNA-seq: dg samples
###   **Region-specific analyses**
###   - R-batch job for detxn of optimal PC space with 'sce.DG' object
################################################################################

library(SingleCellExperiment)
#library(scRNAseq)
#library(EnsDb.Mmusculus.v79)
library(scater)
library(scran)
#library(scry)
#library(DropletUtils)
#library(jaffelab)
#library(gridExtra)
#library(GenomicRanges)
library(ggplot2)
library(tidyverse)
library(bluster)
library(pheatmap)
#library(GSEABase)
load("/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50_clusteredPCs.rda", verbose=T)

sce.DG<-sce.subset[,sce.subset$annotation %in% c('DG.1','DG.2')]
rm(sce.subset)
# ===
### Normalization ============================================================
set.seed(1000)
clusters <- quickCluster(sce.DG)
sce.DG <- computeSumFactors(sce.DG, cluster=clusters)
sce.DG <- logNormCounts(sce.DG)

### NEW Feature Selection ============================================================
#This is based on fitting Poisson model and then calculating Pearson residuals
#Let's plot highly deviant genes first; this might give us an idea about how 
#many genes to retain
sce.DG<-devianceFeatureSelection(sce.DG, assay="counts", 
                                    fam="poisson",sorted=TRUE)

plot(rowData(sce.DG)$poisson_deviance, type="l", 
     xlab="ranked genes",
     ylab="poisson deviance", 
     main="Feature Selection with Deviance",
     ylim=c(0,100000))
abline(v=2000,col='red')
abline(v=3000,col='blue')
abline(v=1000)


hdg<-rownames(counts(sce.DG))[1:2000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce.DG<-nullResiduals(sce.DG, assay = "counts",
                         fam = "poisson", type = "pearson")


### Dimensionality Reduction ============================================================
set.seed(1000) 
sce.DG <- runPCA(sce.DG,exprs_values="poisson_pearson_residuals", 
                    subset_row=hdg,ncomponents=50) 

set.seed(100000)
sce.DG <- runUMAP(sce.DG, dimred="PCA",subset_row=hdg)

reducedDimNames(sce.DG)

#Verify length
ncol(reducedDim(sce.DG, "PCA"))


# getClusteredPCs() to identify working PC space
pc.choice.dg <- getClusteredPCs(reducedDim(sce.DG, "PCA"))

# How many PCs should use in this space?
metadata(pc.choice.dg)$chosen
#47 was the number here

## Plot n Clusters vs. d PCs
pdf("20210915_sceDG_getClusteredPCs.pdf")
plot(pc.choice.dg$n.pcs, pc.choice.dg$n.clusters,
     main="PC Choice via getClusteredPCs()",
abline(v=metadata(pc.choice.dg)$chosen, col="red", lty="dashed", lwd=0.8)) 
dev.off()


set.seed(1000) 
sce.DG <- runPCA(sce.DG,exprs_values="poisson_pearson_residuals", 
                    subset_row=hdg, BSPARAM=BiocSingular::RandomParam(),ncomponents=50) 
reducedDimNames(sce.DG)

set.seed(100000)
sce.DG <- runUMAP(sce.DG, dimred="PCA",subset_row=hdg)

#Verify length
ncol(reducedDim(sce.DG, "PCA"))

g50 <- buildSNNGraph(sce.DG, k=50, use.dimred = 'PCA')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(sce.DG)$k_50_label <- factor(clust50)
table(colData(sce.DG)$k_50_label)

g25 <- buildSNNGraph(sce.DG, k=25, use.dimred = 'PCA')
clust25 <- igraph::cluster_walktrap(g25)$membership
colData(sce.DG)$k_25_label <- factor(clust25)
table(colData(sce.DG)$k_25_label)

#g10 <- buildSNNGraph(sce.DG, k=10, use.dimred = 'PCA')
#clust10 <- igraph::cluster_walktrap(g10)$membership
#colData(sce.DG)$k_10_label <- factor(clust10)
#table(colData(sce.DG)$k_10_label)

g20 <- buildSNNGraph(sce.DG, k=20, use.dimred = 'PCA')
clust20 <- igraph::cluster_walktrap(g20)$membership
colData(sce.DG)$k_20_label <- factor(clust20)
table(colData(sce.DG)$k_20_label)

g15 <- buildSNNGraph(sce.DG, k=15, use.dimred = 'PCA')
clust15 <- igraph::cluster_walktrap(g15)$membership
colData(sce.DG)$k_15_label <- factor(clust15)
table(colData(sce.DG)$k_15_label)
pdf('20210922_violin_DG_ventral.pdf',h=5,w=5)
plotExpression(sce.DG,features=c('Tenm4'),
               x="label", colour_by="label", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3)+
  ggtitle("Ventral granule cells")
dev.off()

pdf('20210922_violin_DG_dorsal1.pdf',h=5,w=5)
plotExpression(sce.DG,features=c('Ntng1'),
               x="label", colour_by="label", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3)+
  ggtitle("Dorsal granule cells 1")
dev.off()

pdf('20210922_violin_DG_dorsal2.pdf',h=5,w=5)
plotExpression(sce.DG,features=c('Inf2'),
               x="label", colour_by="label", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3)+
  ggtitle("Dorsal granule cells 2")
dev.off()

'Trhr','Nr2f2','Atp2b4','Gsg1l','Inf2'
# Save
save(sce.DG, hdg,
     file="/users/enelson/20210922_sceDG_postClusteredPCs.rda")

rm(list=ls())
sessionInfo()


dbl.dens <- computeDoubletDensity(sce.subset, subset.row=hdg, 
                                  d=ncol(reducedDim(sce.subset,'PCA')))
summary(dbl.dens)

sce.sub$cellType<-as.character(sce.sub$cellType)
sce.sub$cellType<-factor(
  ifelse(sce.sub$cellType=='CA1.V' & sce.sub$condition=='seizure','CA1.V.seizure',
         ifelse(sce.sub$cellType=='CA3.D' & sce.sub$condition=='seizure','CA3.D.seizure',
                ifelse(sce.sub$cellType=='DG.control','DG',
                       sce.sub$cellType))))








markers.custom=list(
  "Upregulated in DG"='Arc',
  "Upregulated in CA"='Mfap3l'
)
pdf('20210924_violinPlots_bycellType.pdf',h=5,w=5)
for(i in 1:length(markers.custom)){
  print(
    plotExpression(sce.sub, exprs_values = "logcounts", features=markers.custom[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
                   width = 0.3) + 
      ggtitle(label=(names(markers.custom)[i])) +
      theme(axis.text.x = element_text(angle = 90))
  )}
dev.off()
sce.subset$cellType<-factor(sce.subset$cellType,levels=c("CA1.D",       "CA1.V",       "CA2",         "CA3.D",       "CA3.V",
                                                         "CA4",         "Sub", "DG.control",  "DG.seizure",  "RHP.L2/3",    "RHP.L5/6/6b", "GABA.Sst", "GABA.PV", "GABA.NGF","GABA.CGE") )  


DG.dorsal.seiz<-findMarkers(sce.DG[,sce.DG$DV=='DG.dorsal'],pval.type='all',
                           direction='up',group=sce.DG[,sce.DG$DV=='DG.dorsal']$condition)
DG.ventral.seiz<-findMarkers(sce.DG[,sce.DG$DV=='DG.ventral'],pval.type='all',
                           direction='up',group=sce.DG[,sce.DG$DV=='DG.ventral']$condition)


sce.DG$annotation<-factor(
  ifelse(sce.DG$condition=='control',paste0(sce.DG$axis,'.sham'),paste0(sce.DG$axis,'.seizure')))
         
         ifelse(sce.DG$label %in% c(4,1),'DG.ventral',
                ifelse(sce.DG$label %in% c(8,2),'DG.dorsal.seizure',
                       ifelse(sce.DG$label==7, 'DG.ventral.seizure', 'DG.ventral.mixed')))))

pdf('20210924_umap_sceDG_byCondition.pdf')
plotUMAP(sce.DG,text_by='cellType',colour_by='condition')
dev.off()

jpeg('20210924_umap_sceDG_byCType.jpg')
plotUMAP(sce.DG,text_by='label',colour_by='cellType')
dev.off()

markers.custom=list(
  "Stronger dorsal upregulation"='Nptx2',
  "Stronger ventral upregulation"='Sv2c'
)
pdf('20210924_violinPlots_DG_bycellType.pdf',h=5,w=5)
for(i in 1:length(markers.custom)){
  print(
    plotExpression(sce.DG, exprs_values = "logcounts", features=markers.custom[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
                   width = 0.3) + 
      ggtitle(label=(names(markers.custom)[i])) +
      theme(axis.text.x = element_text(angle = 90))
  )}
dev.off()

hex <- make_hexbin(sce.DG, nbins = 50,
                   dimension_reduction = "UMAP", use_dims=c(1,2))
label_df <- make_hexbin_label(hex, col="cellType")
pp <- plot_hexbin_meta(hex, col="condition", action="majority",
                       xlab='UMAP1',ylab='UMAP2') + ggtitle(" ")

pdf("UMAP_hex_condition.pdf",h=5,w=6)
pp + #ggrepel::geom_label_repel(data = label_df, aes(x=x, y=y), 
  #                           colour="black",  label.size = NA, fill = NA)+
  theme(legend.position='right',text = element_text(size = 20))
dev.off()





marks.filtered<-list()
for (i in 1:length(marks)){
  marks.filtered[[i]]<-rownames(marks[[i]][marks[[i]]$FDR < 0.05 | marks[[i]]$summary.logFC > 2,])
}
for (i in 1:length(marks)){
  marks.filtered[[i]]<-na.omit(marks.filtered[[i]][1:5])
}
marks.filtered<-do.call(c,marks.filtered)

#sce.subset.ID<-sce.subset[,unique(colnames(sce.subset))]

pdf('test.pdf',height=6,width=12)
plotHeatmap(sce.subset,features=features,order_columns_by="k_50_label",
            colour_columns_by="k_50_label",show_rownames=F,cluster_rows=F,center=F)
dev.off()




save(sce.subset,hdg,file="/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50.rda")



pdf('heatmap_IEGs_bySample.pdf')
#features=c('Arc','Bdnf','Nptx2','Jun','Junb','Nr4a1','Baz1a','Rheb',')
plotGroupedHeatmap(sce.subset,features=features,cluster_rows=F,cluster_cols=F,group='sample_name',center=F)
dev.off()


markers.custom=list(
  'CA3' = 'Bves',
  'CA1' = 'Stac',
  'DG' = 'Nptx2'
)
pdf('violinPlots_byARG.pdf',h=6,w=6)
for(i in 1:length(markers.custom)){
  print(
    plotExpression(sce.sub, exprs_values = "logcounts", features=markers.custom[[i]],
                   x="a", colour_by="a", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
                   width = 0.3) + 
      ggtitle(label=(names(markers.custom)[i])) +
      theme(axis.text.x = element_text(angle = 90))
  )}
dev.off()





pdf('20210913_UMAP_byCondition.pdf',h=5,w=5)
plotUMAP(sce.subset,colour_by='condition',text_by='annotation',text_size=3,add_legend=T) 
dev.off()

pdf('20210913_UMAP_byNptx2.pdf',h=7,w=7)
plotUMAP(sce.subset,colour_by='Nptx2',text_by='annotation',text_size=3,add_legend=F) 
dev.off()

pdf('20210913_UMAP_byBdnf.pdf',h=7,w=7)
plotUMAP(sce.subset,colour_by='Bdnf',text_by='annotation',text_size=3,add_legend=F) 
dev.off()

pdf('20210913_UMAP_byStac.pdf',h=7,w=7)
plotUMAP(sce.subset,colour_by='Stac',text_by='annotation',text_size=3,add_legend=F) 
dev.off()

pdf('20210913_UMAP_byBves.pdf',h=7,w=7)
plotUMAP(sce.subset,colour_by='Bves',text_by='annotation',text_size=3,add_legend=F) 
dev.off()

pdf('20210913_UMAP_byDGcluster.pdf',h=7,w=7)
plotUMAP(sce.DG,colour_by='label',text_by='label',text_size=3,add_legend=T) 
dev.off()

pdf('20210913_UMAP_byDGcondition.pdf',h=7,w=7)
plotUMAP(sce.DG,colour_by='condition',text_by='label',text_size=3,add_legend=T) 
dev.off()

pdf('20210913_UMAP_byInf2.pdf',h=7,w=7)
plotUMAP(sce.DG,colour_by='Inf2',text_by='label',text_size=3,add_legend=T) 
dev.off()

pdf('20210913_UMAP_byTenm4.pdf',h=7,w=7)
plotUMAP(sce.DG,colour_by='Tenm4',text_by='label',text_size=3,add_legend=T) 
dev.off()

cellType<-list()

for (i in 1:length(levels(sce.DG$cellType))){
  cellType[[i]]<-sce.DG$cellType==levels(sce.DG$cellType)[i]
}
names(cellType)<-levels(sce.DG$cellType)
dimnames<-list(colnames(sce.DG),names(cellType))
cellType<-matrix(unlist(cellType),ncol=5,nrow=4472,dimnames=dimnames)
map_cellType<-cor(cellType,result)
pheatmap(map_cellType,cluster_cols=T,cluster_rows=T)

label<-list()

for (i in 1:length(levels(sce.DG$label))){
  label[[i]]<-sce.DG$label==levels(sce.DG$label)[i]
}
names(label)<-levels(sce.DG$label)
dimnames<-list(colnames(sce.DG),names(label))
label<-matrix(unlist(label),ncol=9,nrow=4472,dimnames=dimnames)
map_label<-cor(label,result)
pheatmap(map_label,cluster_cols=T,cluster_rows=T)

plotExpression(sce.DG,features=features,
               x="cellType", colour_by="DV", point_alpha=0.5, point_size=.7,add_legend=F,ncol=4)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3)+
  theme(axis.text.x = element_text(angle = 45))

sce.DG$annotation<-ifelse(sce.DG$k_20_label %in% c(1,2,8,9),'DG.ventral',
                          'DG.dorsal')
sce.DG$shared<-factor(ifelse(sce.DG$annotation %in% c('DG.sham.1','DG.seizure.1'), 'group1',
                      ifelse(sce.DG$annotation %in% c('DG.sham.2','DG.seizure.2'),'group2',
                             'group3')))

pdf('umapDG2.pdf',h=4,w=4)
plotUMAP(sce.DG,colour_by='condition') + 
  theme(legend.position='bottom',
        legend.justification='center',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(size = 10))
dev.off()

pdf('umapDG1.pdf',h=4,w=4)
plotUMAP(sce.DG,colour_by='axis') + 
  theme(legend.position='bottom',
        legend.justification='center',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(size = 10))+
  scale_color_manual(values=values)
dev.off()

pdf('DG-features.pdf',h=2.5,w=8)
plotExpression(sce.DG,features=features,
               x="annotation", colour_by="axis", point_alpha=0.5, 
               point_size=.7,add_legend=F,ncol=6)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3)+
  theme(axis.text.x = element_text(angle = 90,vjust=0.6),
        text=element_text(size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_color_manual(values=values)
dev.off()
enrichList<-list()
for(i in 1:length(marksList))
enrichList[[i]] <- enrichGO(gene = marksList[[i]],
                      OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "all",
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

enrich_go2 <- enrichGO(gene = marks2,
                      OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "CC",
                      pAdjustMethod = "BH", pvalueCutoff = 0.001, qvalueCutoff = 0.05)

enrich_go3 <- enrichGO(gene = marks3,
                       OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "CC",
                       pAdjustMethod = "BH", pvalueCutoff = 0.001, qvalueCutoff = 0.05)

colPull<-function(x,v){
pull(as.data.frame(x),var=v)
}
BP<-do.call(pull,go_BP)