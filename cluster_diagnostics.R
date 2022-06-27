##This script will do the following
##1) Take original sce.subset object, rerun feature selection and dimred such
##   that we have 1500 HDGs and 100 PCs. Keep annotation col.
##2) Run getClusteredPCs
##3) Take PC choice from getClusteredPCs and run PCA again
##4) Use multiple values of k for graph-based clustering with new PCA
##5) Spit out UMAPs for all clustering variation
##6) Run two sets of cluster diagnostics

##Code for running getClusteredPCs, sce.subset w/1500 HDGs
##Gonna do this in a PC space of 100 PCs


library(SingleCellExperiment)
library(scRNAseq)
library(scater)
library(scran)
library(jaffelab)
library(SingleCellExperiment)
library(scRNAseq)
library(scater)
library(scran)
library(scry)
library(jaffelab)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(bluster)
library(pheatmap)

# ===


load("/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50_clusteredPCs.rda",verbose=T)

##make sure we have annotation
annotationTab<-data.frame("cluster"=c(1:23),"annotation"=NA)
annotationTab$annotation[c(11,21,8,13)] <- paste0("GABA.", c("1","2","3","4"))
annotationTab$annotation[c(22,17,23,6,19,18)]<-("RHP.L2/3.")# , c("1","2","3","4","5","6"))
annotationTab$annotation[c(9,16,15,1)]<-paste0("RHP.L5/6/6b.", c("1","2","3","4"))
annotationTab$annotation[c(2,14)]<-paste0("CA3.", c("1","2"))
annotationTab$annotation[c(5,3)]<-paste0("DG.", c("1","2"))
annotationTab$annotation[c(7,4,10,12,20)]<-c("CA1.1","CA1.2","Sub","CA2","CA4")
annotationTab$cluster<-factor(annotationTab$cluster)
annotationTab$annotation<-factor(annotationTab$annotation)
sce.subset$annotation <- annotationTab$annotation[match(sce.subset$k_50_label,
                                                        annotationTab$cluster)]

##################################set HDGs and run PCA=========================
hdg<-rownames(counts(sce.subset))[1:4000]
set.seed(1000) 
sce.subset <- runPCA(sce.subset,exprs_values="poisson_pearson_residuals", 
                    subset_row=hdg,ncomponents=50) 
reducedDimNames(sce.subset)

set.seed(100000)
sce.subset <- runUMAP(sce.subset, dimred="PCA",subset_row=hdg)

#Verify length
ncol(reducedDim(sce.subset, "PCA"))

# getClusteredPCs() to identify working PC space
pc.choice.hpc <- getClusteredPCs(reducedDim(sce.subset, "PCA"))

# How many PCs should use in this space?
chosen<-metadata(pc.choice.hpc)$chosen
chosen

## Plot n Clusters vs. d PCs
pdf("plots/newSceSubset_getClusteredPCs.pdf")
plot(pc.choice.hpc$n.pcs, pc.choice.hpc$n.clusters,
     main="PC Choice via getClusteredPCs()",
     abline(v=metadata(pc.choice.hpc)$chosen, col="red", lty="dashed", lwd=0.8)) 
dev.off()


##run PCA again
set.seed(1000) 
sce.subset <- runPCA(sce.subset,exprs_values="poisson_pearson_residuals", 
                    subset_row=hdg,ncomponents=50) 
reducedDimNames(sce.subset)


##run UMAP
set.seed(100000)
sce.subset <- runUMAP(sce.subset, dimred="PCA",subset_row=hdg)



# Save
save(sce.subset, hdg, pc.choice.hpc,
     file="/users/enelson/newSceSubset_postClusterdPCs.rda")

##################################clustering for new PC choice=========================
g40 <- buildSNNGraph(sce.subset, k=40, use.dimred = 'PCA')
clust40 <- igraph::cluster_walktrap(g40)$membership
colData(sce.subset)$k_40_label <- factor(clust40)
table(colData(sce.subset)$k_40_label)

g15 <- buildSNNGraph(sce.subset, k=15, use.dimred = 'PCA')
clust15 <- igraph::cluster_walktrap(g15)$membership
colData(sce.subset)$c_15_label <- factor(clust15)
table(colData(sce.subset)$k_15_label)

g30 <- buildSNNGraph(sce.subset, k=30, use.dimred = 'PCA')
clust30 <- igraph::cluster_walktrap(g30)$membership
colData(sce.subset)$k_30_label <- factor(clust30)
table(colData(sce.subset)$k_30_label)

g20 <- buildSNNGraph(sce.subset, k=20, use.dimred = 'PCA')
clust20 <- igraph::cluster_walktrap(g20)$membership
colData(sce.subset)$k_20_label <- factor(clust20)
table(colData(sce.subset)$k_20_label)

g25 <- buildSNNGraph(sce.subset, k=25, use.dimred = 'PCA')
clust25 <- igraph::cluster_walktrap(g25)$membership
colData(sce.subset)$k_25_label <- factor(clust25)
table(colData(sce.subset)$k_25_label)


##plot UMAPs
pdf('k50_UMAP.pdf',h=7,w=7)
plotUMAP(sce.subset,colour_by='k_50_label',
         text_by='k_50_label',text_size=3,add_legend=F) 
dev.off()

pdf('k40_UMAP.pdf',h=7,w=7)
plotUMAP(sce.subset,colour_by='k_40_label',
         text_by='k_40_label',text_size=3,add_legend=F) 
dev.off()

pdf('k30_UMAP.pdf',h=7,w=7)
plotUMAP(sce.subset,colour_by='k_30_label',
         text_by='k_30_label',text_size=3,add_legend=F) 
dev.off()

pdf('k20_UMAP.pdf',h=7,w=7)
plotUMAP(sce.subset,colour_by='k_20_label',
         text_by='k_20_label',text_size=3,add_legend=F) 
dev.off()

save(sce.subset,file='sce_subset_for_diagnostics.rda')


#######cluster diagnostics!=====================================
##silhouettes
library(bluster)
#sil50 <- as.data.frame(approxSilhouette(reducedDim(sce.subset, "PCA"), clusters=sce.subset$k_50_label))
#sil40 <- as.data.frame(approxSilhouette(reducedDim(sce.subset, "PCA"), clusters=sce.subset$k_40_label))
sil30 <- as.data.frame(approxSilhouette(reducedDim(sce.subset, "PCA"), clusters=sce.subset$k_30_label))
sil25 <- as.data.frame(approxSilhouette(reducedDim(sce.subset, "PCA"), clusters=sce.subset$k_25_label))
sil15 <- as.data.frame(approxSilhouette(reducedDim(sce.subset, "PCA"), clusters=sce.subset$k_15_label))
sil20 <- as.data.frame(approxSilhouette(reducedDim(sce.subset, "PCA"), clusters=sce.subset$k_20_label))

sil30$closest <- factor(ifelse(sil30$width > 0, sce.subset$k_30_label, sil30$other))
sil30$cluster <- sce.subset$k_30_label
sil25$closest <- factor(ifelse(sil25$width > 0, sce.subset$k_25_label, sil25$other))
sil25$cluster <- sce.subset$k_25_label
sil20$closest <- factor(ifelse(sil20$width > 0, sce.subset$k_20_label, sil20$other))
sil20$cluster <- sce.subset$k_20_label
sil15$closest <- factor(ifelse(sil15$width > 0, sce.subset$k_15_label, sil15$other))
sil15$cluster <- sce.subset$k_15_label


library(ggplot2)
pdf('plots/new2_sil25.pdf')
ggplot(sil25, aes(x=cluster, y=width, colour=closest)) +
  ggbeeswarm::geom_quasirandom(method="smiley")      
dev.off()

pdf('plots/new2_sil15.pdf')
ggplot(sil15, aes(x=cluster, y=width, colour=closest)) +
  ggbeeswarm::geom_quasirandom(method="smiley")      
dev.off()

pdf('plots/new2_sil30.pdf')
ggplot(sil30, aes(x=cluster, y=width, colour=closest)) +
  ggbeeswarm::geom_quasirandom(method="smiley")      
dev.off()

pdf('plots/new2_sil20.pdf')
ggplot(sil20, aes(x=cluster, y=width, colour=closest)) +
  ggbeeswarm::geom_quasirandom(method="smiley")      
dev.off()

pdf('sil_annotation.pdf')
ggplot(sil.annotation, aes(x=cluster, y=width, colour=closest)) +
  ggbeeswarm::geom_quasirandom(method="smiley")  + theme(axis.text.x = element_text(angle = 45))    
dev.off()

##graph

combined <- cbind(k.annotation=sce.subset$annotation, 
                  k.50=sce.subset$k_50_label, 
                  k.40=sce.subset$k_40_label,
                  k.30=sce.subset$k_30_label,
                  k.20=sce.subset$k_20_label,)


library(clustree)
set.seed(1111)
pdf('clustree.pdf')
clustree(combined, prefix="k.", edge_arrow=FALSE)

save(sce.subset,file="/users/enelson/newSceSubset_postClusterdPCs.rda")


rm(list=ls())
sessionInfo()