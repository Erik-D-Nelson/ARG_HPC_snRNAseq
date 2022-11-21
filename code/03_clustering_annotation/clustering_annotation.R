library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(schex)
library(viridis)

load("processed_data/sce_subset.rda")


### Clustering ============================================================
g25 <- buildSNNGraph(sce.subset, k=25, use.dimred = 'PCA')
clust25 <- igraph::cluster_walktrap(g25)$membership
colData(sce.subset)$k_25_label <- factor(clust20)
table(colData(sce.subset)$k_25_label)

###Annotation=============================================================
annotationTab<-data.frame("cluster"=c(1:28),"annotation"=NA)
annotationTab$annotation[c(26,15,9)] <- paste0("GABA.", c("3","4","5"))
annotationTab$annotation[c(21,8)] <- paste0("GABA.", c("1","2"))
annotationTab$annotation[c(6,4)]<-paste0("GC.", c("1","2"))
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
markers<-findMarkers(sce.subset,pval.type='all',
                  direction='up',group=sce.subset$annotation)

save(markers,file='processed_data/sceSubset_markerGenes.rda')

##make cellType factor for plots
sce.subset$cellType<-as.character(sce.subset$annotation)
sce.subset$cellType<-ifelse(sce.subset$k_25_label %in% c(1,27,24,10,14,23,16),'L2/3',
         ifelse(sce.subset$k_25_label %in% c(20,25),'L5/Po',
                ifelse(sce.subset$k_25_label %in% c(28,12,7,19),'L6/6b',
                                sce.subset$cellType)))
                                                                        
sce.subset$cellType<-factor(sce.subset$cellType,
                            levels=c('GC.1','GC.2','CA4','CA3.1','CA3.2','CA2',
                                     'CA1','PS.1','PS.2','Sub',
                                     'L2/3','L5/Po','L6/6b','GABA.1',
                                     'GABA.2','GABA.3','GABA.4','GABA.5')
)

###Figure 1 plots===============================================================
##Hexbin plots, by cell type and condition
##Make new sce with hexbins using make_hexbin
hex <- make_hexbin(sce.subset, nbins = 100,
                   dimension_reduction = "UMAP", use_dims=c(1,2))

##Make cluster labels for plot
label_df <- make_hexbin_label(hex, col="cellType")

##custom palette using turbo 
colors<-viridis_pal(1,0,1,1,'H')(18)

#make plot 
pp <- plot_hexbin_meta(hex, col="cellType", action="majority",
                xlab='UMAP1',ylab='UMAP2',colors=colors) + ggtitle("") 

pdf("UMAP_hex_cellType.pdf",h=5,w=6)
pp + theme(legend.position='right',text = element_text(size = 14)) + 
        ggrepel::geom_label_repel(data = label_df, aes(x=x, y=y, label = label), 
                                 colour="black",  label.size = NA, fill = NA)  
dev.off()


#make plot 
pp2 <- plot_hexbin_meta(hex, col="condition", action="majority",
                       xlab='UMAP1',ylab='UMAP2',colors=c('cornflowerblue','orange')) + ggtitle("") 
pdf("UMAP_hex_condition.pdf",h=5,w=6)
pp2 +theme(legend.position='right',text = element_text(size = 14))
dev.off()

##dotplot for figure 1
features=c('Slc17a7','Gad2','Prox1',
           'Calb2','Nectin3','St18',
           'Ntsr2','Man1a','Dcn','Fras1','Ndst4','Fn1','Satb2',
           'Pamr1','Foxp2','Ndnf','Htr3a',
           'Lhx6','Lamp5','Pvalb','Sst')

#pdf('dotPlot_clusterMarkers.pdf',h=6.5,w=12)
plotDots(sce.subset,group='cellType',features=features,color=c('white','red')) + 
  scale_y_discrete(limits=rev(features)) +
  theme(axis.text.x = element_text(angle = 45,vjust=0.75),text = element_text(size = 20))
#dev.off()


##Make plots for fig s2
features=c('Zbtb20','Satb2',
           'Ntng2','Cux1',
           'Dcn', 'Cux2',
           'Nos1', 'Pamr1',
           'Fras1','Sema5b',
           'Ndst4','Foxp2',
           'Cntn6','Cobll1',
           'Fn1', 'Ctgf')
pdf('plots/figS2/figS2.pdf',w=8,h=11)
plotExpression(sce.subset,features=features,
               x="annotation", colour_by="annotation", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3)+
  theme(axis.text.x = element_text(angle = 45))
dev.off()

save(sce.subset, hdg, pc.choice.hpc,marks,file="sce_subset.rda")
##Proceed to figure 2 plots
