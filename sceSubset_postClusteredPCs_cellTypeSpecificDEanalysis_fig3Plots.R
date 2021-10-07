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
library(edgeR)
library(org.Mm.eg.db)
library(clusterProfiler)
library(UpSetR)

load("/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50_clusteredPCs.rda")
load('/users/enelson/mouseHPC_broad_differentialExpressionAnalysis_results.rda')
load('argList.rda') #this was compiled from an excel sheet from a paper from Jesse Gray's lab, linked later in script

###Cell type-specific differential expression analysis==========================
summed$annotation.DG<-as.character(summed$annotation)
summed$annotation.DG<-factor(ifelse(summed$annotation.DG %in% c('DG.1','DG.2'),'DG',summed$annotation.DG))

de.specific <- pseudoBulkDGE(summed,
                             label=summed$annotation.DG,
                             design=~condition,
                             coef=2,
                             condition=summed$condition
)

is.de <- decideTestsPerLabel(de.specific, threshold=0.05)
summarizeTestsPerLabel(is.de)

##Pull out DEGs for 3 labels with most DE genes
is.de<-as.data.frame(is.de)
DEG.DG<-rownames(is.de)[is.de$DG %in% c(1,-1)]
DEG.CA1.2<-rownames(is.de)[is.de$CA1.2 %in% c(1,-1)]
DEG.CA3.1<-rownames(is.de)[is.de$CA3.1 %in% c(1,-1)]
DEG.overall<-rownames(x)

##Pull out DEGs that are also ARGs (as defined by Jesse Gray's group's paper
##Link is here: https://www.cell.com/neuron/fulltext/S0896-6273(18)30285-X?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS089662731830285X%3Fshowall%3Dtrue)
DEG.ARG.DG<-DEG.DG %in% arg_list
DEG.ARG.CA1.2<-DEG.CA1.2 %in% arg_list
DEG.ARG.CA3.1<-DEG.CA3.1 %in% arg_list
DEG.ARG.overall<-rownames(x) %in% arg_list

##make lists of lists
DEG.all<-list(DEG.DG,DEG.CA1.2,DEG.CA3.1,DEG.overall)
DEG.arg.all<-list(DEG.ARG.DG,DEG.ARG.CA1.2,DEG.ARG.CA3.1,DEG.ARG.overall)

##make upset lists
upset.all<-fromList(DEG.all)
upset.arg<-fromList(DEG.arg.all)

##upset plots for fig 3
pdf('upsetPlot_allDEGs.pdf')
upset(upset.all)
dev.off()

pdf('upsetPlot_ARGSonly.pdf')
upset(upset.arg)
dev.off()


##make new SCE object for CA3.D,CA1.V,and DG for the violin plots
sce.fig3<-sce.subset[,sce.subset$cellType %in% c('CA1.V','CA3.D','DG.control','DG.seizure')]

##set up labels for plots
sce.fig3$cellType<-droplevels(sce.fig3$cellType)
sce.fig3$cellType<-as.character(sce.fig3$cellType)
sce.fig3$cellType<-factor(
  ifelse(sce.fig3$condition=='seizure' & sce.fig3$cellType %in% c('CA1.V','CA3.D'), 
         paste0(sce.fig3,'.seizure'),sce.fig3$cellType)
)

##Violin plots
pdf('20210922_seizure_violin_CA1_DG.pdf')
plotExpression(sce.fig3,features=c("Stac"),
               x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3)+
  ggtitle('Upregulated in CA1.2 & DG only')
dev.off()

pdf('20210922_seizure_violin_all.pdf')
plotExpression(sce.fig3,features=c("Fam126b"),
               x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3)+
  ggtitle('Upregulated in all')
dev.off()

pdf('20210922_seizure_violin_DG_CA3.pdf')
plotExpression(sce.fig3,features=c( "Bdnf"),
               x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3) +
  ggtitle('Upregulated in CA3.1 & DG only')
dev.off()

pdf('20210922_seizure_violin_DGOnly.pdf')
plotExpression(sce.fig3,features=c("Arc"),
               x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3) +
  ggtitle('Upregulated in DG only')
dev.off()

pdf('20210922_seizure_violin_CA_only.pdf')
plotExpression(sce.fig3,features=c("Mfap3l"),
               x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3)+
  ggtitle('Upregulated in CA1.2 & CA3.1 only')
dev.off()

save(de.specific,is.de,DEG.all,DEG.arg.all,sce.fig3,file='HPC_fig3Objects.rda')
