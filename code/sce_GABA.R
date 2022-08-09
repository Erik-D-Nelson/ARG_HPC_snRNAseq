sce.subset$level_6<-as.character(sce.subset$level_5)
sce.subset$level_6<-factor(
  ifelse(sce.subset$label %in% c(2,9,10), "GABA",
        ifelse(sce.subset$label %in% c(5,8), "DG",  
            "Pyramidal")))

markers.subset<-findMarkers(sce.GABA,groups='cell_type',pval.type='all',direction='up')

###Subclustering, GABA==========================================================
##ok cool, let's subcluster some stuff
sce.GABA<-sce.subset[,sce.subset$level_6=="GABA"]

##Feature selection for GABA
sce.GABA<-devianceFeatureSelection(sce.GABA, assay="counts", 
                                   fam="poisson",sorted=TRUE)
plot(rowData(sce.GABA)$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,50000))

hdg.GABA<-rownames(counts(sce.GABA))[1:5000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce.GABA<-nullResiduals(sce.GABA, assay = "counts",
                        fam = "poisson", type = "pearson")

##Dimensionality reduction for GABA
set.seed(1000) 
sce.GABA <- runPCA(sce.GABA,exprs_values="poisson_pearson_residuals", subset_row=hdg.GABA) 
reducedDimNames(sce.GABA)

set.seed(100000)
sce.GABA <- runUMAP(sce.GABA, dimred="PCA",subset_row=hdg.GABA,n_neighbors=20,min_dist=0.1)

g <- buildSNNGraph(sce.GABA, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.GABA)$k_20_label <- factor(clust)
table(sce.GABA$k_20_label)

g <- buildSNNGraph(sce.GABA, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.GABA)$k_10_label <- factor(clust)
table(sce.GABA$k_10_label)

markers.GABA <- findMarkers(sce.GABA, pval.type="all", direction="up")

##Let's plot some markers
GABA_sub_list<-list(
  'PV basket' = c('Pvalb' , 'Tac1' , 'Satb1' , 'Syt2'),
  'PV bistratified' = c('Pvalb','Sst' , 'Npy'),
  'PV axo-axonal' = c('Pvalb' , 'C1ql1' ,'Tacr1' ,  'Pthlh'),
  'Sst' = c('Sst' , 'Grm1', 'Elfn1'),
#  'Sst/Npy/Reln' = c('Sst' , 'Npy' , 'Reln' , 
 #                    'Calb1' , 'Grm1'),
  #'Sst O-LM' = c('Sst' , 'Calb1' , 'Reln',
   #              'Elfn1' , 'Serpine2', 'Npas1', 
    #             'Npas3' , 'Myh8'),
  'CGE Neurogliaform/Ivy' = c('Cacna2d1' , 'Lamp5' , 'Ndnf'),
  'MGE Neurogliaform/Ivy' = c('Cacna2d1', 'Lamp5', 'Nos1', 'Lhx6'),
  'CGE long-range projection' = c('Ntng1'),
  'CCK' = c('Cck' , 'Cnr1'),
  'VIP/Interneuron-selective-interneurons' = c('Penk' , 'Calb2' , 'Vip')
)

pdf("20200222_GABA_sub_markers.pdf", height=6,width=12)
for(i in 1:length(GABA_sub_list)){
  print(
    plotExpression(sce.GABA,features=c(GABA_sub_list[[i]]),
                   x="label", colour_by="label", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(GABA_sub_list)[i], " markers"))
  )}
dev.off()


library(BiocSingular)
library(scDblFinder)
set.seed(100)

dbl.dens <- computeDoubletDensity(sce.GABA, subset.row=hdg.GABA, 
                                  d=50)
summary(dbl.dens)
sce.GABA$doubletScore <- dbl.dens
plotUMAP(sce.GABA, colour_by="doubletScore", text_by="k_10_label")

chosen.doublet <- colnames(sce.GABA)[isOutlier(sce.GABA$doubletScore, type="high", log=TRUE)]
chosen.doublet

library(scDblFinder)
dbl.out <- findDoubletClusters(sce.GABA)
dbl.out

chosen.doublet <- rownames(dbl.out)[isOutlier(dbl.out$num.de, type="lower", log=TRUE)]
chosen.doublet

g <- buildSNNGraph(sce.GABA, k=5, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.GABA)$k_5_label <- factor(clust)
table(sce.GABA$k_5_label)

summed<-sce.GABA
summed <- aggregateAcrossCells(summed, 
                               ids=colData(summed)$sample_name)

library(edgeR)
y <- DGEList(counts(summed), samples=colData(summed),
             genes=rowData(summed))
y <- y[filterByExpr(y, group=y$samples$condition),]
y <- calcNormFactors(y)

design <- model.matrix(~condition, y$samples)
v <- voomWithQualityWeights(y, design)
fit <- lmFit(v)
fit <- eBayes(fit, robust=TRUE)

res <- topTable(fit, sort.by="p", n=Inf)
res[1:25,14:19]


sce.GABA$cell_type<-factor(
  ifelse(sce.GABA$label==3, "PV_basket",
         ifelse(sce.GABA$label==7, "PV_bistrat",
                ifelse(sce.GABA$label==11, "PV_axonal",
                       ifelse(sce.GABA$label==4, "Sst_1",
                              ifelse(sce.GABA$label==10, "Sst_2",
                                     ifelse(sce.GABA$label==8, "MGE_NGF/Ivy",
                                            ifelse(sce.GABA$label==5, "CGE_NGF",
                                                   ifelse(sce.GABA$label==9, "Projection",
                                                          ifelse(sce.GABA$label==1, "Cck_1",
                                                                 ifelse(sce.GABA$label==2, "Cck_2",
                                                                      "VIP")))))))))))