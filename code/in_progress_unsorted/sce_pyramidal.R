###Subclustering, Pyramidal==========================================================
##ok cool, let's subcluster some stuff
sce.Pyramidal<-sce.subset[,sce.subset$level_6=="Pyramidal"]

##Feature selection for Pyramidal
sce.Pyramidal<-devianceFeatureSelection(sce.Pyramidal, assay="counts", 
                                   fam="poisson",sorted=TRUE)
plot(rowData(sce.Pyramidal)$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,50000))

hdg.Pyramidal<-rownames(counts(sce.Pyramidal))[1:5000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce.Pyramidal<-nullResiduals(sce.Pyramidal, assay = "counts",
                        fam = "poisson", type = "pearson")

##Dimensionality reduction for Pyramidal
set.seed(1000) 
sce.Pyramidal <- runPCA(sce.Pyramidal,exprs_values="poisson_pearson_residuals", subset_row=hdg.Pyramidal) 
reducedDimNames(sce.Pyramidal)

set.seed(100000)
sce.Pyramidal <- runUMAP(sce.Pyramidal, dimred="PCA",n_neighbors=20,min_dist=0.2)

g <- buildSNNGraph(sce.Pyramidal, k=50, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.Pyramidal)$label <- factor(clust)
table(sce.Pyramidal$label)

g <- buildSNNGraph(sce.Pyramidal, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.Pyramidal)$k_20_label <- factor(clust)
table(sce.Pyramidal$k_20_label)

g <- buildSNNGraph(sce.Pyramidal, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.Pyramidal)$k_10_label <- factor(clust)
table(sce.Pyramidal$k_10_label)

plotUMAP(sce.Pyramidal,colour_by='label',text_by='label')

library(scDblFinder)
dbl.out <- findDoubletClusters(sce.Pyramidal)
dbl.out

chosen.doublet <- rownames(dbl.out)[isOutlier(dbl.out$num.de, type="lower", log=TRUE)]
chosen.doublet

markers.Pyramidal <- findMarkers(sce.Pyramidal, pval.type="all", direction="up")

sce.Pyramidal$pyr_cell_class<-factor(
  ifelse(sce.Pyramidal$label==15, "Mossy",
      ifelse(sce.Pyramidal$label==1, "CA3_D",
          ifelse(sce.Pyramidal$label==11, "CA2",
                ifelse(sce.Pyramidal$label==7, "CA3_V",
                    ifelse(sce.Pyramidal$label==5, "CA1_D",
                          ifelse(sce.Pyramidal$label==6, "CA1_V",
                             ifelse(sce.Pyramidal$label==10, "CA1_V_ProS",
                              ifelse(sce.Pyramidal$label==13, "ProS_Sub",
                                      ifelse(sce.Pyramidal$label %in% c(4,8,12,17,18), 'L5-6_RHP',
                                             'L2/3-4_RHP'))))))))))


summed<-sce.subset[,sce.subset$level_6=='Pyramidal']
summed <- aggregateAcrossCells(summed, 
                               ids=colData(summed)$sample_name)
location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(summed)$GENEID, 
                   column="SEQNAME", keytype="GENEID")
is.mito <- which(location=="MT")
stats <- perCellQCMetrics(summed, subsets=list(Mito=is.mito))
colnames(stats)
colData(summed)$sum<-stats$sum
summed<-summed[!is.na(rowData(summed)$length),]

#Run cqn
cqn.summed <- cqn(counts(summed), x=rowData(summed)$gc_content,lengths=rowData(summed)$length,
                  sizeFactors = summed$sum,
                  verbose = TRUE)
library(edgeR)
y <- DGEList(counts(summed), samples=colData(summed),
             genes=rowData(summed))


design <- model.matrix(~condition, y$samples)
y$offset <- cqn.summed$glm.offset
summed.cqn <- estimateGLMCommonDisp(y, design = design)

efit <- glmFit(summed.cqn, design = design)
elrt <- glmLRT(efit, coef = 2)
x2<-topTags(elrt,n=10000)
x2<-as.data.frame(x2)
#x2<-x2[rownames(x2) %in% arg_list,15:19]
x2<-x2[x2$FDR<=0.05,]
plotSmear(elrt, de.tags = rownames(x2))


pyr<-topTags(elrt, n = 10000)
pyr<-as.data.frame(pyr)
pyr<-pyr[pyr$FDR<=0.05,]
x2[,14:18]

pyramidal_sub_list<-list(
  'CA1+subiculum, all'= c('Fibcd1','Pex5l','Adgrl2'),
  'CA1, dorsal' = c('Zdhhc2','Scn3b'),
  'CA1, ventral' = c('Dcn','Nov','Adra1a'),
  'Prosubiculum' = c('Ndst4'),
  'Subiculum' = c('Fn1','Nts','Slc17a6'),
  'CA2' = c('Ntf3','Cacng5','Adcy5','Ccdc3'),
  'CA3, all' = c('Galnt3','Adamts9','Nectin3'),
  'CA3, dorsal' = c('Cdh24','Ephb1','Chgb'),
  'CA3, ventral' = c('St18','Htr2c'),
  'Mossy cell' = c('Calb2','Csf2rb2')#,
#  'Layer 2/3 cortex' = c('Cux1' , 'Cux2' , 'Satb2','Mef2c'),
#  'Layer 5/6' = c('Foxp2','Tle4')
  )

pdf("20200222_Pyramidal_sub_markers3.pdf", height=6,width=15)
for(i in 1:length(pyramidal_sub_list)){
  print(
    plotExpression(sce.Pyramidal,features=c('Gpc3'),
                   x="k_20_label", colour_by="pyr_cell_class", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(pyramidal_sub_list)[i], " markers"))
  )}
dev.off()

###Pyr subclustering============================================================
sce.Pyramidal$level_7<-factor(
  ifelse(sce.Pyramidal$label %in% c(1,5,6,7,10,11,13,15), 
         "CA/Sub/ProS",
              "RHP/Cortex"))

sce.CA<-sce.Pyramidal[,sce.Pyramidal$level_7=='CA/Sub/ProS',]
sce.RHP<-sce.Pyramidal[,sce.Pyramidal$level_7=='RHP/Cortex',]