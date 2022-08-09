###DG subclustering and comparison, revisited
### Normalization===============================================================
##Just to get logcounts assay
##Don't actually use for dimension reduction or feature selection
set.seed(13700)
clusters <- quickCluster(DG.control)
DG.control <- computeSumFactors(DG.control, cluster=clusters)
DG.control <- logNormCounts(DG.control)

### Dimensionality reduction/feature selection===============================================
#This is based on fitting Poisson model and then calculating Pearson residuals
#Let's plot highly deviant genes first
DG.control<-devianceFeatureSelection(DG.control, assay="counts", 
                                        fam="poisson",sorted=TRUE)
plot(rowData(DG.control)$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,100000))

##take the top 10000 highly deviant genes
hdg<-rownames(counts(DG.control))[1:5000]

#Run nullResiduals to calculate Pearson residuals from model
DG.control<-nullResiduals(DG.control, assay = "counts",
                             fam = "poisson", type = "pearson")

##Now, PCA on residuals to approximate GLM-PCA
set.seed(1000) 
DG.control <- runPCA(DG.control,exprs_values="poisson_pearson_residuals",ncomponents=50, 
                        subset_row=hdg,BSPARAM=BiocSingular::RandomParam()) 
reducedDimNames(DG.control)

##Plot UMAP
set.seed(100000)
DG.control <- runUMAP(DG.control, dimred="PCA",subset_row=hdg)

#Verify length
ncol(reducedDim(DG.control, "PCA"))

# Save into a new data file, which will dedicate for downstream [prelim] analysis
save(DG.control, hdg, DG.seizure,
     file="20210830_DG_sce_objects.rda")

### Clustering=============================================================
g20 <- buildSNNGraph(DG.control, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g20)$membership
colData(DG.control)$label <- factor(clust)
table(colLabels(DG.control))


### Normalization===============================================================
##Just to get logcounts assay
##Don't actually use for dimension reduction or feature selection
set.seed(13700)
clusters <- quickCluster(DG.seizure)
DG.seizure <- computeSumFactors(DG.seizure, cluster=clusters)
DG.seizure <- logNormCounts(DG.seizure)

### Dimensionality reduction/feature selection===============================================
#This is based on fitting Poisson model and then calculating Pearson residuals
#Let's plot highly deviant genes first
DG.seizure<-devianceFeatureSelection(DG.seizure, assay="counts", 
                                     fam="poisson",sorted=TRUE)
plot(rowData(DG.seizure)$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,100000))

##take the top 10000 highly deviant genes
hdg<-rownames(counts(DG.seizure))[1:5000]

#Run nullResiduals to calculate Pearson residuals from model
DG.seizure<-nullResiduals(DG.seizure, assay = "counts",
                          fam = "poisson", type = "pearson")

##Now, PCA on residuals to approximate GLM-PCA
set.seed(1000) 
DG.seizure <- runPCA(DG.seizure,exprs_values="poisson_pearson_residuals",ncomponents=50, 
                     subset_row=hdg,BSPARAM=BiocSingular::RandomParam()) 
reducedDimNames(DG.seizure)

##Plot UMAP
set.seed(100000)
DG.seizure <- runUMAP(DG.seizure, dimred="PCA",subset_row=hdg)

#Verify length
ncol(reducedDim(DG.seizure, "PCA"))

# Save into a new data file, which will dedicate for downstream [prelim] analysis
save(DG.control, DG.seizure,
     file="20210830_DG_sce_objects.rda")

### Clustering=============================================================
g20 <- buildSNNGraph(DG.seizure, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g20)$membership
colData(DG.seizure)$label <- factor(clust)
table(colLabels(DG.seizure))

#Subset for DG only and aggregateAcrossCells to pseudobulk
sce.DG<-sce.subset[,sce.subset$label %in% c('5','8')]
summed <- aggregateAcrossCells(sce.DG, 
                               ids=colData(sce.DG)[c('condition','sample_name')])

#Use perCellQCMetrics to get library sizes for each sample 
#(didn't carry over after pseudobulk, may be an easier way but this worked
#so I ran with it)
location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(summed)$GENEID, 
                   column="SEQNAME", keytype="GENEID")
is.mito <- which(location=="MT")
stats <- perCellQCMetrics(summed, subsets=list(Mito=is.mito))
colnames(stats)
colData(summed)$sum<-stats$sum

#Filter out lowly expressed genes and remove genes with length=NA
summed<-summed[filterByExpr(summed, group=summed$condition),]

#make DGElist
y <- DGEList(counts(summed), samples=colData(summed),
             genes=rowData(summed))
y$samples$condition<-factor(y$samples$condition)

#make design matrix
design <- model.matrix(~condition, y$samples)

###GLMfit, no length normalization===============================================================


#normalization 
y <- calcNormFactors(y)

#estimate dispersion parameter
y <- estimateGLMCommonDisp(y, design = design)

#run DE analysis and look at results
efit <- glmFit(y, design = design)
elrt <- glmLRT(efit, coef = 2)
res<-topTags(elrt,n=Inf)
x<-as.data.frame(res)
x<-x[x$FDR<=0.05,]

set.seed(1000101001)
mnn.out <- fastMNN(sce.DG, batch=sce.DG$Sample, d=50, k=200, subset.row=hdg,
                   assay.type="logcounts",
                   BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
mnn.out


plotExpression(sce.DG, exprs_values = "logcounts", features=c('Ntng1', 'Rasgrp1', 'Nrxn1', 
                                                              'Tenm2', 'Tenm4', 'Cntn4', 
                                                              'Trim2', 'Slc4a4', 'Mast4'),
               x="label", colour_by="label", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
               width = 0.3)

features=c('Ntng1', 'Rasgrp1', 'Nrxn1', 
           'Tenm2', 'Tenm4', 'Cntn4', 
           'Trim2', 'Slc4a4', 'Mast4')