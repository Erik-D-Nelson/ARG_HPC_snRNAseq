sce.CA<-sce.subset[,sce.subset$label %in% c(3,9,8,4,11,15,19)]
sce.CA<-devianceFeatureSelection(sce.CA, assay="counts", 
                                     fam="poisson",sorted=TRUE)

#pdf('rankedGenes_vs_poissonDeviance_sceSubset.pdf')
plot(rowData(sce.CA)$poisson_deviance, type="l", 
     xlab="ranked genes",
     ylab="poisson deviance", 
     main="Feature Selection with Deviance",
     ylim=c(0,50000))
#dev.off()

#hdg<-rownames(counts(sce.CA))[1:5000]
hdg<-rownames(counts(sce.CA))[1:10000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce.CA<-nullResiduals(sce.CA, assay = "counts",
                          fam = "poisson", type = "pearson")

save(sce.CA,hdg,file="/users/enelson/20210907_sceCA_preClustering_preDimRed.rda")


### Dimensionality Reduction ============================================================
set.seed(1000) 
sce.CA <- runPCA(sce.CA,exprs_values="poisson_pearson_residuals", 
                     subset_row=hdg, BSPARAM=BiocSingular::RandomParam(),ncomponents=50) 
reducedDimNames(sce.CA)

#percent.var <- attr(reducedDim(sce.CA), "percentVar")
#chosen.elbow <- PCAtools::findElbowPoint(percent.var)
#chosen.elbow

#plot(percent.var, xlab="PC", ylab="Variance explained (%)")
#abline(v=chosen.elbow, col="red")

set.seed(100000)
sce.CA <- runUMAP(sce.CA, dimred="PCA")

#Verify length
ncol(reducedDim(sce.CA, "PCA"))


g20 <- buildSNNGraph(sce.CA, k=20, use.dimred = 'PCA')
clust20 <- igraph::cluster_walktrap(g20)$membership
colData(sce.CA)$k_20_label <- factor(clust20)
table(colData(sce.CA)$k_20_label)

g50 <- buildSNNGraph(sce.CA, k=50, use.dimred = 'PCA')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(sce.CA)$k_50_label <- factor(clust50)
table(colData(sce.CA)$k_50_label)
