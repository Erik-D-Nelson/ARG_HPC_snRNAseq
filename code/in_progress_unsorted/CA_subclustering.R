##Feature selection for CA
sce.CA<-devianceFeatureSelection(sce.CA, assay="counts", 
                                   fam="poisson",sorted=TRUE)
plot(rowData(sce.CA)$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,50000))

hdg.CA<-rownames(counts(sce.CA))[1:5000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce.CA<-nullResiduals(sce.CA, assay = "counts",
                        fam = "poisson", type = "pearson")

##Dimensionality reduction for CA
set.seed(1000) 
sce.CA <- runPCA(sce.CA,exprs_values="poisson_pearson_residuals", subset_row=hdg.CA) 
reducedDimNames(sce.CA)

set.seed(100000)
sce.CA <- runUMAP(sce.CA, dimred="PCA",min_dist=.5)

g <- buildSNNGraph(sce.CA, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.CA)$k_20_label <- factor(clust)
table(sce.CA$k_20_label)

g <- buildSNNGraph(sce.CA, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.CA)$k_10_label <- factor(clust)
table(sce.CA$k_10_label)

g <- buildSNNGraph(sce.CA, k=50, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.CA)$k_50_label <- factor(clust)
table(sce.CA$k_50_label)

set.seed(1000)
clusters <- quickCluster(sce.CA)
sce.CA <- computeSumFactors(sce.CA, cluster=clusters)
sce.CA <- logNormCounts(sce.CA)

markers.CA <- findMarkers(sce.CA, pval.type="all", direction="up")
