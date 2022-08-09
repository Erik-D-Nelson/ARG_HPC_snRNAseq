##Feature selection for RHP
sce.RHP<-devianceFeatureSelection(sce.RHP, assay="counts", 
                                   fam="poisson",sorted=TRUE)
plot(rowData(sce.RHP)$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,50000))

hdg.RHP<-rownames(counts(sce.RHP))[1:5000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce.RHP<-nullResiduals(sce.RHP, assay = "counts",
                        fam = "poisson", type = "pearson")

##Dimensionality reduction for RHP
set.seed(1000) 
sce.RHP <- runPCA(sce.RHP,exprs_values="poisson_pearson_residuals", subset_row=hdg.RHP) 
reducedDimNames(sce.RHP)

set.seed(100000)
sce.RHP <- runUMAP(sce.RHP, dimred="PCA",subset_row=hdg.RHP,n_neighbors=20,min_dist=0.1)

g <- buildSNNGraph(sce.RHP, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.RHP)$k_20_label <- factor(clust)
table(sce.RHP$k_20_label)

g <- buildSNNGraph(sce.RHP, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.RHP)$k_10_label <- factor(clust)
table(sce.RHP$k_10_label)


set.seed(1000)
clusters <- quickCluster(sce.RHP)
sce.RHP <- computeSumFactors(sce.RHP, cluster=clusters)
sce.RHP <- logNormCounts(sce.RHP)

markers.RHP <- findMarkers(sce.RHP, groups=k_50_label,pval.type="all", direction="up")


g <- buildSNNGraph(sce.RHP, k=50, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.RHP)$k_50_label <- factor(clust)
table(sce.RHP$k_50_label)