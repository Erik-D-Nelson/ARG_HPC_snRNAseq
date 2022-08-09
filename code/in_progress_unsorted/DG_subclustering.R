###Subclustering, seizure and sham separate=====================================
sce.DG.sham<-sce.DG[,sce.DG$condition=="control"]
sce.DG.seizure<-sce.DG[,sce.DG$condition=='seizure']
DG.list<-list(
  "sham"=sce.DG.sham,
  "seizure"=sce.DG.seizure
)

##Feature selection for DG
for(i in 1:length(DG.list)){
DG.list[[i]]<-devianceFeatureSelection(DG.list[[i]], assay="counts", 
                                 fam="poisson",sorted=TRUE)
plot(rowData(DG.list[[i]])$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,50000))
}
hdg.DG.list<-list(
  "hdg.DG.sham"=rownames(counts(DG.list[[1]]))[1:5000],
  "hdg.DG.seizure"=rownames(counts(DG.list[[2]]))[1:5000]
)

#Run nullResiduals to calculate Pearson residuals from poisson model
for(i in 1:length(DG.list)){
DG.list[[i]]<-nullResiduals(DG.list[[i]], assay = "counts",
                      fam = "poisson", type = "pearson")
}
##Dimensionality reduction for DG
for(i in 1:length(DG.list)){
set.seed(1000) 
DG.list[[i]] <- runPCA(DG.list[[i]],exprs_values="poisson_pearson_residuals", subset_row=hdg.DG.list[[i]]) 
reducedDimNames(DG.list[[i]])
}

for(i in 1:length(DG.list)){
  set.seed(1000)
  clusters <- quickCluster(DG.list[[i]])
  DG.list[[i]] <- computeSumFactors(DG.list[[i]], cluster=clusters)
  DG.list[[i]] <- logNormCounts(DG.list[[i]])
  
  set.seed(100000)
  DG.list[[i]] <- runUMAP(DG.list[[i]], dimred="PCA",subset_row=hdg.DG.list[[i]],n_neighbors=20)
}

for(i in 1:length(DG.list)){
g <- buildSNNGraph(DG.list[[i]], k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(DG.list[[i]])$k_20_label <- factor(clust)
table(DG.list[[i]]$k_20_label)

g <- buildSNNGraph(DG.list[[i]], k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(DG.list[[i]])$k_10_label <- factor(clust)
table(DG.list[[i]]$k_10_label)
}

for(i in 1:length(DG.list)){
  DG.list[[i]]$label<-DG.list[[i]]$k_20_label
}

DG.markers.list<-list(
  "general"=c('Prox1','Glis3','Dock10','Maml2'),
  "dorsal"=c('Pdzd2','Lct','Gsg1l','Spata13'),
  "ventral"=c('Tox3','Grp','Trhr','Nr2f2')
)

pdf("DG_seizure_markers.pdf", height=6,width=12)
for(i in 1:length(DG.markers.list)){
  print(
    plotExpression(DG.list[[2]],features=c(DG.markers.list[[i]]),
                   x="label", colour_by="label", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(DG.markers.list)[i], " markers"))
  )}
dev.off()