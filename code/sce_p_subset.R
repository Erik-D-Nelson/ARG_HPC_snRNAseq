###Subclustering, DG==========================================================
##ok cool, let's subcluster some stuff
sce.DG<-sce.subset[,sce.subset$level_6=="DG"]

##Feature selection for DG
sce.DG<-devianceFeatureSelection(sce.DG, assay="counts", 
                                   fam="poisson",sorted=TRUE)
plot(rowData(sce.DG)$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,50000))

hdg.DG<-rownames(counts(sce.DG))[1:5000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce.DG<-nullResiduals(sce.DG, assay = "counts",
                        fam = "poisson", type = "pearson")

##Dimensionality reduction for DG
set.seed(1000) 
sce.DG <- runPCA(sce.DG,exprs_values="poisson_pearson_residuals", subset_row=hdg.DG) 
reducedDimNames(sce.DG)

set.seed(100000)
sce.DG <- runUMAP(sce.DG, dimred="PCA",subset_row=hdg.DG,n_neighbors=20)

g <- buildSNNGraph(sce.DG, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.DG)$k_20_label <- factor(clust)
table(sce.DG$k_20_label)

g <- buildSNNGraph(sce.DG, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.DG)$k_10_label <- factor(clust)
table(sce.DG$k_10_label)

g <- buildSNNGraph(sce.DG, k=50, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.DG)$k_50_label <- factor(clust)
table(sce.DG$k_50_label)

set.seed(1000)
clusters <- quickCluster(sce.DG)
sce.DG <- computeSumFactors(sce.DG, cluster=clusters)
sce.DG <- logNormCounts(sce.DG)

summed<-sce.subset[,sce.subset$level_6=='DG']
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
x3<-topTags(elrt,n=10000)
x3<-as.data.frame(x3)
#x3<-x3[rownames(x3) %in% arg_list,15:19]
x3<-x3[x3$FDR<=0.05,]

dg<-as.data.frame(dg)
dg<-dg[dg$FDR<=0.05,]
x3[,14:18]

x<-topTags(elrt,n=500)
x<-as.data.frame(x)
x<-x[x$Symbol %in% arg_list,]
x[,14:18]

###MNN correction
set.seed(100100100)
mnn.out <- fastMNN(sce.DG, batch=sce.DG$sample_name, 
                    subset.row=hdg.DG,k=50)
colData(mnn.out)<-colData(sce.DG)

set.seed(100000)
mnn.out<-runUMAP(mnn.out,dimred='corrected',exprs_values='reconstructed')
plotUMAP(mnn.out,colour_by='condition',text_by='k_20_label')

set.seed(1000)
clusters <- quickCluster(mnn.out)
mnn.out <- computeSumFactors(mnn.out, cluster=clusters)
mnn.out <- logNormCounts(sce.DG)

###Subclustering, seizure and sham separate=====================================
sce.DG.sham<-sce.DG[,sce.DG$condition=="control"]
sce.DG.seizure<-sce.DG[,sce.DG$condition=='seizure']
DG.list<-list(
  "sham"=sce.DG.sham,
  "seizure"=sce.DG.seizure
)

##Feature selection for DG
sce.DG<-devianceFeatureSelection(sce.DG, assay="counts", 
                                 fam="poisson",sorted=TRUE)
plot(rowData(sce.DG)$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,50000))

hdg.DG<-rownames(counts(sce.DG))[1:5000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce.DG<-nullResiduals(sce.DG, assay = "counts",
                      fam = "poisson", type = "pearson")

##Dimensionality reduction for DG
set.seed(1000) 
sce.DG <- runPCA(sce.DG,exprs_values="poisson_pearson_residuals", subset_row=hdg.DG) 
reducedDimNames(sce.DG)

set.seed(100000)
sce.DG <- runUMAP(sce.DG, dimred="PCA",subset_row=hdg.DG,n_neighbors=20)

g <- buildSNNGraph(sce.DG, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.DG)$k_20_label <- factor(clust)
table(sce.DG$k_20_label)

g <- buildSNNGraph(sce.DG, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.DG)$k_10_label <- factor(clust)
table(sce.DG$k_10_label)

g <- buildSNNGraph(sce.DG, k=50, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.DG)$k_50_label <- factor(clust)
table(sce.DG$k_50_label)

set.seed(1000)
clusters <- quickCluster(sce.DG)
sce.DG <- computeSumFactors(sce.DG, cluster=clusters)
sce.DG <- logNormCounts(sce.DG)