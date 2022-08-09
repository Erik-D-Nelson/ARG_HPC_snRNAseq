################################################################################
### Martinowich mouse 10x snRNA-seq samples pre-processing
### ** OSCA methods - using reference script:
###    /dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot/10x-pilot_all-Frankensteins_n13_step01_MNT28Oct2019.R
### STEP 01: QC datasets, separately, and integrate
###    PROJECT DIR: /dcl01/ajaffe/data/lab/singleCell/mouse_10x/analysis_MNT/
### Initiated: MNT 06Jan2020   ==#==   modified: EDN 17Nov2020
### Modification notes: Just HPC samples
################################################################################

library(SingleCellExperiment)
library(scRNAseq)
library(batchelor)
library(EnsDb.Mmusculus.v79)
library(scater)
library(scran)
library(scry)
library(uwot)
library(DropletUtils)
library(jaffelab)
library(gridExtra)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(bluster)
library(pheatmap)
library(GSEABase)
library(AUCell)
library(SingleR)

setwd("/dcl02/lieber/ajaffe/ErikNelson/")

### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
             "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
             "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
             "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

## Read in raw UMI x barcode matrix - **use pre-mRNA-aligned reads
## 4 mouse HPC samples (all NeuN-sorts); 2 ECS (2987, 2988), 2 sham (2985, 2986)
path.2985hpc.neun <- file.path("/dcl01/ajaffe/data/lab/singleCell/mouse_10x/premRNA/tag2985_hpc_premRNA/outs/raw_feature_bc_matrix")
tag2985hpc.neun <- read10xCounts(path.2985hpc.neun, col.names=T)

path.2986hpc.neun <- file.path("/dcl01/ajaffe/data/lab/singleCell/mouse_10x/premRNA/tag2986_hpc_premRNA/outs/raw_feature_bc_matrix")
tag2986hpc.neun <- read10xCounts(path.2986hpc.neun, col.names=T)

path.2987hpc.neun <- file.path("/dcl01/ajaffe/data/lab/singleCell/mouse_10x/premRNA/tag2987_hpc_premRNA/outs/raw_feature_bc_matrix")
tag2987hpc.neun <- read10xCounts(path.2987hpc.neun, col.names=T)

path.2988hpc.neun <- file.path("/dcl01/ajaffe/data/lab/singleCell/mouse_10x/premRNA/tag2988_hpc_premRNA/outs/raw_feature_bc_matrix")
tag2988hpc.neun <- read10xCounts(path.2988hpc.neun, col.names=T)

####Create new integrated SCE object#############################################################

##remove unneeded files, save memory 
rm(list=ls(pattern="path"))

## Gene annotation (from scater)
sceList.mouse10x <- list(tag2985hpc.neun, tag2986hpc.neun, 
                         tag2987hpc.neun, tag2988hpc.neun)
object.size(sceList.mouse10x)

## remove stuff to save memory
rm(tag2985hpc.neun, tag2986hpc.neun, tag2987hpc.neun, tag2988hpc.neun)

## set names
names(sceList.mouse10x) <- c("tag2985hpc.neun", "tag2986hpc.neun", 
                             "tag2987hpc.neun", "tag2988hpc.neun")

## Get unique gene symbols
# From Ensembl ID > Gene Symbol
for(i in 1:length(sceList.mouse10x)){
  rownames(sceList.mouse10x[[i]]) <- uniquifyFeatureNames(rowData(sceList.mouse10x[[i]])$ID,
                                                          rowData(sceList.mouse10x[[i]])$Symbol)
}
##Add feature annotation
for(i in 1:length(sceList.mouse10x)){
  sceList.mouse10x[[i]]<-getBMFeatureAnnos(sceList.mouse10x[[i]], 
                                           ids=rowData(sceList.mouse10x)[[i]]$ID)}


##Add length rowData column for downstream analysis
for(i in 1:length(sceList.mouse10x)){
  rowData(sceList.mouse10x)[[i]]$length<-
    rowData(sce.total)$end_position-rowData(sce.total)$start_position
}

##Add sample_name colData columns
for(i in 1:length(sceList.mouse10x)){
  colData(sceList.mouse10x[[i]])$sample_name<-names(sceList.mouse10x[i])
}

##Bind sce.total
sce.total<-cbind(sceList.mouse10x[[1]],sceList.mouse10x[[2]],sceList.mouse10x[[3]],sceList.mouse10x[[4]])

##Add mean_counts column for downstream analysis
rowData(sce.total)$mean_counts<-rowMeans(counts(sce.total))

##Add condition metadata column to colData
colData(sce.total)$condition<-
  ifelse(colData(sce.total)$sample_name %in% c("tag2985hpc.neun","tag2986hpc.neun"),
         "control","seizure")
##Make sample_name a factor
sce.total$sample_name<-factor(sce.total$sample_name)
# In case mess anything up
#save(sce.total, file="mouse10x_Dec2019_n7-processed_prelimAnalysis_EDNNov2020.rda")

### Quality control ============================================================
## - Going to ignore the adaptiveNMAD-approach to outlier detection for UMI/feature count
#    because this hasn't been as straightforward in past experience 
#       - Tends to toss the upper neuronal cells even with lower 'nmads' settings
## - Vignette for the 10x PBMC dataset (OSCA Ch.24) only does mito & droplet QC anyhow
## ** This approach used for human data processing


## Cell detection (empty droplet exclusion, rather)
# Can use UMI count vs barcode rank (knee/inflection plot) to decide threshold, but
#      "this unnecessarily discards libraries derived from cell types with low RNA content" (OSCA, Ch. 6)
#      -> Instead should prefer this Monte Carlo-simulation-based empty droplet test:
# Additionally:
# For any Sig==FALSE & Limited==TRUE, may need to increase n iterations (default = 10000) with 'niters='
#   - this field = whether "the computed p-value for a...barcode is bounded by the number of iterations"
  set.seed(109)
  e.out <- emptyDrops(counts(sce.total), niters=15000)##total is good too
  # Check if powered enough; if not, double n iterations
  ## will have to just do hands-on for those that are insufficiently powered...
  ## OR can just increase the 'niter', preemptively
  
  #e.out <- ifelse(sum(e.out$FDR > 0.001 & e.out$Limited) > 0,
                      # emptyDrops(counts(sce.total), niters=20000), e.out)


# Push save bc doing interactively...
save(sce.total, e.out, file="mouse10x_Dec2019_n7-processed_prelimAnalysis_EDNDec2020_total.rda")

# If picking up later
load("/dcl02/lieber/ajaffe/ErikNelson/mouse10x_Dec2019_n7-processed_prelimAnalysis_EDNDec2020_total.rda")

str(e.out, max.level=1)
print(table(Signif = e.out$FDR <= 0.001, Limited = e.out$Limited))


##For HPC samples, Matt got:
  #[1] "tag2985hpc.neun"
  #           Limited
  #Signif  FALSE  TRUE
  #FALSE 62968     0
  #TRUE      1  4562

  #[1] "tag2986hpc.neun"
  #           Limited
  #Signif  FALSE TRUE
  #FALSE  8587    0
  #TRUE     11 4943

  #[1] "tag2987hpc.neun"
  #           Limited
  #Signif  FALSE  TRUE
  #FALSE 61122     0
  #TRUE      0  5144

  #[1] "tag2988hpc.neun"
  #           Limited
  #Signif  FALSE TRUE
  #FALSE    95    0
  #TRUE     11 4956         - all samples good MNT 06Jan2020
##I got:
  #[1] "tag2985hpc.neun"
#Limited
#Signif  FALSE  TRUE
#FALSE 62968     0
#TRUE      1  4562
#[1] "tag2986hpc.neun"
#Limited
#Signif  FALSE TRUE
#FALSE  8587    0
#TRUE     11 4943
#[1] "tag2987hpc.neun"
#Limited
#Signif  FALSE  TRUE
#FALSE 61122     0
#TRUE      0  5144
#[1] "tag2988hpc.neun"
#Limited
#Signif  FALSE TRUE
#FALSE    95    0
#TRUE     11 4956

# Subset in for-loop:
sce.total <- sce.total[ ,which(e.out$FDR <= 0.001)]

# Check
dim(sce.total)

## Mito rate QC
#table(rownames(sce.total[[4]])==rownames(sce.total[[1]]))  # and checked various other pairs

location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(sce.total)$ID, 
                   column="SEQNAME", keytype="GENEID")
is.mito <- which(location=="MT")


# ID & remove those with high mito alignment
stats <- perCellQCMetrics(sce.total, subsets=list(Mito=is.mito))
colnames(stats)

# apply MAD threshold detection
high.mito <- isOutlier(stats$subsets_Mito_percent, nmads=3, type="higher")
high.mito.table <- table(high.mito)

## look at %mito
high.mito.table["TRUE"]/sum(high.mito.table)  
attributes(high.mito)

# Save original for comparison/plotting (optional)  -  N/A bc no high-mito-TRUE
#sce.total.unfiltered <- sce.total   #[7] "subsets_Mito_sum"      "subsets_Mito_detected" "subsets_Mito_percent"   * <- these three added by that 'subsets' arg
#[10] "total"
## Bind stats to each SCE
colData(sce.total) <- cbind(colData(sce.total), stats)
colData(sce.total) <- cbind(colData(sce.total), high.mito)


mitoCutoffs <- attributes(high.mito)$thresholds["higher"]
mitoCutoffs <- round(mitoCutoffs, 2)



## Plot metrics
pdf("plots/n7_mouse10x_totalQCmetrics_EDN19Dec2020.pdf", height=5)
  grid.arrange(
    plotColData(sce.total, y="sum", colour_by="high.mito") +
      scale_y_log10() + ggtitle("Total count"),
    plotColData(sce.total, y="detected", colour_by="high.mito") +
      scale_y_log10() + ggtitle("Detected features"),
    plotColData(sce.total, y="subsets_Mito_percent",
                colour_by="high.mito") + ggtitle(paste0("Mito % (cutoff = ", mitoCutoffs,")")),
    ncol=3
  )
dev.off()



# Save
save(sce.total, e.out, high.mito, file="mouse10x_Dec2019_n7-processed_prelimAnalysis_EDNDec2020_total.rda")

###Looking at high.mito hpc plots, some have quite a few cells with low features detected
###probably need to add both length qc step and something related to library size
###Let's try removing cells with low expressed features first?
qc.lib <- isOutlier(sce.total$sum, log=TRUE, type="lower")
qc.nexprs <- isOutlier(sce.total$detected, log=TRUE, type="lower")

discard <- qc.lib | qc.nexprs | high.mito


##Bind library size/expressed feature logical to sce
colData(sce.total) <- cbind(colData(sce.total), qc.lib)

colData(sce.total) <- cbind(colData(sce.total), qc.nexprs)

colData(sce.total) <- cbind(colData(sce.total), discard)


##Plot results of library size/expressed features
pdf("plots/n7_mouse10x_totalQCmetrics_sizeandmito_EDN19Dec2020.pdf", height=5)
  grid.arrange(
    plotColData(sce.total, y="sum", colour_by="discard") +
      scale_y_log10() + ggtitle("Total count"),
    plotColData(sce.total, y="detected", colour_by="discard") +
      scale_y_log10() + ggtitle("Detected features"),
    plotColData(sce.total, y="subsets_Mito_percent",
                colour_by="high.mito") + ggtitle(paste0("Mito % (cutoff = ", mitoCutoffs,")")),
    ncol=3
  )

dev.off()

save(sce.total, e.out, high.mito, discard, file="mouse10x_Dec2019_n7-processed_prelimAnalysis_EDNDec2020_total.rda")

##discard low expressed features/small libraries
sce.total<-sce.total[,!sce.total$discard]
dim(sce.total)



### Normalization ============================================================
set.seed(1000)
clusters <- quickCluster(sce.total)
sce.total <- computeSumFactors(sce.total, cluster=clusters)
sce.total <- logNormCounts(sce.total)


##do exploratory data analysis to see if these seizure cells have lowly expressed genese

##plot library size vs deconvolution factors
pdf("plots/n7_mouse10x_totalsizefactors_EDNDec2020.pdf", height=5)
    plot(librarySizeFactors(sce.total), 
         sizeFactors(sce.total), pch=16,
          xlab="Library size factors", ylab="Deconvolution factors", log="xy") + 
          title("Normalization factors")
dev.off()

##make a lognorm_mean_counts column in rowData
rowData(sce.total)$lognorm_mean_counts<-
    rowMeans(assay(sce.total,"logcounts"))

assay(sce.total,'lognorm_lengthnorm_counts')<-assay(sce.total,'logcounts')/rowData(sce.total)$length


###can use break and cut to make bins with equal numbers of cells
##cut provides a factor (you tell it where to cut)
##make sure to note how many cells are in each bin

##now should be able to plot counts vs gene length?
#binning for boxplots, adapted from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5428526/
q<-quantile(sqrt(rowData(sce.total)$length),probs=seq(0.1,1,0.1),na.rm=TRUE)
decile <- rep(NA,nrow(rowData(sce.total)))
decile[sqrt(rowData(sce.total)$length)<=q[1]] <- 1
for(x in 2:10) decile[sqrt(rowData(sce.total)$length)>q[x-1] & sqrt(rowData(sce.total)$length)<=q[x]] <- x
par(mar=c(8.5,4.5,3,2))
par(mfrow=c(1,3))
par(mgp=c(3,1,0))
q2 <- c(0,q^2)
labels <- rep(NA,10)
for(x in 1:10) labels[x] <-paste(round(q2[x]),"-",round(q2[x+1]),sep="")


##Make box plot
pdf("plots/boxplots_genelengthbias_EDNDec2020.pdf", height=5)
boxplot(rowData(sce.total)$lognorm_mean_counts ~ decile,ylab="LogAvgCounts",names=labels,cex.lab=1,cex.axis=1,las=2,ylim=c(0,5)) +
    title("counts by gene length",cex.main=1.5) +
    title(xlab="Gene length",cex.lab=1.5,line=7)
dev.off()


##Length is a problem, so let's divide lognorm_mean_counts by length and re boxplot
rowData(sce.total)$lognorm_lengthnorm_mean_counts<-rowData(sce.total)$lognorm_mean_counts/rowData(sce.total)$length
#rowData(sce.total)$lognorm_lengthnorm_mean_counts<-
 # rowMeans(assay(sce.total,"logcounts_lengthnorm"))

##make boxplot for comparison
pdf("plots/n7_mouse10x_boxplots_genelengthbias_EDNDec2020_lengthnorm.pdf", height=5)
boxplot(rowData(sce.total)$lognorm_lengthnorm_mean_counts ~ decile,ylab="LogAvgCounts_lengthnorm",names=labels,cex.lab=1.5,cex.axis=1.2,las=2,ylim=c(0,0.002)) +
  title("counts by gene length,length normalized",cex.main=1.5) +
  title(xlab="Gene length",cex.lab=1.5,line=7)
dev.off()

save(sce.total,file="mouse10x_Dec2019_n7-processed_prelimAnalysis_EDNDec2020_total_justsce.rda")

### INITIAL Feature Selection ============================================================
#INITIAL HVG SELECTION, moving to poisson distribution method now from 
#var.mod <- modelGeneVar(sce.total)

# Visualizing the fit:
#fit.var.mod <- metadata(var.mod)
#pdf("plots/n7_mouse10x_meanvariance_EDNDec2020.pdf", height=5)
#plot(fit.var.mod$mean, fit.var.mod$var, xlab="Mean of log-expression",
#     ylab="Variance of log-expression")
#curve(fit.var.mod$trend(x), col="dodgerblue", add=TRUE, lwd=2)
#dev.off()

#Get HVGs
#hvg <- getTopHVGs(var.mod, var.threshold=0)
#length(hvg)

#save(sce.total,hvg,file="mouse10x_Dec2019_n7-processed_prelimAnalysis_EDNDec2020_total_postnorm.rda")
### NEW Feature Selection ============================================================
#This is based on fitting Poisson model and then calculating Pearson residuals
#Let's plot highly deviant genes first; this might give us an idea about how 
#many genes to retain
sce<-devianceFeatureSelection(sce, assay="counts", 
                                    fam="poisson",sorted=TRUE)
plot(rowData(sce.total)$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,50000))
hdg<-rownames(counts(sce))[1:5000]
hdg<-rownames(counts(sce.total))[1:10000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce<-nullResiduals(sce, assay = "counts",
  fam = "poisson", type = "pearson")
#Probably going to be top 2000 genes
#Plot poisson and color by top residuals?

### Dimensionality Reduction ============================================================
set.seed(1000) 
sce.total <- runPCA(sce,exprs_values="poisson_pearson_residuals", subset_row=hdg) 
reducedDimNames(sce)

set.seed(100000)
sce.total <- runUMAP(sce, dimred="PCA",n_neighbors=20,subset_row=hdg)

#Verify length
ncol(reducedDim(sce.total, "PCA"))


### Dimensionality Reduction, minus ARGs ============================================================
arg <- read_excel("ARG_List_JesseGray.xlsx")
arg<-as.data.frame.matrix(arg)
arg <- arg[(arg$Gene_Class=="NA"),]
arg <- arg[(arg$Gene_Class=="Control"),]
arg_list<-arg$gene_name
hdg.arg.edit<-hdg[!(hdg %in% arg_list)]

sce.total.arg<-sce.total
set.seed(1000) 
sce.total.arg <- runPCA(sce.total.arg,exprs_values="poisson_pearson_residuals", subset_row=hdg.arg.edit) 
reducedDimNames(sce.total)

set.seed(100000)
sce.total.arg <- runTSNE(sce.total.arg, dimred="PCA")

#Verify length
ncol(reducedDim(sce.total.arg, "PCA"))

### Clustering before MNN============================================================
g10 <- buildSNNGraph(sce.total, k=10, use.dimred = 'PCA')
clust10 <- igraph::cluster_walktrap(g10)$membership
colData(sce.total)$k_10_label <- factor(clust10)
table(colData(sce.total)$k_10_label)

g20 <- buildSNNGraph(sce.total, k=20, use.dimred = 'PCA')
clust20 <- igraph::cluster_walktrap(g20)$membership
colData(sce.total)$k_20_label <- factor(clust20)
table(colData(sce.total)$k_20_label)

g50 <- buildSNNGraph(sce.total, k=50, use.dimred = 'PCA')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(sce.total)$k_50_label <- factor(clust50)
table(colData(sce.total)$k_50_label)

#Using a seurat-style clustering method
g.jaccard <- buildSNNGraph(sce.total, use.dimred = 'PCA', type='jaccard',k=20)
clust <- igraph::cluster_walktrap(g.jaccard)$membership
colData(sce.total)$jaccard_label_2 <- factor(clust)
table(colData(sce.total)$jaccard_label)

pdf("plots/n7_mouse10x_TSNE_colorbycluster_residuals_EDNDec2020.pdf", height=5)
plotTSNE(sce.total, colour_by="label", text_by="label")
dev.off()

pdf("plots/n7_mouse10x_TSNE_colorbysample_residuals_EDNDec2020.pdf", height=5)
plotTSNE(sce.total, colour_by="sample_name", text_by="label")
dev.off()

pdf("plots/n7_mouse10x_TSNE_colorbycondition_residuals_EDNDec2020.pdf", height=5)
plotTSNE(sce.total, colour_by="condition", text_by="label")
dev.off()
#Not sure what to make of this, but we're gonna go ahead and do marker detection
#Maybe the separation is due to IEGs, which may show up in markers?

pdf("plots/n7_mouse10x_PCA_colorbycluster_residuals_EDNDec2020.pdf", height=5)
plotReducedDim(sce.total, dimred="PCA", colour_by="label")
dev.off()

pdf("plots/n7_mouse10x_PCA_colorbysample_residuals_EDNDec2020.pdf", height=5)
plotReducedDim(sce.total, dimred="PCA", colour_by="sample_name")
dev.off()

pdf("plots/n7_mouse10x_PCA_colorbycondition_residuals_EDNDec2020.pdf", height=5)
plotReducedDim(sce.total, dimred="PCA", colour_by="condition")
dev.off()

pdf("plots/n7_mouse10x_PCAx4_colorbycluster_residuals_EDNDec2020.pdf", height=5)
plotReducedDim(sce.total, dimred="PCA",ncomponents=4, colour_by="label")
dev.off()

pdf("plots/n7_mouse10x_PCAx4_colorbysample_residuals_EDNDec2020.pdf", height=5)
plotReducedDim(sce.total, dimred="PCA",ncomponents=4, colour_by="sample_name")
dev.off()

pdf("plots/n7_mouse10x_PCAx4_colorbycondition_residuals_EDNDec2020.pdf", height=5)
plotReducedDim(sce.total, dimred="PCA",ncomponents=4, colour_by="condition")
dev.off()

save(sce.total,hdg,file="mouse10x_Dec2019_n7-processed_prelimAnalysis_EDNDec2020_total_postnorm.rda")

### Cluster diagnostics=========================================================
##Looking at transcript stats per cluster
cellIdx <- splitit(sce.total$label)
sapply(cellIdx, function(x){quantile(sce.total[ ,x]$sum)})

#Silhouette evaluation
#for k=5
sil.approx <- approxSilhouette(reducedDim(sce.total, "PCA"), clusters=clust5)

sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, clust, sil.data$other))
sil.data$cluster <- factor(clust5)

pdf("plots/n7_mouse10x_silhouette_k=5_colorbycondition_residuals_EDNDec2020.pdf", height=5)
ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
  ggbeeswarm::geom_quasirandom(method="smiley")
dev.off()

#for k=10
sil.approx <- approxSilhouette(reducedDim(sce.total, "PCA"), clusters=clust10)

sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, clust, sil.data$other))
sil.data$cluster <- factor(clust10)

pdf("plots/n7_mouse10x_silhouette_k=10_colorbycondition_residuals_EDNDec2020.pdf", height=5)
ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
  ggbeeswarm::geom_quasirandom(method="smiley")
dev.off()

#for k=20
sil.approx <- approxSilhouette(reducedDim(sce.pbmc, "PCA"), clusters=clust)

sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, clust, sil.data$other))
sil.data$cluster <- factor(clust)

pdf("plots/n7_mouse10x_silhouette_k=20_colorbycondition_residuals_EDNDec2020.pdf", height=5)
ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
  ggbeeswarm::geom_quasirandom(method="smiley")
dev.off()

#evaluate cluster separation
#k=5
ratio5 <- pairwiseModularity(g5, clust5, as.ratio=TRUE)
dim(ratio5)

library(pheatmap)
pheatmap(log2(ratio5+1), cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white", "blue"))(100))

#k=10
ratio10 <- pairwiseModularity(g10, clust10, as.ratio=TRUE)
dim(ratio10)

pheatmap(log2(ratio10+1), cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white", "blue"))(100))

#k=20
ratio <- pairwiseModularity(g20, clust, as.ratio=TRUE)
dim(ratio20)

pheatmap(log2(ratio20+1), cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white", "blue"))(100))


##checking if k=10 cluster 29 is similar to k=20 cluster 4
library(clustree)
combined <- cbind(k_10=clust10, k_5=clust5)
clustree(combined, prefix="k.", edge_arrow=FALSE)
#confirmed, it is

### Clustering before MNN,minus ARGs============================================================
g <- buildSNNGraph(sce.total.arg, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.total.arg)$clusters_ARG <- factor(clust)
table(colData(sce.total.arg)$clusters_ARG)

pdf("plots/n7_mouse10x_TSNE_ARG_colorbycluster_residuals_EDNDec2020.pdf", height=5)
plotTSNE(sce.total.arg, colour_by="clusters_ARG", text_by="clusters_ARG")
dev.off()

pdf("plots/n7_mouse10x_TSNE_ARG_colorbysample_residuals_EDNDec2020.pdf", height=5)
plotTSNE(sce.total.arg, colour_by="sample_name", text_by="clusters_ARG")
dev.off()

pdf("plots/n7_mouse10x_TSNE_ARG_colorbycondition_residuals_EDNDec2020.pdf", height=5)
plotTSNE(sce.total.arg, colour_by="condition", text_by="clusters_ARG")
dev.off()
#Not sure what to make of this, but we're gonna go ahead and do marker detection
#Maybe the separation is due to IEGs, which may show up in markers?

pdf("plots/n7_mouse10x_PCA_ARG_colorbycluster_residuals_EDNDec2020.pdf", height=5)
plotReducedDim(sce.total.arg, dimred="PCA", colour_by="clusters_ARG")
dev.off()

pdf("plots/n7_mouse10x_PCA_ARG_colorbysample_residuals_EDNDec2020.pdf", height=5)
plotReducedDim(sce.total.arg, dimred="PCA", colour_by="clusters_ARG")
dev.off()

pdf("plots/n7_mouse10x_PCA_ARG_colorbycondition_residuals_EDNDec2020.pdf", height=5)
plotReducedDim(sce.total.arg, dimred="PCA", colour_by="condition")
dev.off()

pdf("plots/n7_mouse10x_PCAx4_ARG_colorbycluster_residuals_EDNDec2020.pdf", height=5)
plotReducedDim(sce.total.arg, dimred="PCA",ncomponents=4, colour_by="clusters_ARG")
dev.off()

pdf("plots/n7_mouse10x_PCAx4_ARG_colorbysample_residuals_EDNDec2020.pdf", height=5)
plotReducedDim(sce.total.arg, dimred="PCA",ncomponents=4, colour_by="sample_name")
dev.off()

pdf("plots/n7_mouse10x_PCAx4_ARG_colorbycondition_residuals_EDNDec2020.pdf", height=5)
plotReducedDim(sce.total.arg, dimred="PCA",ncomponents=4, colour_by="condition")
dev.off()

save(sce.total.arg,hdg.arg.edit,file="mouse10x_Dec2019_n7-processed_prelimAnalysis_EDNDec2020_total_postnorm_ARG.rda")

cellIdx <- splitit(sce.total.arg$clusters_ARG)
sapply(cellIdx, function(x){quantile(sce.total.arg[ ,x]$sum)})

### Clustering after MNN============================================================
#g.mnn <- buildSNNGraph(sce.total, k=50, use.dimred = 'PCA')
#clust.mnn <- igraph::cluster_walktrap(g.mnn)$membership
#colData(sce.total)$clust.mnn <- factor(clust.mnn)
#colLabels(sce.total) <- factor(clust)
#table(colLabels(sce.total))

#pdf("plots/n7_mouse10x_TSNE_colorbycluster_EDNDec2020.pdf", height=5)
#plotTSNE(sce.total, colour_by="label")
#dev.off()

#pdf("plots/n7_mouse10x_TSNE_colorbysample_EDNDec2020.pdf", height=5)
#plotTSNE(sce.total, colour_by="sample_name")
#dev.off()

#pdf("plots/n7_mouse10x_TSNE_colorbycondition_EDNDec2020.pdf", height=5)
#plotTSNE(sce.total, colour_by="condition")
#dev.off()
#Not sure what to make of this, but we're gonna go ahead and do marker detection
#Maybe the separation is due to IEGs, which may show up in markers?

#save(sce.total,hvg,file="mouse10x_Dec2019_n7-processed_prelimAnalysis_EDNDec2020_total_postnorm.rda")



### Marker gene detection ============================================================
markers <- findMarkers(sce.total)
markers

chosen <- "6"
interesting <- markers[[chosen]]
best.set <- interesting[interesting$Top <= 4,]
logFCs <- getMarkerEffects(best.set)

library(pheatmap)
pdf("plots/n7_mouse10x_cluster6_heatmap_residuals_EDNDec2020.pdf")
pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))
dev.off()

save(sce.total,hdg,markers,file="mouse10x_Dec2019_n7-processed_prelimAnalysis_EDNDec2020_total_postnorm.rda")

### Marker gene detection, k=5==================================================
markers <- findMarkers(sce.total,groups=sce.total$k_5_label)
markers
### Marker gene detection after arg removal ============================================================
#markers.arg <- findMarkers(sce.total.arg)
#markers.arg

#save(sce.total.arg,hdg.arg.edit,markers.arg,file="mouse10x_Dec2019_n7-processed_prelimAnalysis_EDNDec2020_total_postnorm_ARG.rda")


### Cell type annotation,k=20============================================================
#First let's plot some general marker genes (e.g. GABA vs glutamate)
general_markers<-list(
  'neuron' = c('Syt1', 'Snap25', 'Grin1'),
  'excitatory_neuron' = c('Camk2a', 'Nrgn','Slc17a5', 'Slc17a6', 'Slc17a8'),
  'inhibitory_neuron' = c('Gad1', 'Gad2', 'Slc32a1'))

pdf("plots/test.pdf", height=6,width=12)
  print(
    plotExpression(x,features=c("GAD2","SLC17A6","SLC17A7","POU4F1","GPR151","TAC3"),
               x="Sex", colour_by="Sex", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3))
dev.off()


##Looks like 10,11,14,18,22,29,31 are definitely gabaergic, 4 is ambig
colData(sce.total)$level_1<-
  ifelse(colData(sce.total)$label %in% c(10,11,14,18,22,29,31),
         "GABAergic","Glutamatergic")
colData(sce.total)$level_1<-
  ifelse(colData(sce.total)$label==4,
         "Ambiguous",colData(sce.total)$level_1)
colData(sce.total)$level_1<-factor(colData(sce.total)$level_1)

##Now let's see if we can annotate the GABAergics
sce.gaba<-sce.total[,colData(sce.total)$label %in% c(4,6,10,11,14,18,22,29,31)]
pdf("plots/n7_mouse10x_interneuron_markers_residuals_EDNDec2020.pdf", height=6,width=12)
print(
  plotExpression(sce.gaba,features=c('Pvalb','Sst','Cck','Vip','Npy','Htr3a','Kit'),
               x="label", colour_by="label", point_alpha=0.5, point_size=.7,add_legend=F) +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3) + 
  ggtitle(label="interneuron markers"))
dev.off()

pdf("plots/n7_mouse10x_ARG_markers_residuals_EDNDec2020.pdf", height=6,width=12)
  print(
    plotExpression(sce.total,features=c('Egr1', 'Fos', 'Bdnf', 'Npas4','Arc','Nptx2'),
                   x="label", colour_by="label", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label="ARG Expression by Cluster")
  )
dev.off()

pdf("plots/n7_mouse10x_ARG_condition_markers_residuals_EDNJan2020.pdf", height=6,width=12)
print(
  plotExpression(sce.total,features=c('Arc', 'Fos', 'Bdnf','Ntrk2', 'Nptx2','Junb','Jun','Vgf'),
                 x="jaccard_label", colour_by="level_1", point_alpha=0.5, point_size=.7,add_legend=T)+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                 width = 0.3) +
    ggtitle(label="ARG Expression by Cluster")
)
dev.off()

pdf("plots/n7_mouse10x_ARG_condition_markers_residuals_seizures_only_EDNJan2020.pdf", height=6,width=12)
  plotExpression(sce.total$cell_type[,sce.total$condition=='seizure'],features=c('Arc', 'Fos', 'Bdnf','Ntrk2', 'Nptx2','Junb','Jun','Vgf'),
                 x="jaccard_label", colour_by="jaccard_label", point_alpha=0.5, point_size=.7,add_legend=F)+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                 width = 0.3) +
    ggtitle(label="ARG Expression by Cluster, Seizure Cells Only")
dev.off()

##2 ways to account for gene length: 1) during normalization just divide by gene length or 2) include as covariate in DE modeling
##ID HVGs
##dim reduction
##plot pcs, color by mouse, color by condition
##Point of this is to ID batch effect problems. If we see batch, need to correct usingin MNN or something

##For DE, look at differences length/no length norm
##Do DE at count level, incorporate library size, incorporate w/without gene length as covariate
### Cell type annotation,k=10===================================================
#First let's plot some general marker genes (e.g. GABA vs glutamate)
general_markers<-list(
  'excitatory_neuron' = c('Nrgn' , 'Slc17a6' , 'Slc17a7' , 'Slc17a8'),
  'inhibitory_neuron' = c('Gad1' , 'Gad2' , 'Slc32a1')
  )

pdf("plots/n7_mouse10x_general_markers_residuals_k=10_EDNDec2020.pdf", height=6,width=12)
for(i in 1:length(general_markers)){
  print(
    plotExpression(sce.total,features=c(general_markers[[i]]),
                   x="k_10_label", colour_by="k_10_label", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(general_markers)[i], " markers"))
  )}
dev.off()

##Looks like 5,7,10,15,maybe 29,31,32,33,36,42
colData(sce.total)$level_1<-
  ifelse(colData(sce.total)$k_10_label %in% c(5,7,10,15,29,31,32,33,36,42),
         "GABAergic","Glutamatergic")
colData(sce.total)$level_1<-
  ifelse(colData(sce.total)$label==29,
         "Ambiguous",colData(sce.total)$level_1)
colData(sce.total)$level_1<-factor(colData(sce.total)$level_1)


##Now lets look at the GABAergics (markers based on https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2006387#sec013)
sce.gaba<-sce.total[,sce.total$k_10_label %in% c(5,7,10,15,29,31,32,33,36,42)]
GABA_markers<-list(
  'PV basket' = c('Kit' , 'Tac1' , 'Pvalb' , 
                    'Kcnc1', 'Nr4a2', 'Satb1' , 
                        'Lhx6' , 'Erbb4' , 'Syt2'),
  'PV bistratified' = c('Kit' , 'Tac1' , 'Pvalb',
                          'Sst' , 'Npy' , 'Akr1c18',
                            'Satb1', 'Timp3' , 'Thsd7a'),
  'PV axo-axonal' = c('Kit' , 'C1ql1' , 'Pvalb' ,
                        'Tacr1' , 'Snca', 'Cck', 
                          'Pthlh' , 'Pcdh10' , 'Cpne5'),
  'Sst/Nos1' = c('Sst' , 'Nos1' , 'Npy',
                    'Pcp4' , 'Penk' , 'Chodl',
                      'Grm1', 'Chrm2' , 'Ntn1'),
  'Sst/Npy hippocamposeptal' = c('Sst' , 'Npy' , 'Calb1'),
  'Sst/Npy/Reln' = c('Sst' , 'Npy' , 'Reln' , 
                       'Calb1' , 'Grm1'),
  'Sst O-LM' = c('Sst' , 'Calb1' , 'Reln',
                   'Elfn1' , 'Serpine2', 'Npas1', 
                    'Npas3' , 'Myh8'),
  'Neurogliaform/Ivy' = c('Cacna2d1' , 'Lhx6' , 'Reln',
                            'Npy' , 'Lamp5' , 'Ndnf'),
  'Long-range projection' = c('Ntng1' , 'Rgs4'),
  'CCK cells' = c('Cck' , 'Cnr1' , 'Chrm1',
                    'Chrm3', 'Vip', 'Calb1',
                      'Slc17a8' , 'Cxcl14'),
  'Interneuron-selective' = c('Penk' , 'Calb2' , 'Vip',
                                'Cntnap5a' , 'Crh')
)

pdf("plots/n7_mouse10x_GABAergic_markers_residuals_k=10_EDNDec2020.pdf", height=6,width=12)
for(i in 1:length(interneuron_markers)){
  print(
    plotExpression(sce.gaba,features=c(interneuron_markers[[i]]),
                   x="k_10_label", colour_by="k_10_label", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(interneuron_markers)[i], " markers"))
  )}
dev.off()
##Well this is fun...
#36 is likely trilaminar long-range projection neurons (Ntng1/Chrm2/Rgs4)
#7 is likely a population of CCK neurons (Cck/Cnr1/Chrm2/Vip/Slc17a8)
#31 is likely Nos1+ neurons(including Sst/Nos1 and Ngf/Nos1) plus some Cck cells
#15 is likely PV basket cells
#I THINK 42 is PV/SST bistratified cells
#I THINK 10 is PV AAC cells + some MGE NGF/Ivy cells
#5 is some kind of SST cells (or all?)
#29 is interneuron selective interneurons
#Not sure about 32
#Not sure about 33
#still need 32, 33
save(sce.total,hdg,sce.gaba,file="mouse10x_Dec2019_n7-processed_prelimAnalysis_EDNDec2020_total_postnorm.rda")
### Cell type annotation,k=5====================================================
general_markers<-list(
  'excitatory_neuron' = c('Nrgn' , 'Slc17a6' , 'Slc17a7' , 'Slc17a8'),
  'inhibitory_neuron' = c('Gad1' , 'Gad2' , 'Slc32a1')
)

pdf("plots/n7_mouse10x_general_markers_residuals_k=5_EDNDec2020.pdf", height=6,width=20)
for(i in 1:length(general_markers)){
  print(
    plotExpression(sce.total,features=c(general_markers[[i]]),
                   x="k_5_label", colour_by="k_5_label", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(general_markers)[i], " markers"))
  )}
dev.off()

##Looks like for k=5 the definite GABAergics are: 2,15,16,24,28,32,40,41,42,48,51,54,60,61
##14 is that weird ambiguous one
##Maybe it's serotonergic or cholinergic?
plotExpression(sce.total,features=c('Chat','Slc18a3','Ache',
                                      'Tph1','Slc6a4','Fev'),
               x="k_5_label", colour_by="k_5_label", point_alpha=0.5, point_size=.7,add_legend=F)
#LOL they're not, but 27 definitely expresses Chat! That's fun!
#We'll keep working on 
colData(sce.total)$level_1<-
  ifelse(sce.total$k_5_label %in% c(2,14,15,16,24,28,32,40,41,42,48,51,54,60,61),
         "GABAergic","Glutamatergic")
colData(sce.total)$level_1_k_5<-factor(colData(sce.total)$level_1)


##Now lets look at the GABAergics (markers based on https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2006387#sec013)
sce.gaba<-sce.total[,colData(sce.total)$k_5_label %in% c(2,14,15,16,24,28,32,40,41,42,48,51,54,60,61)]
interneuron_markers<-list(
  'PV basket' = c('Kit' , 'Tac1' , 'Pvalb' , 
                  'Kcnc1', 'Nr4a2', 'Satb1' , 
                  'Lhx6' , 'Erbb4' , 'Syt2'),
  'PV bistratified' = c('Kit' , 'Tac1' , 'Pvalb',
                        'Sst' , 'Npy' , 'Akr1c18',
                        'Satb1', 'Timp3' , 'Thsd7a'),
  'PV axo-axonal' = c('Kit' , 'C1ql1' , 'Pvalb' ,
                      'Tacr1' , 'Snca', 'Cck', 
                      'Pthlh' , 'Pcdh10' , 'Cpne5'),
  'Sst/Nos1' = c('Sst' , 'Nos1' , 'Npy',
                 'Pcp4' , 'Penk' , 'Chodl',
                 'Grm1', 'Chrm2' , 'Ntn1'),
  'Sst/Npy hippocamposeptal' = c('Sst' , 'Npy' , 'Calb1'),
  'Sst/Npy/Reln' = c('Sst' , 'Npy' , 'Reln' , 
                     'Calb1' , 'Grm1'),
  'Sst O-LM' = c('Sst' , 'Calb1' , 'Reln',
                 'Elfn1' , 'Serpine2', 'Npas1', 
                 'Npas3' , 'Myh8'),
  'Neurogliaform/Ivy' = c('Cacna2d1' , 'Lhx6' , 'Reln',
                          'Npy' , 'Lamp5' , 'Ndnf'),
  'Long-range projection' = c('Ntng1' , 'Rgs4'),
  'CCK cells' = c('Cck' , 'Cnr1' , 'Chrm1',
                  'Chrm3', 'Vip', 'Calb1',
                  'Slc17a8' , 'Cxcl14'),
  'Interneuron-selective' = c('Penk' , 'Calb2' , 'Vip',
                              'Cntnap5a' , 'Crh')
)

pdf("plots/n7_mouse10x_GABAergic_markers_residuals_k=5_EDNDec2020.pdf", height=6,width=12)
for(i in 1:length(interneuron_markers)){
  print(
    plotExpression(sce.gaba,features=c(interneuron_markers[[i]]),
                   x="k_5_label", colour_by="k_5_label", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(interneuron_markers)[i], " markers"))
  )}
dev.off()
###This is fun
#Kit does not mark PV cells in HPC!!! Just striatum!!!
#Now that we've figured that out:
#2 is mostly PV basket cells with maybe some PV/SST- bistratified cells
#48 is PV/SST bistratified + axo-axonics (low Satb1/Erbb4 compared with cluster 2/60)
#60 is a small population (8 cells) of PV+/Satb1+/Erbb4+/Sst- cells. 
#Figure out the main markers that differentiate from 2
#42 is all Nos1 expressing INs, likely composed mostly of Lhx6+/Ndnf- NGF/Ivy cells
#Paper calls these MGE NGF/Ivy cells. Might also have some Nos1/Sst neurons.
#15 is a population of Cck neurons expressing Cxcl14,Reln,Serpine2.
#Paper calls these sr/slm Cck neurons
#32 and 40 are both Cck neurons too, expressing Dpy19l1/Sel1l3/Pbx1/Serpini1.
#Paper calls these "any layer" Cck neurons
#40 is Npas1+/Vip+/Slc17a8+, potentially higher in Cnr1. These are probably the "any layer" ones
#32 is Cck, maybe including some Cck projections?
#
#still need 14,16,24,28,41,51,54,61
###Cell type annotation with hipposeq markers===================================
#download and read in hipposeq markers
download.file('https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMTQ5OTcvZWxpZmUtMTQ5OTctc3VwcDEtdjEudHh0/elife-14997-supp1-v1.txt?_hash=9xZRz6wcgqh%2F81c12z8aqlirf8PCQwK%2FTvS9ZxaiVik%3D','hippo_atlas.txt')
geneset_hippo<-read.table('hippo_atlas.txt',header=TRUE,sep='\t')

#now pull out vector of markers for each cell type for use in aucell
#start with dentate granule and not dentate granule(first sets in hipposeq)
dentate_granule_cells<-geneset_hippo$gene_short_name[geneset_hippo$enriched == 'dg_d-dg_v']
not_granule_cells<-geneset_hippo$gene_short_name[geneset_hippo$enriched == 'ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v']

#We'll also add a GABAergic option from table_s1 from https://science.sciencemag.org/content/353/6302/925
#It's in an odd format--two data tables in same spreadsheet
#so, uploaded to cluster, pre-deleted the irrelevant table (second table) and saved new file to cluster
download.file('https://science.sciencemag.org/highwire/filestream/682212/field_highwire_adjunct_files/1/Table_S1.xls','markers_table.xls')
GABA<-readxl::read_excel("markers_table.xls")
GABA_markers<-GABA$'Gene names'[356:478]
level_1_ID<-list(dentate_granule_cells,not_granule_cells,GABA_markers)
names(level_1_ID)<-c('DG_cells','CA_cells','GABA_cells')

#found some duplicate colnames in sce.total...lets remove those for now
sce.total.ID<-sce.total[,unique(colnames(sce.total))]

#now let's set up the genesets
sets <- lapply(names(level_1_ID), function(x) {
  GeneSet(level_1_ID[[x]], setName=x)        
})
sets <- GeneSetCollection(sets)

#now use AUCell to build rankings 
rankings <- AUCell_buildRankings(counts(sce.total.ID),
                                 plotStats=FALSE, verbose=FALSE)
aucs <- AUCell_calcAUC(sets, rankings)
results <- t(assay(aucs))
head(results)

#Check out the results
new.labels <- colnames(results)[max.col(results)]
tab <- table(new.labels,sce.total.ID$jaccard_label_2)
tab

##This was ok, but I think there may be issues because the GABA list is so long
##Let's shorten it to just the top 40 markers
GABA<-GABA[,1:6]
GABA<-GABA[356:478,]
GABA<-as.data.frame(GABA)
rownames(GABA)<-GABA[,1]
GABA<-GABA[,-1]
GABA<-as.data.frame.matrix(GABA)

##get the top non-GABA value
GABA$max_non_GABA<-do.call(pmax,GABA[1:4])
GABA$best_markers<-GABA$'GABAergic interneurons'-GABA$max_non_GABA
GABA_markers<-rownames(GABA)[1:40,]

#remake the list with new GABA markers
level_1_ID<-list(dentate_granule_cells,not_granule_cells,GABA_markers)
names(level_1_ID)<-c('DG_cells','CA_cells','GABA_cells')
#now let's set up the genesets
sets <- lapply(names(level_1_ID), function(x) {
  GeneSet(level_1_ID[[x]], setName=x)        
})
sets <- GeneSetCollection(sets)

#now use AUCell to build rankings 
rankings <- AUCell_buildRankings(counts(sce.total.ID),
                                 plotStats=FALSE, verbose=FALSE)
aucs <- AUCell_calcAUC(sets, rankings)
results <- t(assay(aucs))
head(results)

#Check out the results
new.labels <- colnames(results)[max.col(results)]
tab <- table(new.labels,sce.total.ID$jaccard_label)
tab

##mostly straightforward.
#9,11 are DG cells
#3,4,5,7,8,10,12,13,14,15,16,17,18,19,21,22,23,24,25,26,28,29,30,31,33 are CA cells (mossy or pyramidal)
#1,2,6,20,27,32 are interneurons
#8 is ambiguous
#Let's plot expression of some general markers(including glia and ependymal just in case)

general_markers<-list(
  'excitatory_neuron' = c('Nrgn' , 'Slc17a6' , 'Slc17a7' , 'Slc17a8'),
  'inhibitory_neuron' = c('Nxph1' , 'Gad2' ,'Alk','Dlx1','Dlx6os1','Arx'),
  'CA cell' = c('Ablim3','Akap7','Arhgap20','Btg2','C1ql2','Cald1'),
  'DG cell' = c('Cnr1','Cpe','Dkk3','Dock9','Elavl4','Fhl2')
)

pdf("plots/n7_mouse10x_general_markers_residuals_jaccard_2_EDNDec2020.pdf", height=6,width=12)
for(i in 1:length(general_markers)){
  print(
    plotExpression(sce.total.ID,features=c(general_markers[[i]]),
                   x="jaccard_label", colour_by="level", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(general_markers)[i], " markers"))
  )}
dev.off()

pdf("plots/n7_mouse10x_general_markers_residuals_k=20_EDNDec2020.pdf", height=6,width=12)

pdf("plots/n7_mouse10x_TSNE_colorbyinitiallabels_residuals_EDNDec2020.pdf", height=5)
plotTSNE(sce.total.ID, colour_by="initial_cell_labels", text_by="jaccard_label")
dev.off()for(i in 1:length(general_markers)){
  print(
    plotExpression(sce.total.ID,features=c(general_markers[[i]]),
                   x="label", colour_by="label", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(general_markers)[i], " markers"))
  )}
dev.off()


#8 are kind of ambiguous...the rest I'll stick with the inital labels
colData(sce.total)$level_1<-
  factor(
    ifelse(colData(sce.total)$jaccard_label %in% c(3,4,5,7,10,12,13,
                                                      14,15,16,17,18,19,21,
                                                      22,23,24,25,26,28,29,
                                                      30,31,33),
      'CA_cells',
         ifelse(colData(sce.total)$jaccard_label %in% c(9,11),
         'DG_cells',
            ifelse(sce.total$jaccard_label %in% c(1,2,6,20,27,32),
                   'GABA_cells',
                      'Ambiguous'))))
colData(sce.total)$level_1<-factor(colData(sce.total)$level_1,levels=c('CA_cells','DG_cells','GABA_cells','Ambiguous'))
pdf("plots/n7_mouse10x_TSNE_colorbylevel1_residuals_EDNDec2020.pdf", height=5)
plotTSNE(sce.total.ID, colour_by="level_1", text_by="jaccard_label")
dev.off()


#On to level 2--mossy cells
mossy_cells<-geneset_hippo$gene_short_name[geneset_hippo$enriched == 'ca4']
not_mossy_cells<-geneset_hippo$gene_short_name[geneset_hippo$enriched == 'ca3_d-ca3_v-ca2-ca1_d-ca1_v']
level_2_ID<-list(mossy_cells,not_mossy_cells)
names(level_2_ID)<-c('Mossy_cells','Pyramidal_cells')
#Build sets and run AUCell
sets <- lapply(names(level_2_ID), function(x) {
  GeneSet(level_2_ID[[x]], setName=x)        
})
sets <- GeneSetCollection(sets)

#now use AUCell to build rankings 
rankings.mossy <- AUCell_buildRankings(counts(sce.total.ID)[,sce.total.ID$level_1=='CA_cells'],
                                 plotStats=FALSE, verbose=FALSE)
aucs.mossy <- AUCell_calcAUC(sets, rankings)
results.mossy <- t(assay(aucs))
head(results.mossy)

#Check out the results
mossy.labels <- colnames(results.mossy)[max.col(results.mossy)]
mossy.tab <- table(labels2,sce.total$jaccard_label[level_2=='CA_cells'])
mossy.tab


download.file('https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMTQ5OTcvZWxpZmUtMTQ5OTctc3VwcDItdjEudHh0/elife-14997-supp2-v1.txt?_hash=aS7i%2B0blDQBkfwCQeJC3jgO9%2FHlyznwhIdF%2FqGiwVHc%3D' , 'mossy_atlas.txt')
geneset_mossy<-read.table('mossy_atlas.txt',header=TRUE,sep='\t')
### Cell type annotation after arg removal ============================================================
#First let's plot some general marker genes (e.g. GABA vs glutamate)
#pdf("plots/n7_mouse10x_general_neuronmarkers_residuals_EDNDec2020.pdf", height=6,width=12)
#plotExpression(sce.total,features=c('Camk2a', 'Nrgn','Slc17a7', 'Slc17a6','Gad1', 'Gad2','Slc32a1'),
#              x="label", colour_by="condition", point_alpha=0.5, point_size=.7) +
#  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
#               width = 0.3, length(general.markers[[i]])) + 
#  ggtitle(label=paste0(names(general.markers)[i], " markers"))
#dev.off()

#pdf("plots/n7_mouse10x_markers_justArc_EDNDec2020.pdf", height=5)
#plotExpression(sce.total, features="Arc", x="label", colour_by="label")
#dev.off()

###NOTES=====================================================
##2 ways to account for gene length: 1) during normalization just divide by gene length or 2) include as covariate in DE modeling
##ID HVGs
##dim reduction
##plot pcs, color by mouse, color by condition
##Point of this is to ID batch effect problems. If we see batch, need to correct usingin MNN or something

##For DE, look at differences length/no length norm
##Do DE at count level, incorporate library size, incorporate w/without gene length as covariate/
rownames(interneurons) <- uniquifyFeatureNames(rowData(interneurons)$GENEID,
                                                        rowData(interneurons)$SYMBOL)
nams<-rownames(interneurons)
rownames(interneurons) = make.names(nams, unique=TRUE)

###Glutamate marker annotation==================================================
general_markers<-list(
  "granule" = c('Prox1','Pdzd2','Tox3'),
  "mossy" = c('Calb2'),
  "pyramidal" = c('Ociad2'),
  "CA2 pyramid" = c('Cacng5'),
  "CA3 pyramid" = c('Iyd','Coch'),
  "CA1 pyramid" = c('Fibcd1','Wfs1','Dcn')
)
pdf("plots/n7_mouse10x_glutamate_markers_residuals_k=20_EDNDec2020.pdf", height=6,width=20)
for(i in 1:length(general_markers)){
  print(
    plotExpression(sce.glutamate,features=c(general_markers[[i]]),
                   x="label", colour_by="label", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(general_markers)[i], " markers"))
  )}
dev.off()


###Batch correction using MNN===================================================
##First, show that there are batch effects
tab.mnn<-table(Cluster=sce.total$jaccard_label,Sample=sce.total$sample_name)
tab.mnn

#Lots of clusters with unequal #of cells from at least one sample
##Now, let's look at variation in proportion of cells per sample in each cluster
##This is from OSCA Chapter 13.6.1
# Avoid minor difficulties with the 'table' class.
tab.mnn <- unclass(tab.mnn)

# Using a large pseudo.count to avoid unnecessarily
# large variances when the counts are low.
norm <- normalizeCounts(tab.mnn, pseudo.count=10)

# Ranking clusters by the largest variances.
rv <- rowVars(norm)
DataFrame(Batch=tab.mnn, var=rv)[order(rv, decreasing=TRUE),]
#Indeed there are some major differences in variation (besides obvious between DG cell clusters)

##So, let's split the sce into seizure and control and perform MNN
sce.5558<-s3e.hb[,s3e.hb$sample_name=='Br5558']
sce.1204<-s3e.hb[,s3e.hb$sample_name=='Br1204']

#remove sce.total for now
rm(sce.total)

##First, for both, cluster and look at membership by sample
g.jaccard <- buildSNNGraph(sce.total.control, use.dimred = 'PCA', type='jaccard',k=20)
clust <- igraph::cluster_walktrap(g.jaccard)$membership
colData(sce.total.control)$jaccard_label <- factor(clust)
table(colData(sce.total.control)$jaccard_label)
sce.total.control$sample_name<-factor(sce.total.control$sample_name)
tab.mnn.control<-table(Cluster=sce.total.control$jaccard_label,Sample=sce.total.control$sample_name)
tab.mnn

set.seed(12323123)
sce.total <- runUMAP(sce.total, dimred="PCA",n_neighbors=20)

##do fastMNN() on both new SCE objects
set.seed(1000101001)
mnn.out <- fastMNN(sce.5558, sce.1204, d=50, k=20, subset.row=hdg,
          assay.type="poisson_pearson_residuals",
            BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
mnn.out

###Subclustering, DG================================================================
##First, set up sce.total dataset for k_50_label (will move up in script later)
sce.total$level_3<-as.character(sce.total$level_3)
sce.total$level_4<-factor(
  ifelse(sce.total$k_50_label %in% c(1,2,3,9,15,17,18,22,23), "CA1",
    ifelse(sce.total$k_50_label==21, "CA2",
      ifelse(sce.total$k_50_label %in% c(4,16),"CA3",
        ifelse(sce.total$k_50_label==11, "ChAT",
          ifelse(sce.total$k_50_label %in% c(6,7), "DG",
            ifelse(sce.total$k_50_label %in% c(8,10,12,14,19), "GABA",
              ifelse(sce.total$k_50_label==20, "Mossy",
                ifelse(sce.total$k_50_label==13, "Subiculum",
                  "Thalamic"))))))))
)

##ok cool, let's subcluster some stuff
sce.DG<-sce.total[,sce.total$level_4=="DG"]

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
sce.DG <- runUMAP(sce.DG, dimred="PCA",subset_row=hdg.DG)

g <- buildSNNGraph(sce.DG, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.DG)$label <- factor(clust)
table(colData(sce.DG)$label)

markers.DG <- findMarkers(sce.DG, pval.type="all", direction="up")

##Let's plot some markers
DG_sub_list<-list(
  "Hipposeq_dorsal_recs"=c("Lct","Gsg1l","Spata13","Stra6","Mgll","Cadm2"),
  "Hipposeq_ventral_recs"=c("Trhr","Cpne7","Grp","Nr2f2","Efnb2","Resp18"),
  "More_dorsal"=c("Cotl1","Brip1","Cd47","Frmpd1","Robo3","Gm11423","Akap5","Ifitm2"),
  "More_ventral"=c("Tnfaip8l3","Zfp423","Tenm4","Tox3","Ntm","Gsta4","Rnd2","Nr2f2")
)

pdf("n7_mouse10x_DG_sub_markers_residuals_EDNDec2020.pdf", height=6,width=12)
for(i in 1:length(DG_sub_list)){
  print(
    plotExpression(sce.DG,features=c(DG_sub_list[[i]]),
                   x="label", colour_by="condition", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(DG_sub_list)[i], " markers"))
  )}
dev.off()
###Subclustering, GABA==========================================================
##ok cool, let's subcluster some stuff
sce.GABA<-sce.total[,sce.total$level_4=="GABA"]

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
sce.GABA <- runUMAP(sce.GABA, dimred="PCA",subset_row=hdg.GABA)

g <- buildSNNGraph(sce.GABA, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.GABA)$label <- factor(clust)
table(colData(sce.total)$label)

markers.GABA <- findMarkers(sce.GABA, pval.type="all", direction="up")

##Let's plot some markers
GABA_sub_list<-list(
  'PV basket' = c('Kit' , 'Tac1' , 'Pvalb' , 
                  'Kcnc1', 'Nr4a2', 'Satb1' , 
                  'Lhx6' , 'Erbb4' , 'Syt2'),
  'PV bistratified' = c('Kit' , 'Tac1' , 'Pvalb',
                        'Sst' , 'Npy' , 'Akr1c18',
                        'Satb1', 'Timp3' , 'Thsd7a'),
  'PV axo-axonal' = c('Kit' , 'C1ql1' , 'Pvalb' ,
                      'Tacr1' , 'Snca', 'Cck', 
                      'Pthlh' , 'Pcdh10' , 'Cpne5'),
  'Sst/Nos1' = c('Sst' , 'Nos1' , 'Npy',
                 'Pcp4' , 'Penk' , 'Chodl',
                 'Grm1', 'Chrm2' , 'Ntn1'),
  'Sst/Npy hippocamposeptal' = c('Sst' , 'Npy' , 'Calb1'),
  'Sst/Npy/Reln' = c('Sst' , 'Npy' , 'Reln' , 
                     'Calb1' , 'Grm1'),
  'Sst O-LM' = c('Sst' , 'Calb1' , 'Reln',
                 'Elfn1' , 'Serpine2', 'Npas1', 
                 'Npas3' , 'Myh8'),
  'Neurogliaform/Ivy' = c('Cacna2d1' , 'Lhx6' , 'Reln',
                          'Npy' , 'Lamp5' , 'Ndnf'),
  'Long-range projection' = c('Ntng1' , 'Rgs4'),
  'CCK cells' = c('Cck' , 'Cnr1' , 'Chrm1',
                  'Chrm3', 'Vip', 'Calb1',
                  'Slc17a8' , 'Cxcl14'),
  'Interneuron-selective' = c('Penk' , 'Calb2' , 'Vip',
                              'Cntnap5a' , 'Crh')
)

pdf("n7_mouse10x_GABA_sub_markers_residuals_EDNDec2020.pdf", height=6,width=12)
for(i in 1:length(GABA_sub_list)){
  print(
    plotExpression(sce.GABA,features=c(GABA_sub_list[[i]]),
                   x="label", colour_by="condition", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(GABA_sub_list)[i], " markers"))
  )}
dev.off()
#who knows if k=20 is right but if so:
#1=CCK neurons
#2=PV bistratified
#3=
sce.total$level_5<-as.character(sce.total$level_4)
sce.total$level_5<-factor(
  ifelse(sce.total$k_20_label %in% c(33,23,29,27,22,11,14,17,30,26,19,9,25,16,10,3,6), "CA1",
         ifelse(sce.total$k_20_label==21, "CA2",
                ifelse(sce.total$k_20_label %in% c(2,20),"CA3",
                       ifelse(sce.total$k_20_label %in% c(5,7), "DG",
                              ifelse(sce.total$k_20_label %in% c(15,32,12,18,8), "GABA",
                                     ifelse(sce.total$k_20_label==28, "Mossy",
                                            ifelse(sce.total$k_20_label==13, "Subiculum",
                                                ifelse(sce.total$k_20_label==31,"Habenula",
                                                      ifelse(sce.total$k_20_label==34, "TRN",
                                                   "Thalamus"))))))))))

sce.subset<-sce.total[,!sce.total$label %in% c(1,4,24,31,34)]
sce.subset<-devianceFeatureSelection(sce.subset, assay="counts", 
                                     fam="poisson",sorted=TRUE)
plot(rowData(sce.subset)$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,250000))
hdg<-rownames(counts(sce.subset))[1:5000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce.subset<-nullResiduals(sce.subset, assay = "counts",
                          fam = "poisson", type = "pearson")

set.seed(1000) 
sce.subset <- runPCA(sce.subset,exprs_values="poisson_pearson_residuals", subset_row=hdg) 
reducedDimNames(sce.subset)

set.seed(100000)
sce.subset <- runUMAP(sce.subset, dimred="PCA",n_neighbors=50,subset_row=hdg,min_dist=0.2)

#Verify length
ncol(reducedDim(sce.subset, "PCA"))

g <- buildSNNGraph(sce.subset, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.subset)$k_20_label <- factor(clust)
table(colData(sce.subset)$k_20_label)

sce.subset2<-devianceFeatureSelection(sce.subset2, assay="counts", 
                                      fam="poisson",sorted=TRUE)
plot(rowData(sce.subset2)$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,50000))
hdg<-rownames(counts(sce.subset2))[1:15000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce.subset2<-nullResiduals(sce.subset2, assay = "counts",
                           fam = "poisson", type = "pearson")

set.seed(1000) 
sce.subset2 <- runPCA(sce.subset2,exprs_values="poisson_pearson_residuals", subset_row=hdg) 
reducedDimNames(sce.subset2)

set.seed(100000)
sce.subset <- runUMAP(sce.subset, dimred="PCA",n_neighbors=20,min_dist=0.1)

#Verify length
ncol(reducedDim(sce.subset2, "PCA"))

g <- buildSNNGraph(sce.subset2, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colData(sce.subset2)$label <- factor(clust)
table(colData(sce.subset2)$label)

set.seed(1000)
clusters <- quickCluster(x)
sce.subset <- computeSumFactors(x, cluster=clusters)
sce.subset <- logNormCounts(x)

sce.subset$
  
  sce.subset$annotation<-factor(
    ifelse(,
           ifelse(sce.subset$label %in% c(1,3,4,14,15,17,18), paste0("RHP.", c("1","2","3","4","5","6","7")),
                  ifelse(sce.subset$label %in% c(7,12), paste0("CA3.", c("1","2")),
                         ifelse(sce.subset$label %in% c(5,8) , paste0("DG.", c("1","2")),
                                ifelse(sce.subset$label==16,"MC",
                                       ifelse(sce.subset$label==13,"CA2",
                                              ifelse(sce.subset$label==6,"CA1",
                                                     'Sub')))))))
  )

annotationTab<-data.frame("cluster"=c(1:18),"annotation"=NA)
annotationTab$annotation[c(2,9,10)] <- paste0("GABA.", c("1","2","3"))
annotationTab$annotation[c(1,3,4,14,15,17,18)]<-paste0("RHP.", c("1","2","3","4","5","6","7"))
annotationTab$annotation[c(7,12)]<-paste0("CA3.", c("1","2"))
annotationTab$annotation[c(5,8)]<-paste0("DG.", c("1","2"))
annotationTab$annotation[c(6,11,13,16)]<-c("CA1","Sub","CA2","CA4")
annotationTab$cluster<-factor(annotationTab$cluster)
annotationTab$annotation<-factor(annotationTab$annotation)

sce.subset$annotation <- annotationTab$annotation[match(sce.subset$label,
                                                         annotationTab$cluster)]


plotGroupedHeatmap(sce.subset,features=c('Fibcd1','Ntsr2','Galnt3','Csf2rb2','Glis3','Gad2','Nts'),cluster_rows=F,cluster_cols=F,group='annotation',center=F)
plotGroupedHeatmap(sce.DG,features=c('Nptx2','Rasgrp1','Scg2','Sv2b','Kitl','Sv2c','Nr4a3','Slc6a17','Rheb','Maml3','Fam126b'),cluster_rows=F,cluster_cols=F,group='label',center=F)


plotDots(
  sce.subset,
  features=c('Fibcd1','Ntsr2','Galnt3','Csf2rb2','Glis3','Gad2','Nts'),
  group = 'annotation',
  block = NULL,
  exprs_values = "logcounts",
  detection_limit = 0,
  low_color = "white",
  high_color = "blue",
  max_ave = NULL,
  max_detected = NULL,
  swap_rownames = NULL,
  cluster_rows=F
)

###Set up data==================================================================
#Subset for DG only and aggregateAcrossCells to pseudobulk
summed <- aggregateAcrossCells(sce.subset, 
                               ids=colData(sce.subset)[c('condition','sample_name','annotation')])

#Use perCellQCMetrics to get library sizes for each sample 
#(didn't carry over after pseudobulk, may be an easier way but this worked
#so I ran with it)
location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(summed)$GENEID, 
                   column="SEQNAME", keytype="GENEID")
is.mito <- which(location=="MT")
stats <- perCellQCMetrics(summed, subsets=list(Mito=is.mito))
colnames(stats)
colData(summed)$sum<-stats$sum

#filter out small samples
summed <- summed[,summed$ncells >= 20]

#Filter out lowly expressed genes and remove genes with length=NA
summed<-summed[filterByExpr(summed, group=summed$sample_name),]

#make DGElist
y <- DGEList(counts(summed), samples=colData(summed),
             genes=rowData(summed))
y$samples$condition<-factor(y$samples$condition)

#make design matrix
design <- model.matrix(~condition+annotation, y$samples)

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

ids <- c("Nptx2" , 'Bdnf' , 'Apoe')
gene.labels <- x[x$Symbol %in% ids,]
text(x=gene.labels$logCPM, y=gene.labels$logFC,
labels=rownames(gene.labels), cex=0.7, pos=1)

pdf('20210913_overall_MA_plot.pdf',height=5,width=7)
plotSmear(elrt, de.tags = rownames(x))
text(x=gene.labels$logCPM, y=gene.labels$logFC,
     labels=rownames(gene.labels), cex=0.7, pos=1)

dev.off()

de.specific <- pseudoBulkDGE(summed.spec,
                             label=summed.spec$annotation,
                             design=~condition,
                             coef=2,
                             condition=summed.spec$condition
)

de.specific <- pseudoBulkDGE(summed,
                             label=summed$annotation,
                             design=~condition,
                             coef=2,
                             condition=summed$condition
)

is.de <- decideTestsPerLabel(de.specific, threshold=0.05)
summarizeTestsPerLabel(is.de)
