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
library(Rtsne)
library(GenomicRanges)
library(ggplot2)
library(dplyr)

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

####Create new object, cbind counts, rowbind coldata######

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
  rownames(sceList.mouse10x)[[i]] <- uniquifyFeatureNames(rowData(sceList.mouse10x)[[i]]$ID,
                                                          rowData(sceList.mouse10x)[[i]]$Symbol)
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
for(i in 1:length(sce.total)){
  colData(sce.total)$sample_name<-names(sce.total[i])
}

##Bind sce.total
sce.total<-cbind(sce.total[[1]],sce.total[[2]],sce.total[[3]],sce.total[[4]])

##Add mean_counts column for downstream analysis
rowData(sce.total)$mean_counts<-rowMeans(counts(sce.total))

##Add condition metadata column to colData
colData(sce.total)$condition<-
  ifelse(colData(sce.total)$sample_name %in% c("tag2985hpc.neun","tag2986hpc.neun"),
         "control","seizure")

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

#####################
##Matt stopped here##
#####################

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

##########Gene length vs counts########

##first, simple scatter plot
#pdf("plots/n7_mouse10x_plotRowData_genelengthbias_EDNNov2020.pdf", height=5)
#for(i in 1:length(sce.total)){
#  grid.arrange(
#    plotRowData(sce.total, y = "mean_counts",x = "length") +
#      scale_y_log10() + scale_x_log10() +
#      ggtitle(paste0("Mean counts by gene length:",names(sce.total)[[i]])),
#      ncol=1)
#  }
#dev.off()

##From these scatter plots, pretty clear there's a gene length bias.
##let's do normalization now so we can make box plots

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
pdf("plots/n7_mouse10x_boxplots_genelengthbias_EDNNov2020.pdf", height=5)
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
sce.total<-devianceFeatureSelection(sce.total, assay="counts", 
                                    fam="poisson",sorted=TRUE)
plot(rowData(sce.total)$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,50000))
hdg<-rownames(counts(sce.total))[1:5000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce.total<-nullResiduals(sce.total, assay = "counts",
  fam = "poisson", type = "pearson")
#Probably going to be top 2000 genes
#Plot poisson and color by top residuals?

### Dimensionality Reduction ============================================================
set.seed(1000) 
sce.total <- runPCA(sce.total,exprs_values="poisson_pearson_residuals", subset_row=hdg) 
reducedDimNames(sce.total)

set.seed(100000)
sce.total <- runTSNE(sce.total, dimred="PCA")

#Verify length
ncol(reducedDim(sce.total, "PCA"))


### Dimensionality Reduction, minus ARGs ============================================================
arg <- read_excel("ARG_List_JesseGray.xlsx")
arg<-as.data.frame.matrix(arg)
arg <- arg[!(arg$Gene_Class=="NA"),]
arg <- arg[!(arg$Gene_Class=="Control"),]
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
g <- buildSNNGraph(sce.total, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce.total) <- factor(clust)
table(colLabels(sce.total))

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

cellIdx <- splitit(sce.total$label)
sapply(cellIdx, function(x){quantile(sce.total[ ,x]$sum)})
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

### Marker gene detection after arg removal ============================================================
#markers.arg <- findMarkers(sce.total.arg)
#markers.arg

#save(sce.total.arg,hdg.arg.edit,markers.arg,file="mouse10x_Dec2019_n7-processed_prelimAnalysis_EDNDec2020_total_postnorm_ARG.rda")


### Cell type annotation ============================================================
#First let's plot some general marker genes (e.g. GABA vs glutamate)
markers<-list(
  'neuron' = c('Syt1', 'Snap25', 'Grin1'),
  'excitatory_neuron' = c('Camk2a', 'Nrgn','Slc17a5', 'Slc17a6', 'Slc17a8'),
  'inhibitory_neuron' = c('Gad1', 'Gad2', 'Slc32a1'))

pdf("plots/n7_mouse10x_general_markers_residuals_EDNDec2020.pdf", height=6,width=12)
for(i in 1:length(markers)){
  print(
    plotExpression(sce.total,features=c(markers[[i]]),
               x="label", colour_by="label", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(markers)[i], " markers"))
)}
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
  plotExpression(sce.gaba,features=c('Pvalb','Sst','Cck','Vip','Npy','Htr3a'),
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

pdf("plots/n7_mouse10x_ARG_condition_markers_residuals_EDNDec2020.pdf", height=6,width=12)
print(
  plotExpression(sce.total,features=c('Egr1', 'Fos', 'Bdnf', 'Npas4','Arc','Nptx2'),
                 x="label", colour_by="condition", point_alpha=0.5, point_size=.7,add_legend=F)+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                 width = 0.3) +
    ggtitle(label="ARG Expression by Cluster")
)
dev.off()

##2 ways to account for gene length: 1) during normalization just divide by gene length or 2) include as covariate in DE modeling
##ID HVGs
##dim reduction
##plot pcs, color by mouse, color by condition
##Point of this is to ID batch effect problems. If we see batch, need to correct usingin MNN or something

##For DE, look at differences length/no length norm
##Do DE at count level, incorporate library size, incorporate w/without gene length as covariate
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
##Do DE at count level, incorporate library size, incorporate w/without gene length as covariate