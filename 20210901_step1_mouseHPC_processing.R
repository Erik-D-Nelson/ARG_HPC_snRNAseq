################################################################################
### Martinowich mouse 10x snRNA-seq samples pre-processing
### ** OSCA methods - using reference script:
###   /dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot/10x-pilot_all-Frankensteins_n13_step01_MNT28Oct2019.R
###   PROJECT DIR: /dcl01/ajaffe/data/lab/singleCell/mouse_10x/analysis_MNT/
###   STEP 1: Code for preprocessing data, clustering, and Fig 1 plots
### Initiated: MNT 06Jan2020   ==#==   modified: EDN 17Nov2020
### Modification notes: Just HPC samples
################################################################################
library(SingleCellExperiment)
library(scRNAseq)
library(EnsDb.Mmusculus.v79)
library(scater)
library(scran)
library(scry)
library(DropletUtils)
library(jaffelab)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(bluster)
library(pheatmap)


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

##Change Sample metadata column so it's not a long file name anymore
for(i in 1:length(sceList.mouse10x)){
  colData(sceList.mouse10x[[i]])$sample_name<-names(sceList.mouse10x[i])
}

##Bind sce.total
sce.total<-cbind(sceList.mouse10x[[1]],sceList.mouse10x[[2]],
                 sceList.mouse10x[[3]],sceList.mouse10x[[4]])

##Add condition metadata column to colData
colData(sce.total)$condition<-factor(
  ifelse(
    colData(sce.total)$sample_name %in% c("tag2985hpc.neun","tag2986hpc.neun"),
         "control","seizure")
)
# In case mess anything up
#save(sce.total, file="mouse_hpc_ECS_preprocessing.rda")

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
#save(sce.total, e.out, file="mouse_hpc_ECS_preprocessing.rda")

# If picking up later
#load("mouse_hpc_ECS_preprocessing.rda")

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

# Save original for comparison/plotting (optional)
#sce.total.unfiltered <- sce.total   #[7] "subsets_Mito_sum"      "subsets_Mito_detected" "subsets_Mito_percent"   * <- these three added by that 'subsets' arg
#[10] "total"
## Bind stats to each SCE
colData(sce.total) <- cbind(colData(sce.total), stats, high.mito)


mitoCutoffs <- attributes(high.mito)$thresholds["higher"]
mitoCutoffs <- round(mitoCutoffs, 2)



# Save
save(sce.total, e.out, high.mito, file="mouse_hpc_ECS_preprocessing.rda")

###Looking at high.mito hpc plots, some have quite a few cells with low features detected
###probably need to add both length qc step and something related to library size
###Let's try removing cells with low expressed features first?
qc.lib <- isOutlier(sce.total$sum, log=TRUE, type="lower")
qc.nexprs <- isOutlier(sce.total$detected, log=TRUE, type="lower")

discard <- qc.lib | qc.nexprs | high.mito


##Bind library size/expressed feature logical to sce
colData(sce.total) <- cbind(colData(sce.total), qc.lib, qc.nexprs, discard)


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

save(sce.total, e.out, high.mito, discard, file="mouse_hpc_ECS_preprocessing.rda")

##discard low expressed features/small libraries
sce.total<-sce.total[,!sce.total$discard]
dim(sce.total)

### Normalization ============================================================
set.seed(1000)
clusters <- quickCluster(sce.total)
sce.total <- computeSumFactors(sce.total, cluster=clusters)
sce.total <- logNormCounts(sce.total)

### NEW Feature Selection ============================================================
#This is based on fitting Poisson model and then calculating Pearson residuals
#Let's plot highly deviant genes first; this might give us an idea about how 
#many genes to retain
sce.total<-devianceFeatureSelection(sce.total, assay="counts", 
                        fam="poisson",sorted=TRUE)

plot(rowData(sce.total)$poisson_deviance, type="l", 
     xlab="ranked genes",
     ylab="poisson deviance", 
     main="Feature Selection with Deviance",
     ylim=c(0,50000))
dev.off()

hdg<-rownames(counts(sce.total))[1:5000]
hdg<-rownames(counts(sce.total))[1:10000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce.total<-nullResiduals(sce.total, assay = "counts",
      fam = "poisson", type = "pearson")


### Dimensionality Reduction ============================================================
set.seed(1000) 
sce.total <- runPCA(sce.total,exprs_values="poisson_pearson_residuals", 
                    subset_row=hdg, BSPARAM=BiocSingular::RandomParam(),ncomponents=100) 
reducedDimNames(sce.total)

set.seed(100000)
sce.total <- runUMAP(sce.total, dimred="PCA",subset_row=hdg)

#Verify length
ncol(reducedDim(sce.total, "PCA"))

### Clustering before MNN============================================================
g50 <- buildSNNGraph(sce.total, k=50, use.dimred = 'PCA')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(sce.total)$k_50_label <- factor(clust50)
table(colData(sce.total)$k_50_label)

pdf('violinPlots_thalamusQC_Tcf7l2.pdf')
    plotExpression(sce.total, exprs_values = "logcounts", features='Tcf7l2',
                   x="k_50_label", colour_by="k_50_label", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
                   width = 0.3) 
dev.off()

##save sce.total
save(sce.total,file='sceTotal_postClustering_preThalamusRemoval.rda')

##remove thalamic clusters (express Tcf7l2)
sce.subset<-sce.total[,!sce.total$k_50_label %in% c(6,7,19)]
rm(sce.total)

##re-run analysis from devianceFeatureSelection()
sce.subset<-devianceFeatureSelection(sce.subset, assay="counts", 
                                     fam="poisson",sorted=TRUE)

pdf('rankedGenes_vs_poissonDeviance_sceSubset.pdf')
plot(rowData(sce.subset)$poisson_deviance, type="l", 
     xlab="ranked genes",
     ylab="poisson deviance", 
     main="Feature Selection with Deviance",
     ylim=c(0,50000))
dev.off()

hdg<-rownames(counts(sce.subset))[1:10000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce.subset<-nullResiduals(sce.subset, assay = "counts",
                          fam = "poisson", type = "pearson")

save(sce.subset,hdg,file="/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50.rda")
###Move on to sceSubset_getClusteredPCs.R