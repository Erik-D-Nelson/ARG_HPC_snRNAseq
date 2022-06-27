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
tag.2985 <- read10xCounts(path.2985hpc.neun, col.names=T)

path.2986hpc.neun <- file.path("/dcl01/ajaffe/data/lab/singleCell/mouse_10x/premRNA/tag2986_hpc_premRNA/outs/raw_feature_bc_matrix")
tag.2986 <- read10xCounts(path.2986hpc.neun, col.names=T)

path.2987hpc.neun <- file.path("/dcl01/ajaffe/data/lab/singleCell/mouse_10x/premRNA/tag2987_hpc_premRNA/outs/raw_feature_bc_matrix")
tag.2987 <- read10xCounts(path.2987hpc.neun, col.names=T)

path.2988hpc.neun <- file.path("/dcl01/ajaffe/data/lab/singleCell/mouse_10x/premRNA/tag2988_hpc_premRNA/outs/raw_feature_bc_matrix")
tag.2988 <- read10xCounts(path.2988hpc.neun, col.names=T)


##remove unneeded files, save memory 
rm(list=ls(pattern="path"))

#make sce list
sceList.mouse10x <- list(tag.2985, tag.2986, 
                         tag.2987, tag.2988)
object.size(sceList.mouse10x)

## remove stuff to save memory
rm(tag.2985, tag.2986, tag.2987, tag.2988)

## set names
names(sceList.mouse10x) <- c("2985", "2986", 
                             "2987", "2988")

## Get unique gene symbols
# From Ensembl ID > Gene Symbol
for(i in 1:length(sceList.mouse10x)){
  rownames(sceList.mouse10x[[i]]) <- uniquifyFeatureNames(rowData(sceList.mouse10x[[i]])$ID,
                                                          rowData(sceList.mouse10x[[i]])$Symbol)
}

##Change Sample metadata column so it's not a long file name anymore
for(i in 1:length(sceList.mouse10x)){
  colData(sceList.mouse10x[[i]])$Sample<-names(sceList.mouse10x[i])
}

##Bind sce.total
sce.total<-cbind(sceList.mouse10x[[1]],sceList.mouse10x[[2]],
                 sceList.mouse10x[[3]],sceList.mouse10x[[4]])

##Add condition metadata column to colData
colData(sce.total)$condition<-factor(
  ifelse(
    colData(sce.total)$Sample %in% c("2985","2986"),
    "control","seizure")
)
# In case mess anything up
#save(sce.total, file="mouse_hpc_ECS_preprocessing.rda")

### Quality control ============================================================
##Start with barcode ranks plot
#First, use barcodeRanks() to compute ranks for each barcode (each putative nucleus)

br.out <- barcodeRanks(counts(sceList.mouse10x))


pdf('plots/barcodeRanksPlot.pdf')
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
       legend=c("knee", "inflection"))
abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
       legend=c("knee", "inflection"))
dev.off()

##run emptyDrops
#Doing this by sample as per recommendation by MNT
e.out<-list()

for(i in 1:length(sceList.mouse10x)){
set.seed(73812)
e.out[[i]] <- emptyDrops(counts(sceList.mouse10x[[i]]),niters=15000)

str(e.out[[i]], max.level=1)
print(table(Signif = e.out[[i]]$FDR <= 0.001, Limited = e.out[[i]]$Limited))
}

#Limited
#Signif  FALSE  TRUE
#FALSE 62972     0
#TRUE      1  4558
#Formal class 'DFrame' [package "S4Vectors"] with 6 slots
#Limited
#Signif  FALSE TRUE
#FALSE  8587    0
#TRUE     16 4938
#Formal class 'DFrame' [package "S4Vectors"] with 6 slots
#Limited
#Signif  FALSE  TRUE
#FALSE 61244     0
#TRUE      0  5022
#Formal class 'DFrame' [package "S4Vectors"] with 6 slots
#Limited
#Signif  FALSE TRUE
#FALSE    95    0
#TRUE     18 4949        all samples look good


# Subset:
for(i in 1:length(sceList.mouse10x)){
sceList.mouse10x[[i]] <- sceList.mouse10x[[i]][ ,which(e.out[[i]]$FDR <= 0.001)]
}
# Check
dim(sceList.10x)

##Bind sce.total
sce.total<-cbind(sceList.mouse10x[[1]],sceList.mouse10x[[2]],
                 sceList.mouse10x[[3]],sceList.mouse10x[[4]])

##Add condition metadata column to colData
colData(sce.total)$condition<-factor(
  ifelse(
    colData(sce.total)$Sample %in% c("2985","2986"),
    "control","seizure")
)

## Mito rate QC
#table(rownames(sce.total[[4]])==rownames(sce.total[[1]]))  # and checked various other pairs

location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(sce.total)$ID, 
                   column="SEQNAME", keytype="GENEID")
is.mito <- which(location=="MT")


# ID & remove those with high mito alignment
stats <- perCellQCMetrics(sce.total, subsets=list(Mito=is.mito))
colnames(stats)

# apply MAD threshold detection
high.mito <- isOutlier(stats$subsets_Mito_percent, 
                       nmads=3, type="higher",batch = sce.total$Sample)
high.mito.table <- table(high.mito)

## look at %mito
high.mito.table["TRUE"]/sum(high.mito.table)  
attributes(high.mito)


plotColData(sce.total, x = "Sample", y="subsets_Mito_percent",
            colour_by="high.mito.sample") + geom_hline(yintercept = attr(high.mito, "thresholds")["higher"])
# Save original for comparison/plotting (optional)
#sce.total.unfiltered <- sce.total   #[7] "subsets_Mito_sum"      "subsets_Mito_detected" "subsets_Mito_percent"   * <- these three added by that 'subsets' arg
#[10] "total"
## Bind stats to each SCE
colData(sce.total) <- cbind(colData(sce.total), stats)


mitoCutoffs <- attributes(high.mito)$thresholds["higher"]
mitoCutoffs <- round(mitoCutoffs, 2)



# Save
save(sce.total, e.out, high.mito, file="mouse_hpc_ECS_preprocessing.rda")

###Looking at high.mito hpc plots, some have quite a few cells with low features detected
###probably need to add both length qc step and something related to library size
###Let's try removing cells with low expressed features first?
qc.lib <- isOutlier(sce.total$sum, log=TRUE, type="lower",batch = sce.total$Sample)
qc.nexprs <- isOutlier(sce.total$detected, log=TRUE, type="lower",batch = sce.total$Sample)

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
