
library(SingleCellExperiment)
library(scater)
library(scran)
library(scry)
library(DropletUtils)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(bluster)
library(pheatmap)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(EnsDb.Mmusculus.v79)
library(scDblFinder)


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
names(sceList.mouse10x) <- c("Sham1", "Sham2", 
                             "ECS1", "ECS2")

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

# In case mess anything up
#save(sce.total, file="mouse_hpc_ECS_preprocessing.rda")

### Quality control ============================================================

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
lapply(sceList.mouse10x,dim)

##Bind sce.total
sce.total<-cbind(sceList.mouse10x[[1]],sceList.mouse10x[[2]],
                 sceList.mouse10x[[3]],sceList.mouse10x[[4]])

##Add condition metadata column to colData
colData(sce.total)$condition<-factor(
  ifelse(
    colData(sce.total)$Sample %in% c("Sham1","Sham2"),
    "Sham","ECS"),
  levels=c('Sham','ECS')
)

save(sce.total,file='sce_total_preThalamusRemoval.rda')

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
colData(sce.total) <- cbind(colData(sce.total), stats)
sce.total$high.mito<-high.mito


#plotColData(sce.total, x = "Sample", y="subsets_Mito_percent",
#            colour_by="high.mito") + geom_hline(yintercept = attr(high.mito, "thresholds")["higher"])
# Save original for comparison/plotting (optional)
#sce.total.unfiltered <- sce.total   #[7] "subsets_Mito_sum"      "subsets_Mito_detected" "subsets_Mito_percent"   * <- these three added by that 'subsets' arg
#[10] "total"
## Bind stats to each SCE



mitoCutoffs <- attributes(high.mito)$thresholds["higher"]
mitoCutoffs <- round(mitoCutoffs, 2)



# Save
#save(sce.total, e.out, high.mito, file="processed_data/mouse_hpc_ECS_preprocessing.rda")

###Looking at high.mito hpc plots, some have quite a few cells with low features detected
###probably need to add both length qc step and something related to library size
###Let's try removing cells with low expressed features first?
qc.lib <- isOutlier(sce.total$sum, log=TRUE, type="lower",batch = sce.total$Sample)
qc.nexprs <- isOutlier(sce.total$detected, log=TRUE, type="lower",batch = sce.total$Sample)

discard <- qc.lib | qc.nexprs | high.mito


##Bind library size/expressed feature logical to sce
colData(sce.total) <- cbind(colData(sce.total), qc.lib, qc.nexprs, discard)


##Plot results of library size/expressed features
pdf("plots/totalQCmetrics_sizeandmito.pdf", height=5)
grid.arrange(
  plotColData(sce.total, x='Sample', y="sum", colour_by="discard") +
    scale_y_log10() + ggtitle("Total count"),
  plotColData(sce.total, x='Sample', y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(sce.total, x='Sample', y="subsets_Mito_percent",
              colour_by="high.mito"),
  ncol=1
)

dev.off()

#save(sce.total, e.out, high.mito, discard, file="mouse_hpc_ECS_preprocessing.rda")

##discard low expressed features/small libraries
sce.total<-sce.total[,!sce.total$discard]
dim(sce.total)

##remove all genes with 0 counts
keep <- rowSums(counts(sce.total) > 0) > 0
sce.total <- sce.total[keep, ]

### Normalization ============================================================
set.seed(1000)
clusters <- quickCluster(sce.total)
sce.total <- computeSumFactors(sce.total, cluster=clusters)
sce.total <- logNormCounts(sce.total)

###Feature Selection ============================================================
#This is based on fitting Poisson model and then calculating Pearson residuals
#Let's plot highly deviant genes first; this might give us an idea about how 
#many genes to retain
sce.total<-devianceFeatureSelection(sce.total, assay="counts", 
                                    fam="poisson",sorted=TRUE)

pdf("plots/rank_vs_deviance_plot_sceTotal.pdf")
plot(rowData(sce.total)$poisson_deviance, type="l", 
     xlab="ranked genes",
     ylab="poisson deviance", 
     main="Feature Selection with Deviance, all nuclei",
     ylim=c(0,500000))
abline(v=5000,lty='dashed',col="red")
dev.off()

hdg<-rownames(counts(sce.total))[1:5000]


#Run nullResiduals to calculate Pearson residuals from poisson model
sce.total<-nullResiduals(sce.total, assay = "counts",
                         fam = "poisson", type = "pearson")


### Dimensionality Reduction ============================================================
set.seed(1000) 
sce.total <- runPCA(sce.total,exprs_values="poisson_pearson_residuals", 
                    subset_row=hdg,ncomponents=100) 
reducedDimNames(sce.total)

set.seed(100000)
sce.total <- runUMAP(sce.total, dimred="PCA",subset_row=hdg)

#Verify length
ncol(reducedDim(sce.total, "PCA"))

### Clustering============================================================
g50 <- buildSNNGraph(sce.total, k=50, use.dimred = 'PCA')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(sce.total)$label <- factor(clust50)
table(colData(sce.total)$label)

##doublet detection and removal
sample_id_names<- names(table(sce.total$Sample))
names(sample_id_names)<-sample_id_names

sample_id_rse<- map(sample_id_names,~sce.total[,sce.total$Sample==.x])
## To speed up, run on sample-level top-HVGs - just take top 1000 ===
pilot.data.normd <- map(sample_id_rse, ~logNormCounts(.x))
#geneVar.samples <- map(pilot.data.normd, ~modelGeneVar(.x))
#topHVGs <- map(geneVar.samples, ~getTopHVGs(.x, n = 2000))

# Generate doublet density scores
set.seed(109)
dbl.dens.focused <- map(names(pilot.data.normd), ~computeDoubletDensity(pilot.data.normd[[.x]], subset.row=hdg,dims=50))
names(dbl.dens.focused) <- names(pilot.data.normd)

for(i in names(sample_id_rse)){
  pilot.data.normd[[i]]$doublet.score <- dbl.dens.focused[[i]]
}
ref<-rbind(colData(pilot.data.normd[[1]]),colData(pilot.data.normd[[2]]),
           colData(pilot.data.normd[[3]]),colData(pilot.data.normd[[4]]))
identical(colnames(sce.total),rownames(ref))

for(i in names(sample_id_rse)){
  pilot.data.normd[[i]]$doubletScore <- dbl.dens.focused[[i]]
}

dbl.dens<-lapply(d,c)

summary(dbl.dens)

quantile(dbl.dens, probs=seq(0,1,by=0.01),3)

#0%         1%         2%         3%         4%         5%         6%
#0.0000000  0.0000000  0.0000000  0.0088900  0.0090360  0.0096320  0.0174280
#7%         8%         9%        10%        11%        12%        13%
#0.0177800  0.0180720  0.0192640  0.0261420  0.0266700  0.0288960  0.0288960
#14%        15%        16%        17%        18%        19%        20%
#0.0348560  0.0355600  0.0385280  0.0435700  0.0444500  0.0451800  0.0481600
#21%        22%        23%        24%        25%        26%        27%
#0.0522840  0.0533400  0.0577920  0.0609980  0.0622300  0.0632520  0.0697120
#28%        29%        30%        31%        32%        33%        34%
#0.0711200  0.0770560  0.0784260  0.0813240  0.0866880  0.0889000  0.0958540
#35%        36%        37%        38%        39%        40%        41%
#0.0963200  0.0993960  0.1059520  0.1084320  0.1155700  0.1174680  0.1244600
#42%        43%        44%        45%        46%        47%        48%
#0.1265040  0.1348480  0.1394240  0.1444800  0.1481380  0.1541120  0.1568520
#49%        50%        51%        52%        53%        54%        55%
#0.1637440  0.1716840  0.1742800  0.1829940  0.1897560  0.1955800  0.2022720
#56%        57%        58%        59%        60%        61%        62%
#0.2091360  0.2168640  0.2259000  0.2311680  0.2408000  0.2504320  0.2600640
#63%        64%        65%        66%        67%        68%        69%
#0.2696960  0.2793280  0.2889600  0.2985920  0.3111500  0.3224180  0.3343320
#70%        71%        72%        73%        74%        75%        76%
#0.3485600  0.3659880  0.3795120  0.3949120  0.4141760  0.4337280  0.4527040
#77%        78%        79%        80%        81%        82%        83%
#0.4711700  0.4889500  0.5141260  0.5393920  0.5664100  0.5956300  0.6274080
#84%        85%        86%        87%        88%        89%        90%
#0.6622640  0.7023100  0.7407162  0.7898240  0.8365440  0.8945640  0.9578160
#91%        92%        93%        94%        95%        96%        97%
#1.0223500  1.1076800  1.2017880  1.3070730  1.4579600  1.7530240  2.5094994
#98%        99%       100%
#4.5573262  9.0132210 27.0081280




save(sce.total,file='sce_total_preThalamusRemoval.rda')

##remove putative doublets
dubs<-do.call(cbind,pilot.data.normd)
dubs<-dubs[,match(colnames(sce.total),colnames(dubs))]
sce.total$doubletScore<-dubs$doubletScore
rm(dubs,pilot.data.normd)
colData(sce.total)$doublet<-ifelse(sce.total$doubletScore > 5, T, F)
sce.total<-sce.total[,sce.total$doublet==FALSE]

##find marker genes
markers<-findMarkers(sce.total,pval.type='all',
                     direction='up',group=sce.total$label)

save(markers,file='processed_data/sceTotal_markerGenes.rda')
##thalamus detection and removal
features<-c('Snap25','Tcf7l2','Zfhx3','Shox2','Ano1','Six3')
umap<-list()
for(i in 1:length(features)){
  umap[[i]]<-plotUMAP(sce.total,text_by='label',colour_by=features[[i]])
}


pdf('plots/figS1_umaps_thalamusQC.pdf',h=9,w=8)
grid.arrange(grobs=umap, ncol=2)
dev.off()

##save sce.total
save(sce.total,file='sce_total_preThalamusRemoval.rda')

##remove thalamic clusters (express Tcf7l2)
sce.subset<-sce.total[,!sce.total$label %in% c(5,9,14)]
rm(sce.total)

##re-run analysis from devianceFeatureSelection()
sce.subset<-devianceFeatureSelection(sce.subset, assay="counts", 
                                     fam="poisson",sorted=TRUE)

pdf('plots/rank_vs_poissonDeviance_sceSubset.pdf')
plot(rowData(sce.subset)$poisson_deviance, type="l", 
     xlab="ranked genes",
     ylab="poisson deviance", 
     main="Feature Selection with Deviance",
     ylim=c(0,500000))
abline(v=3000,col='red',lty='dashed')
dev.off()
hdg<-rownames(counts(sce.subset))[1:3000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce.subset<-nullResiduals(sce.subset, assay = "counts",
                          fam = "poisson", type = "pearson")

#PCA
set.seed(2014) 
sce.subset <- runPCA(sce.subset,exprs_values="poisson_pearson_residuals", 
                     subset_row=hdg,ncomponents=50) 
reducedDimNames(sce.subset)

set.seed(202014)
sce.subset <- runUMAP(sce.subset, dimred="PCA")

#Verify length
ncol(reducedDim(sce.subset, "PCA"))
sce.subset

pdf('umap.pdf')
plotUMAP(sce.subset,colour_by='k_25_label',text_by='k_25_label')
dev.off()

save(sce.subset,hdg,file="sce_subset.rda")



#Make sce objects for GitHub release
# assay(sce.total,'poisson_pearson_residuals')<-NULL
# assay(sce.total,'logcounts')<-NULL
# reducedDim(sce.total,'PCA')<-NULL
# reducedDim(sce.total,'UMAP')<-NULL
# rowData(sce.total)<-rowData(sce.total)[,c(1,2)]
# colData(sce.total)<-colData(sce.total)[,c(1,2,3,15)]
# #save(sce.total,compress='xz',file='processed_data/sce_total.rda.xz')
# 
# 
# assay(sce.subset,'poisson_pearson_residuals')<-NULL
# assay(sce.subset,'logcounts')<-NULL
# reducedDim(sce.subset,'PCA')<-NULL
# reducedDim(sce.subset,'UMAP')<-NULL
# rowData(sce.subset)<-rowData(sce.subset)[,c(1,2)]
# colData(sce.subset)<-colData(sce.subset)[,c(1,2,3,21,24,25)]
# save(sce.subset,compress='xz',file='processed_data/sce_subset.rda.xz')



###move on to clustering/annotation/figure 1 plots 
