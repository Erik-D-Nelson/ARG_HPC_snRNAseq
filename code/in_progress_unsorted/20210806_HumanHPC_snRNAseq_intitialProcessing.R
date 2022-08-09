################################################################################
###Human HPC analysis, initial processing
##Sample information: 6 samples from 3 dissections of postmortem human HPC (2 per brain)
##3 NeuN sorted (SCE objects labeled to reflect)
##Objective: Build 1 comprehensive SCE object, quality control steps, reduce 
##dimensions, feature selection, MNN batch correction, graph based clustering
##Packages used================================================================
library(SingleCellExperiment)
library(scRNAseq)
library(batchelor)
library(EnsDb.Hsapiens.v86)
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
library(rtracklayer)

##Build SCE object==============================================================
path.5212hpc.neun <- file.path("/dcl02/lieber/ajaffe/ErikNelson/ARG_HPC/Br5212/outs/raw_feature_bc_matrix")
tag5212hpc.neun <- read10xCounts(path.5212hpc.neun, col.names=T)

path.5287hpc.neun <- file.path("/dcl02/lieber/ajaffe/ErikNelson/ARG_HPC/Br5287/outs/raw_feature_bc_matrix/")
tag5287hpc.neun <- read10xCounts(path.5287hpc.neun, col.names=T)

path.5161hpc.neun <- file.path("/dcl02/lieber/ajaffe/ErikNelson/ARG_HPC/Br5161/outs/raw_feature_bc_matrix/")
tag5161hpc.neun <- read10xCounts(path.5161hpc.neun, col.names=T)

path.5212hpc <- file.path("/dcl02/lieber/ajaffe/ErikNelson/ARG_HPC/Br5212_MNT/outs/raw_feature_bc_matrix")
tag5212hpc <- read10xCounts(path.5212hpc, col.names=T)

path.5287hpc <- file.path("/dcl02/lieber/ajaffe/ErikNelson/ARG_HPC/Br5287_MNT/outs/raw_feature_bc_matrix/")
tag5287hpc <- read10xCounts(path.5287hpc, col.names=T)

path.5161hpc <- file.path("/dcl02/lieber/ajaffe/ErikNelson/ARG_HPC/Br5161_MNT/outs/raw_feature_bc_matrix/")
tag5161hpc <- read10xCounts(path.5161hpc, col.names=T)

rm(list=ls(pattern="path."))

##make SCE list
sceList <- list(tag5212hpc.neun, tag5287hpc.neun, tag5161hpc.neun,
                tag5212hpc, tag5287hpc, tag5161hpc)
object.size(sceList)

## remove stuff to save memory
rm(tag5212hpc.neun, tag5287hpc.neun, tag5161hpc.neun,
   tag5212hpc, tag5287hpc, tag5161hpc)

## set names
names(sceList) <- c('tag5212hpc.neun', 'tag5287hpc.neun', 'tag5161hpc.neun',
                    'tag5212hpc', 'tag5287hpc', 'tag5161hpc')

##Make new Sample column
for(i in 1:length(sceList)){
  colData(sceList[[i]])$Sample<-names(sceList[i])
}

##Bind sce.human.hpc
sce.human.hpc<-do.call("cbind",sceList)

rm(sceList)

#Unique feature names
rownames(sce.human.hpc) <- uniquifyFeatureNames(rowData(sce.human.hpc)$ID, rowData(sce.human.hpc)$Symbol)

#Make Sample a factor
sce.human.hpc$Sample<-factor(sce.human.hpc$Sample)

### Quality control ============================================================
##Empty droplet removal using Monte Carlo-based simulation
cat(paste0("Simulating empty drops for: sce.human.hpc... \n"))
set.seed(109)
e.out <- emptyDrops(counts(sce.human.hpc),niters=30000)

print(table(Signif = e.out$FDR <= 0.001, Limited = e.out$Limited))
#Signif   FALSE   TRUE
#FALSE 105976      0
#TRUE    1285  21037

sce.human.hpc <- sce.human.hpc[ ,which(e.out$FDR <= 0.001)]

save(sce.human.hpc, e.out,
     file="sce_humHPC_processing.rda")

### Mito rate QC ===
# MAD approach for mito rate thresholding
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce.human.hpc)$ID, 
                   column="SEQNAME", keytype="GENEID")
## Warning message: Unable to map 2611 of 36601 requested IDs.

## ID  mito genes
stats <- perCellQCMetrics(sce.human.hpc, subsets=list(Mito=which(location=="MT")))

##Find outliers
high.mito <- isOutlier(stats$subsets_Mito_percent, nmads=3, type="higher")
table(high.mito)
#high.mito
#FALSE  TRUE 
#15904  6418 
attributes(high.mito) 
#$thresholds
#lower    higher 
#-Inf 0.9618738
mitoCutoffs <- attributes(high.mito)$thresholds["higher"]
mitoCutoffs <- round(mitoCutoffs, 2)


# Bind stats to colData
colData(sce.human.hpc) <- cbind(colData(sce.human.hpc), stats, high.mito)

colnames(colData(sce.human.hpc))[9] <- "high.mito"

# Let's also look at low library size and low detected features
qc.lib <- isOutlier(sce.human.hpc$sum, log=TRUE, type="lower")
qc.nexprs <- isOutlier(sce.human.hpc$detected, log=TRUE, type="lower")

discard <- qc.lib | qc.nexprs | high.mito

##Bind library size/expressed feature logical to sce
colData(sce.human.hpc) <- cbind(colData(sce.human.hpc), qc.lib)
colData(sce.human.hpc) <- cbind(colData(sce.human.hpc), qc.nexprs)
colData(sce.human.hpc) <- cbind(colData(sce.human.hpc), discard)

##Plot results of library size/expressed features
pdf("human_HPC_QC_metrics_highmito.pdf", height=5)
grid.arrange(
  plotColData(sce.human.hpc, y="sum", colour_by="discard") +
    scale_y_log10() + ggtitle("Total count"),
  plotColData(sce.human.hpc, y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(sce.human.hpc, y="subsets_Mito_percent",
              colour_by="high.mito") + ggtitle(paste0("Mito % (cutoff = ", mitoCutoffs,")")),
  ncol=3
)

dev.off()

##discard low expressed features/small libraries
sce.human.hpc<-sce.human.hpc[,!sce.human.hpc$discard]
dim(sce.human.hpc)

save(sce.human.hpc, e.out, high.mito, file="sce_humHPC_processing.rda")
### Gene annotation ============================================================
## Add rowRanges data to SCE object
#         (see https://github.com/LieberInstitute/HumanPilot/blob/c8a3a31b991081d656ededee59da45aa0494b334/Analysis/Layer_Notebook.R#L78-L87)
library(rtracklayer)      

## get annotation
map = read.delim("/dcl02/lieber/ajaffe/ErikNelson/ARG_HPC/Br5212/outs/raw_feature_bc_matrix/features.tsv.gz",
                 as.is=TRUE, header=FALSE,col.names=c("EnsemblID", "Symbol", "Type"))
## get GTF
gtf = import("/dcl01/ajaffe/data/lab/singleCell/cellranger_reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf")
## of length 2565061
gtf = gtf[gtf$type == "gene"]
## of length 36601
names(gtf) = gtf$gene_id
table(names(gtf) == map$EnsemblID)
#TRUE
#36601      - (so the reordering below isn't actually necessary, but safe)
gtf = gtf[map$EnsemblID]
seqlevels(gtf)[1:25] = paste0("chr", seqlevels(gtf)[1:25])
mcols(gtf) = mcols(gtf)[,c(5:9)]

table(gtf$gene_type)
#IG_C_gene IG_C_pseudogene       IG_D_gene       IG_J_gene IG_J_pseudogene 
#14               9              37              18               3 
#IG_V_gene IG_V_pseudogene          lncRNA  protein_coding       TR_C_gene 
#144             188           16562           19394               6 
#TR_D_gene       TR_J_gene TR_J_pseudogene       TR_V_gene TR_V_pseudogene 
#4              79               4             106              33

# Then to add:
table(gtf$gene_id == rowData(sce.human.hpc)$ID)  # all 36601 TRUE

rowData(sce.human.hpc) <- cbind(rowData(sce.human.hpc), mcols(gtf))

## Save
save(sce.human.hpc, e.out, high.mito, file="sce_humHPC_processing.rda")

### Normalization===============================================================
##Just to get logcounts assay
##Don't actually use for dimension reduction or feature selection
set.seed(13700)
clusters <- quickCluster(sce.human.hpc)
sce.human.hpc <- computeSumFactors(sce.human.hpc, cluster=clusters)
sce.human.hpc <- logNormCounts(sce.human.hpc)

### Dimensionality reduction/feature selection===============================================
#This is based on fitting Poisson model and then calculating Pearson residuals
#Let's plot highly deviant genes first
sce.human.hpc<-devianceFeatureSelection(sce.human.hpc, assay="counts", 
                                 fam="poisson",sorted=TRUE)
plot(rowData(sce.human.hpc)$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,100000))

##take the top 10000 highly deviant genes
hdg<-rownames(counts(sce.human.hpc))[1:10000]

#Run nullResiduals to calculate Pearson residuals from model
sce.human.hpc<-nullResiduals(sce.human.hpc, assay = "counts",
                      fam = "poisson", type = "pearson")

##Now, PCA on residuals to approximate GLM-PCA
set.seed(1000) 
sce.human.hpc <- runPCA(sce.human.hpc,exprs_values="poisson_pearson_residuals",ncomponents=50, 
                 subset_row=hdg,BSPARAM=BiocSingular::RandomParam()) 
reducedDimNames(sce.human.hpc)

##Plot UMAP
set.seed(100000)
sce.human.hpc <- runUMAP(sce.human.hpc, dimred="PCA",subset_row=hdg)

#Verify length
ncol(reducedDim(sce.human.hpc, "PCA"))

# Save into a new data file, which will dedicate for downstream [prelim] analysis
save(sce.human.hpc, hdg,
     file="sce_humHPC_postdimred.rda")

### Clustering=============================================================
g20 <- buildSNNGraph(sce.human.hpc, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g20)$membership
colData(sce.human.hpc)$label <- factor(clust)
table(colLabels(sce.human.hpc))

#1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
#336 1063  881  241  142  601 1513  353  193  367  870  418 1000  215  195  322 

#17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
#278  192 1620  350   87  352  209  269  531  529  107   75  759   82  262  131 

#33   34   35   36   37   38   39   40   41   42   43   44   45 
#149   59  103   37  352   44  185   99   63   28   43   32   18 

##Let's plot some of our favorite markers from the mouse data
markers.custom = list(
  'neuron' = c('SNAP25','SYT1' , 'GAD2','SLC17A7'),
  #'excitatory_neuron' = c('SLC17A6', 'SLC17A7', 'SLC17A8'),
  #'inhibitory_neuron' = c('GAD1', 'GAD2', 'SLC32A1'),
  #'inhibitory subtypes' = c('SST','CCK','PVALB','VIP'),
  #'DG cell'= c('PROX1','GLIS3','DOCK10','MAML2'),
  #'CA1' = c('FIBCD1','ADGRL2'),
  #'CA2' = c('NTF3','CCDC3'),
  #'CA3' = c('GALNT3','NECTIN3'),
  #'Sub' = c('NTS','FN1'),
  #'Mossy' = c('CALB2','CSF2RB'),
  #'oligodendrocyte' = c('MOBP'),# 'MBP', 'PLP1'),
  'oligodendrocyte_precursor' = c('PDGFRA','VCAN')#,# 'VCAN', 'CSPG4', 'GPR17'),
  #'microglia' = c('C3'),# 'CSF1R', 'C3'),
  #'astrocyte' = c('GFAP')#,# 'AQP4')
  )

###Plots, before MNN correction=================================================
##Plot those markers
pdf("human_hpc_violin_mouse_markers.pdf", height=6, width=12)
for(i in 1:length(markers.custom)){
  print(
    plotExpression(sce.human.hpc, exprs_values = "logcounts", features=c(markers.custom[[i]]),
                   x="label", colour_by="label", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
                   width = 0.3) + 
      ggtitle(label=paste0(names(markers.custom)[i], " markers"))
  )
}
dev.off()

##A LOT of these are really bad...also maybe we also have batch effects?
##Let's plot UMAPs coloured by Sample and by cluster

pdf("human_hpc_UMAP_uncorrected.pdf")
plotUMAP(sce.human.hpc,text_by='label',colour_by='Sample')
dev.off()

pdf("human_hpc_UMAP_uncorrected_byCluster.pdf")
plotUMAP(sce.human.hpc,text_by='label',colour_by='label')
dev.off()

##Oh yeah, definitely need to do MNN correction


###MNN correction for batch effects=============================================
##run fastMNN()
set.seed(1000101001)
mnn.out <- fastMNN(sce.human.hpc, batch=sce.human.hpc$Sample, d=50, k=50, subset.row=hdg,
                   assay.type="logcounts",
                   BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
mnn.out

#class: SingleCellExperiment 
#dim: 10000 15755 
#metadata(2): merge.info pca.info
#assays(1): reconstructed
#rownames(10000): MALAT1 KCNIP4 ... AC023424.3 POLG2
#rowData names(1): rotation
#colnames(15755): AAACCCAAGCGTTCCG-1 AAACCCATCTCGTCAC-1 ...
#TTTGTTGAGGAGTCTG-1 TTTGTTGTCAAGCCAT-1
#colData names(1): batch
#reducedDimNames(1): corrected
#altExpNames(0):

##make UMAP
set.seed(100000)
mnn.out <- runUMAP(mnn.out, dimred="corrected")

##Clustering using the mnn-corrected object
g20 <- buildSNNGraph(mnn.out, k=20, use.dimred = 'corrected')
clust <- igraph::cluster_walktrap(g20)$membership
colData(mnn.out)$label <- factor(clust)
table(colLabels(mnn.out))

#1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
#356  357  606  632  202  100  253 1490   65  243  330  153  291  409  204  977 
#17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
#197 1109   70  373  382 2149  159  427  113   92  121  345  205  551  248  163 
#33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48 
#171  152   84   29   30  147   99   96  866   91   86   46  103   94   91   82 
#49   50   51 
#32   32   52 

g50 <- buildSNNGraph(mnn.out, k=50, use.dimred = 'corrected')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(mnn.out)$k_50_label <- factor(clust50)
table(colData(mnn.out)$k_50_label)

#1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
#584  806 1105  539 2093  334  984  486 3177 1613  678  107  223  212  703  951 
#17   18   19   20   21   22   23   24 
#199  357  203   65   41   77  143   75 

##make new cluster labels in the original SCE object
sce.human.hpc$mnn_20_label<-mnn.out$label
sce.human.hpc$mnn_50_label<-mnn.out$k_50_label

#Save objects
save(mnn.out,sce.human.hpc, file='human_hpc_sce_objects.rda')

###Plots, post MNN batch correction=============================================
pdf("20210811_human_hpc_violin_markers_mnn50.pdf", height=6, width=12)
for(i in 1:length(markers.custom)){
  print(
    plotExpression(sce.human.hpc, exprs_values = "logcounts", features=c(markers.custom[[i]]),
                   x="mnn_50_label", colour_by="mnn_50_label", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
                   width = 0.3) + 
      ggtitle(label=paste0(names(markers.custom)[i], " markers"))
  )
}
dev.off()

pdf("human_hpc_UMAP_corrected20.pdf")
plotUMAP(mnn.out,text_by='label',colour_by='batch')
dev.off()

pdf("human_hpc_UMAP_corrected50.pdf")
plotUMAP(mnn.out,text_by='k_50_label',colour_by='batch')
dev.off()

pdf("human_hpc_UMAP_corrected20_byCluster.pdf")
plotUMAP(mnn.out,text_by='label',colour_by='label')
dev.off()

pdf("human_hpc_UMAP_corrected50_byCluster.pdf")
plotUMAP(mnn.out,text_by='k_50_label',colour_by='k_50_label')
dev.off()

###Neurons analysis=============================================================
##First, make a new SCE object with only SNAP25+ clusters
sce.human.hpc$neuron<-ifelse(
  sce.human.hpc$mnn_50_label %in% c(1,2,4,5,6,7,13,14,
                                    15,17,18,19,20,23),
                                    T,F
  )

sce.human.neuron<-sce.human.hpc[,sce.human.hpc$neuron==TRUE]
sce.human.neuron

### Normalization, neurons===============================================================
##Just to get logcounts assay
##Don't actually use for dimension reduction or feature selection
set.seed(13700)
clusters <- quickCluster(sce.human.neuron)
sce.human.neuron <- computeSumFactors(sce.human.neuron, cluster=clusters)
sce.human.neuron <- logNormCounts(sce.human.neuron)

### Dimensionality reduction/feature selection, neurons===============================================
#This is based on fitting Poisson model and then calculating Pearson residuals
#Let's plot highly deviant genes first
sce.human.neuron<-devianceFeatureSelection(sce.human.neuron, assay="counts", 
                                           fam="poisson",sorted=TRUE)
plot(rowData(sce.human.neuron)$poisson_deviance, type="l", xlab="ranked genes",
     ylab="poisson deviance", main="Feature Selection with Deviance",ylim=c(0,100000))

##take the top 10000 highly deviant genes
hdg<-rownames(counts(sce.human.neuron))[1:10000]

#Run nullResiduals to calculate Pearson residuals from model
sce.human.neuron<-nullResiduals(sce.human.neuron, assay = "counts",
                                fam = "poisson", type = "pearson")

##Now, PCA on residuals to approximate GLM-PCA
set.seed(1000) 
sce.human.neuron <- runPCA(sce.human.neuron,exprs_values="poisson_pearson_residuals",ncomponents=50, 
                           subset_row=hdg,BSPARAM=BiocSingular::RandomParam()) 
reducedDimNames(sce.human.neuron)

##Plot UMAP
set.seed(100000)
sce.human.neuron <- runUMAP(sce.human.neuron, dimred="PCA",subset_row=hdg)

#Verify length
ncol(reducedDim(sce.human.neuron, "PCA"))

# Save into a new data file, which will dedicate for downstream [prelim] analysis
save(sce.human.neuron, hdg,
     file="sce_humHPC_neurons_postdimred.rda")

### Clustering, neurons=============================================================
g20 <- buildSNNGraph(sce.human.neuron, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g20)$membership
colData(sce.human.neuron)$label <- factor(clust)
table(colLabels(sce.human.neuron))

##Let's plot some of our favorite markers from the mouse data
markers.custom = list(
  'neuron' = c('SNAP25','SYT1','SLC17A7','GAD2'),
  #'excitatory_neuron' = c('SLC17A6', 'SLC17A7', 'SLC17A8'),
  #'inhibitory_neuron' = c('GAD1', 'GAD2', 'SLC32A1'),
  'inhibitory subtypes' = c('SST','VIP'),
  'DG cell'= c('PROX1','GLIS3','DOCK10','MAML2'),
  'CA1' = c('FIBCD1','ADGRL2','COL5A2','CLMP')
  #'CA2' = c('NTF3','CCDC3'),
  #'CA3' = c('GALNT3','NECTIN3'),
  #'Sub' = c('NTS','FN1'),
  #'Mossy' = c('CALB2','CSF2RB'),
  #'oligodendrocyte' = c('MOBP'),# 'MBP', 'PLP1'),
  #'oligodendrocyte_precursor' = c('PDGFRA'),# 'VCAN', 'CSPG4', 'GPR17'),
  #'microglia' = c('C3'),# 'CSF1R', 'C3'),
  #'astrocyte' = c('GFAP')#,# 'AQP4')
)

###Plots, before MNN correction, neurons=================================================
##Plot those markers
pdf("human_hpc_neurons_violin_mouse_markers.pdf", height=6, width=12)
for(i in 1:length(markers.custom)){
  print(
    plotExpression(sce.human.neuron, exprs_values = "logcounts", features=c(markers.custom[[i]]),
                   x="label", colour_by="label", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
                   width = 0.3) + 
      ggtitle(label=paste0(names(markers.custom)[i], " markers"))
  )
}
dev.off()

##A LOT of these are really bad...also maybe we also have batch effects?
##Let's plot UMAPs coloured by Sample and by cluster

pdf("human_hpc_neurons_UMAP_uncorrected.pdf")
plotUMAP(sce.human.neuron,text_by='label',colour_by='Sample')
dev.off()

pdf("human_hpc_neurons_UMAP_uncorrected_byCluster.pdf")
plotUMAP(sce.human.neuron,text_by='label',colour_by='label')
dev.off()

##Oh yeah, definitely need to do MNN correction


###MNN correction for batch effects, neurons=============================================
##run fastMNN()
set.seed(1000101001)
mnn.out <- fastMNN(sce.human.neuron, batch=sce.human.neuron$Sample, d=50, k=50, subset.row=hdg,
                   assay.type="logcounts",
                   BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
mnn.out

##make UMAP
set.seed(100000)
mnn.out <- runUMAP(mnn.out, dimred="corrected")

##Clustering using the mnn-corrected object
#g20 <- buildSNNGraph(mnn.out, k=20, use.dimred = 'corrected')
#clust <- igraph::cluster_walktrap(g20)$membership
#colData(mnn.out)$label <- factor(clust)
#table(colLabels(mnn.out))


g50 <- buildSNNGraph(mnn.out, k=50, use.dimred = 'corrected')
clust50 <- igraph::cluster_walktrap(g50)$membership
colData(mnn.out)$k_50_label <- factor(clust50)
table(colData(mnn.out)$k_50_label)

##make new cluster labels in the original SCE object
#sce.human.neuron$mnn_20_label<-mnn.out$label
sce.human.neuron$mnn_50_label<-mnn.out$k_50_label

##look at distribution of transcripts by cluster
cellIdx <- splitit(sce.human.neuron$mnn_50_label)
sapply(cellIdx, function(x){quantile(sce.human.neuron[ ,x]$sum)})

##Make data frame with transcripts and clusters columns
data<-data.frame("transcripts"=sce.human.neuron$sum,"cluster"=as.numeric(sce.human.neuron$mnn_50_label))

##boxplot of transcripts by cluster
pdf("human_hpc_neurons_boxplot_transcriptsPerCluster.pdf", height=5)
boxplot(data$transcripts ~ data$cluster,main="Transcripts per cluster",xlab='Cluster',ylab='Transcripts') +
dev.off()

#Save objects
save(mnn.out,sce.human.neuron, file='humHPC_sce_objects_neuron.rda')

###Plots, post MNN batch correction, neurons=============================================
pdf("human_hpc_neurons_violin_markers_mnn50.pdf", height=6, width=12)
for(i in 1:length(markers.custom)){
  print(
    plotExpression(sce.human.neuron, exprs_values = "logcounts", features=c(markers.custom[[i]]),
                   x="mnn_50_label", colour_by="mnn_50_label", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
                   width = 0.3) + 
      ggtitle(label=paste0(names(markers.custom)[i], " markers"))
  )
}
dev.off()

mnn.out$level_1<-factor(
  ifelse(mnn.out$k_50_label %in% c(1,4,5,8,15,16),"Inhib","Excit")
)
pdf("human_hpc_neurons_UMAP_corrected20.pdf")
plotUMAP(mnn.out,text_by='label',colour_by='batch')
dev.off()

pdf("human_hpc_neurons_UMAP_corrected50.pdf")
plotUMAP(mnn.out,text_by='k_50_label',colour_by='batch')
dev.off()

pdf("human_hpc_neurons_UMAP_corrected20_byCluster.pdf")
plotUMAP(mnn.out,text_by='label',colour_by='label')
dev.off()

pdf("human_hpc_neurons_UMAP_corrected50_byCluster.pdf")
plotUMAP(mnn.out,text_by='k_50_label',colour_by='k_50_label')
dev.off()

###Marker gene detection, neurons========================================================
markers<-findMarkers(sce.human.neuron,pval.type='all',direction='up',
                     groups=sce.human.neuron$mnn_50_label)
markers
save(markers,file='human_hpc_neurons_markers.rda')

##Make list of top 5 markers for each cluster
for(i in 1:length(markers.list)){
  markers.list[[i]]<-rownames(markers[[i]][1:5,])
}

##Make one list with all markers
features<-do.call('c',markers.list)

##removing a few duplicates
s3e.human.neuron<-s3e.human.neuron[,unique(colnames(s3e.human.neuron))]

pdf("human_hpc_neurons_markerHeatmap_draft.pdf")
plotHeatmap(s3e.human.neuron,features=features,order_columns_by="mnn_50_label",
            colour_columns_by="mnn_50_label",main="Marker genes",cluster_rows=F)
dev.off()


