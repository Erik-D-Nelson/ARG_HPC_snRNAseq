#######################################
###Comparison of DE methods, DG only###
#######################################
library(edgeR)
library(cqn)
library(UpSetR)
library(SingleCellExperiment)
library(scRNAseq)
library(GGally)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)
library(scRNAseq)
library(EnsDb.Mmusculus.v79)
library(scater)
library(scran)
library(jaffelab)
load("/users/enelson/sce_subset_take2.rda")

###Set up data==================================================================
#Subset for DG only and aggregateAcrossCells to pseudobulk
summed<-sce.subset[,sce.subset$level_6=='DG']
summed <- aggregateAcrossCells(summed, 
                               ids=colData(summed)$sample_name)

#Use perCellQCMetrics to get library sizes for each sample 
#(didn't carry over after pseudobulk, may be an easier way but this worked
#so I ran with it)
location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(summed)$GENEID, 
                   column="SEQNAME", keytype="GENEID")
is.mito <- which(location=="MT")
stats <- perCellQCMetrics(summed, subsets=list(Mito=is.mito))
colnames(stats)
colData(summed)$sum<-stats$sum

#Filter out lowly expressed genes and remove genes with length=NA
summed<-summed[filterByExpr(summed, group=summed$condition),]
summed<-summed[!is.na(rowData(summed)$length),]

#make DGElist
y <- DGEList(counts(summed), samples=colData(summed),
             genes=rowData(summed))
y$samples$condition<-factor(y$samples$condition)

#make design matrix
design <- model.matrix(~condition, y$samples)


###Run cqn, default======================================================================
cqn.summed <- cqn(counts(summed), x=rowData(summed)$gc_content,
                  lengths=rowData(summed)$length,
                  sizeFactors = summed$sum,
                  verbose = TRUE)
###Run cqn, lengthMethod=fixed==================================================
cqn.fixed <- cqn(counts(summed), x=rowData(summed)$gc_content,
                  lengths=rowData(summed)$length,
                  sizeFactors = summed$sum,
                  verbose = TRUE,
                  lengthMethod="fixed")
###A note on how I used cqn for length/gc only==================================

#Normally, cqn adjusts for length AND removes the systematic influence of a 
#covariate x on counts. X is usually gc content. Both variables are required.
#However, cqn vignette says: "It is also possible to just fit a model like 
#log2 (RPKM) = s(x). In this model gene length is included as a known offset. 
#This is done by using the cqn(lengthMethod = "fixed"). If this is done, and 
#lengths is equal to 1000, it is equivalent to not using gene length at all."
#SO, to adjust for length ONLY, I set lengthMethod="fixed" and x=length. I did 
#the same to adjust for gc only, except I set x=gc.
#My rationale was this effectively eliminates the lengths argument, but allows
#me to adjust for either gc or length individually.

###Run cqn, length only==================================================
cqn.length <- cqn(counts(summed), x=rowData(summed)$length,
                  lengths=1000,
                  sizeFactors = summed$sum,
                  verbose = TRUE,
                  lengthMethod="fixed")
###Run cqn, GC only==================================================
cqn.gc <- cqn(counts(summed), x=rowData(summed)$gc_content,
                  lengths=1000,
                  sizeFactors = summed$sum,
                  verbose = TRUE,
                  lengthMethod="fixed")
###GLMfit with cqn==============================================================
#set DGEList offset to offset calculated with cqn
y1<-y
y1$offset <- cqn.summed$glm.offset
y.cqn <- estimateGLMCommonDisp(y1, design = design)

#run DE analysis
efit <- glmFit(y.cqn, design = design)
elrt <- glmLRT(efit, coef = 2)
res1<-topTags(elrt,n=Inf)
x1<-as.data.frame(res1)
x1<-x1[x1$FDR<=0.05,]

###limma, no length normalization===============================================
y2<-y

#normalization 
y2 <- calcNormFactors(y2)

#run DE analysis and look at results
v <- voomWithQualityWeights(y2, design)
fit <- lmFit(v)
fit <- eBayes(fit, robust=TRUE)
res2 <- topTable(fit, sort.by="P", n=Inf)
x2<-as.data.frame(res2)
x2<-x2[x2$adj.P.Val<=0.05,]


###GLMfit, no length normalization===============================================================
y3<-y

#normalization 
y3 <- calcNormFactors(y3)

#estimate dispersion parameter
y3 <- estimateGLMCommonDisp(y3, design = design)

#run DE analysis and look at results
efit <- glmFit(y3, design = design)
elrt <- glmLRT(efit, coef = 2)
res3<-topTags(elrt,n=Inf)
x3<-as.data.frame(res3)
x3<-x3[x3$FDR<=0.05,]

###limma,with cqn===============================================================
y4<-y

#set DGEList offset to offset calculated by cqn
y4$offset <- cqn.summed$offset

#run DE analysis and look at results
v <- voomWithQualityWeights(y4, design)
fit <- lmFit(v)
fit <- eBayes(fit, robust=TRUE)
res4 <- topTable(fit, sort.by="P", n=Inf)
x4<-as.data.frame(res4)
x4<-x4[x4$adj.P.Val<=0.05,]

###GLMfit with cqn lengthMethod=fixed=======================================================
y5<-y

#set DGEList offset to offset calculated with cqn
y5$offset <- cqn.fixed$glm.offset
y.cqn <- estimateGLMCommonDisp(y5, design = design)

#run DE analysis
efit <- glmFit(y.cqn, design = design)
elrt <- glmLRT(efit, coef = 2)
res5<-topTags(elrt,n=Inf)
x5<-as.data.frame(res5)
x5<-x5[x5$FDR<=0.05,]

###GLMfit with cqn, length only=======================================================
y6<-y

#set DGEList offset to offset calculated with cqn
y6$offset <- cqn.length$glm.offset
y.cqn <- estimateGLMCommonDisp(y6, design = design)

#run DE analysis
efit <- glmFit(y.cqn, design = design)
elrt <- glmLRT(efit, coef = 2)
res6<-topTags(elrt,n=Inf)
x6<-as.data.frame(res6)
x6<-x6[x6$FDR<=0.05,]

###GLMfit with cqn, gc only=======================================================
y7<-y

#set DGEList offset to offset calculated with cqn
y7$offset <- cqn.gc$glm.offset
y.cqn <- estimateGLMCommonDisp(y7, design = design)

#run DE analysis
efit <- glmFit(y.cqn, design = design)
elrt <- glmLRT(efit, coef = 2)
res7<-topTags(elrt,n=Inf)
x7<-as.data.frame(res7)
x7<-x7[x7$FDR<=0.05,]

###make UpSetR plot=============================================================
#make list of DEGs for first 4 methods
deg_list<-list(
  'GLM_with_cqn'= rownames(x1),
  'GLM_no_cqn'= rownames(x3),
  'limma_with_cqn' = rownames(x4),
  'limma_no_cqn' = rownames(x2)
)

#First Upset plot
#pdf('2020_02_12_upsetplot_DEmethods.pdf', height=10,width=15)
upset(fromList(deg_list),order.by='freq')
#dev.off()

#Make list of DEGs for edgeR methods only
deg_list2<-list(
  'GLM_cqn'= rownames(x1),
  'GLM_no_cqn'= rownames(x3),
  'GLM_cqn_fixed' = rownames(x5),
  'GLM_cqn_length_only' = rownames(x6),
  'GLM_cqn_gc_only' = rownames(x7)
)

#make second Upset plot
#pdf('2020_02_12_upsetplot2_DEmethods.pdf', height=10,width=15)
upset(fromList(deg_list2),order.by='freq')
#dev.off()

#make list of DEGs for all methods used
deg_list3<-list(
  'GLM_cqn'= rownames(x1),
  'GLM_no_cqn'= rownames(x3),
  'GLM_cqn_fixed' = rownames(x5),
  'GLM_cqn_length_only' = rownames(x6),
  'GLM_cqn_gc_only' = rownames(x7),
  'limma_cqn' = rownames(x4),
  'limma_no_cqn' = rownames(x2)
)

#make 3rd Upset plot
pdf('2020_02_12_upsetplot3_DEmethods.pdf',height=10,width=15)
upset(fromList(deg_list3),nsets=7,order.by='freq')
dev.off()

###Make pairwise scatter plots==================================================
# match rownames so we can make a data frame with logFCs and FDR/adjusted P for 
# each method 
# (when I clean up code I'll make these functions so it's less redundant)
res2<-res2[match(rownames(res1),rownames(res2)),]
res3<-res3[match(rownames(res1),rownames(res3)),]
res4<-res4[match(rownames(res1),rownames(res4)),]
res5<-res5[match(rownames(res1),rownames(res5)),]
res6<-res6[match(rownames(res1),rownames(res6)),]
res7<-res7[match(rownames(res1),rownames(res7)),]

#coerce to data frame
# (will also make this a function to eliminate redundancy)
res1<-as.data.frame(res1)
res2<-as.data.frame(res2)
res3<-as.data.frame(res3)
res4<-as.data.frame(res4)
res5<-as.data.frame(res5)
res6<-as.data.frame(res6)
res7<-as.data.frame(res7)

#make data frame with logFCs from each method
df<-data.frame(GLM_cqn_logFC=res1$logFC,
               GLM_nocqn_logFC=res3$logFC,
               limma_cqn_logFC=res4$logFC,
               limma_nocqn_logFC=res2$logFC,
               GLM_cqn_fixed_logFC=res5$logFC,
               GLM_cqn_length_logFC=res6$logFC,
               GLM_cqn_gc_logFC=res7$logFC,
               length=res1$length,
               genes=rownames(res1))

#get DE genes exclusive to and shared between edgeR with/without cqn offset
DE_GLMcqn_only<-setdiff(rownames(x1),rownames(x3))
DE_GLMnocqn_only<-setdiff(rownames(x3),rownames(x1))
DE_both_GLM_methods<-intersect(rownames(x1),rownames(x3))

#make DE genes factor in order to add colors to scatter plots
df$DE<-factor(
  ifelse(df$genes %in% DE_GLMcqn_only,"GLM_cqn_only",
         ifelse(df$genes %in% DE_GLMnocqn_only,"GLM_no_cqn_only",
                ifelse(df$genes %in% DE_both_GLM_methods, "both",
                       "neither")))
)

#same as above--get DE genes exclusive to edgeR with/without cqn and limma with cqn
#also get DE genes shared between these 3 methods
#didn't include limma without cqn offset...according to upset plots, there are
# no DE genes unique to limma without cqn only
DE_GLMcqn_only<-setdiff(rownames(x1),c(rownames(x4),rownames(x3)))
DE_GLMnocqn_only<-setdiff(rownames(x3),c(rownames(x4),rownames(x1)))
DE_limmacqn_only<-setdiff(rownames(x4), c(rownames(x1),rownames(x3)))
DE_all<-intersect(rownames(x1),intersect(rownames(x4),rownames(x3)))

#same as above--make DE genes factor to add colors to scatter plots
df$DE2<-factor(
  ifelse(df$genes %in% DE_GLMcqn_only,"GLM_cqn_only",
         ifelse(df$genes %in% DE_GLMnocqn_only,"GLM_no_cqn_only",
              ifelse(df$genes %in% DE_limmacqn_only, "limma_cqn_only",
                ifelse(df$genes %in% DE_all, "all",
                       "none"))))
)

#make 4 different logFC pairwise scatter plots:
#1 is first 4 methods (edgeR with/without cqn and limma with/without cqn,
# colored by first DE factor (DE genes exclusive to edgeR with cqn, DE genes 
# exclusive to edgeR without cqn, DE genes shared by both, and non-DE genes)
#2 is all methods, colored by the first DE factor
#3 is first 4 methods colored by second DE factor (DE genes exclusive to 
#edgeR with cqn, DE genes exclusive to edgeR without cqn, DE genes exclusive to 
#limma with cqn, DE genes shared by all, and non-DE genes)
#4 is all methods colored by second DE factor
pdf('2020_02_14_pairwise_logFC_DEmethods.pdf',width=10)
ggpairs(df,columns=1:4,aes(colour=df$DE,size=df$length))
dev.off()

pdf('2020_02_14_pairwise_logFC_DEmethods2.pdf',width=15)
ggpairs(df,columns=1:7,aes(colour=df$DE,size=df$length))
dev.off()

pdf('2020_02_14_pairwise_logFC_DEmethods3.pdf',width=10)
ggpairs(df,columns=1:4,aes(colour=df$DE2,size=df$length))
dev.off()

pdf('2020_02_14_pairwise_logFC_DEmethods4.pdf',height=10,width=15)
ggpairs(df,columns=1:7,aes(colour=df$DE2,size=df$length))
dev.off()


#make data frame with FDR/adjusted P values from each method
df2<-data.frame(GLM_cqn_FDR=res1$FDR,
               GLM_nocqn_FDR=res3$FDR,
               limma_cqn_adjustedP=res4$adj.P.Val,
               limma_nocqn_adjustedP=res2$adj.P.Val,
               GLM_cqn_fixed_FDR=res5$FDR,
               GLM_cqn_length_FDR=res6$FDR,
               GLM_cqn_gc_FDR=res7$FDR,
               length=res1$length,
               genes=rownames(res1),
               DE=df$DE,
               DE2=df$DE2)

#make 4 different FDR/adjusted P value pairwise scatter plots:
#1 is first 4 methods (edgeR with/without cqn and limma with/without cqn,
# colored by first DE factor (DE genes exclusive to edgeR with cqn, DE genes 
# exclusive to edgeR without cqn, DE genes shared by both, and non-DE genes)
#2 is all methods, colored by the first DE factor
#3 is first 4 methods colored by second DE factor (DE genes exclusive to 
#edgeR with cqn, DE genes exclusive to edgeR without cqn, DE genes exclusive to 
#limma with cqn, DE genes shared by all, and non-DE genes)
#4 is all methods colored by second DE factor
pdf('2020_02_14_pairwise_FDR_DEmethods.pdf',width=10)
ggpairs(df2,columns=1:4,aes(colour=df2$DE,size=df2$length))
dev.off()

pdf('2020_02_14_pairwise_FDR_DEmethods2.pdf',height=10,width=15)
ggpairs(df2,columns=1:7,aes(colour=df2$DE,size=df2$length))
dev.off()

pdf('2020_02_14_pairwise_FDR_DEmethods3.pdf',width=10)
ggpairs(df2,columns=1:4,aes(colour=df2$DE2,size=df2$length))
dev.off()

pdf('2020_02_14_pairwise_FDR_DEmethods4.pdf',height=10,width=15)
ggpairs(df2,columns=1:7,aes(colour=df2$DE2,size=df2$length))
dev.off()

