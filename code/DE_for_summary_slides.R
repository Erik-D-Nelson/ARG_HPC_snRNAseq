###DG===========================================================================
#Subset for DG only and aggregateAcrossCells to pseudobulk
summed<-sce.subset[,sce.subset$level_6=='DG']
summed <- aggregateAcrossCells(sce.subset, 
                               ids=colData(sce.subset)[,c("sample_name", "annotated_cell_type")])

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
y <- calcNormFactors(y)

#make design matrix
design <- model.matrix(~factor(condition)+factor(annotated_cell_type), y$samples)

#Run cqn
#cqn.summed <- cqn(counts(summed), x=rowData(summed)$gc_content,
#                  lengths=rowData(summed)$length,
#                  sizeFactors = summed$sum,
#                  verbose = TRUE)
#y$offset <- cqn.summed$glm.offset

#estimate dispersion parameter
y <- estimateDisp(y, design)

#run DE analysis and look at results
efit <- glmQLFit(y, design, robust=TRUE)
elrt <- glmLRT(efit, coef = 2)
res_DG<-topTags(elrt,n=Inf)
FDR_filtered_DG<-as.data.frame(res_DG)
FDR_filtered_DG<-FDR_filtered_DG[FDR_filtered_DG$FDR<=0.05,]

#Get DEGs
FDR_filtered_DG<-as.data.frame(res_DG)
FDR_filtered_DG<-FDR_filtered_DG[FDR_filtered_DG$FDR<=0.05,]
DG_DEGs<-rownames(FDR_filtered_DG)

#Get DE IEGs
DG_DE_IEGs<-intersect(DG_DEGs,arg_list)

#MA plot 
pdf('DG_MA_plot.pdf')
plotSmear(elrt, de.tags = DG_DEGs)
dev.off()



###GABA===========================================================================
#Subset for GABA only and aggregateAcrossCells to pseudobulk
summed<-sce.subset[,sce.subset$level_6=='GABA']
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

#Run cqn
cqn.summed <- cqn(counts(summed), x=rowData(summed)$gc_content,
                  lengths=rowData(summed)$length,
                  sizeFactors = summed$sum,
                  verbose = TRUE)
y$offset <- cqn.summed$glm.offset

#estimate dispersion parameter
y <- estimateGLMCommonDisp(y, design = design)

#run DE analysis and look at results
efit <- glmFit(y, design = design)
elrt <- glmLRT(efit, coef = 2)
res_GABA<-topTags(elrt,n=Inf)
FDR_filtered_GABA<-as.data.frame(res_GABA)
FDR_filtered_GABA<-FDR_filtered_GABA[FDR_filtered_GABA$FDR<=0.05,]

#Get DEGs
FDR_filtered_GABA<-as.data.frame(res_GABA)
FDR_filtered_GABA<-FDR_filtered_GABA[FDR_filtered_GABA$FDR<=0.05,]
GABA_DEGs<-rownames(FDR_filtered_GABA)

#Get DE IEGs
GABA_DE_IEGs<-intersect(GABA_DEGs,arg_list)

#MA plot 
pdf('GABA_MA_plot.pdf')
plotSmear(elrt, de.tags = GABA_DEGs)
dev.off()
###Pyramidal===========================================================================
#Subset for Pyramidal only and aggregateAcrossCells to pseudobulk
summed<-sce.subset[,sce.subset$level_6=='Pyramidal']
summed <- aggregateAcrossCells(summed, 
                               ids=c(colData(summed)$sample_name,colData(summed))

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

#Run cqn
cqn.summed <- cqn(counts(summed), x=rowData(summed)$gc_content,
                  lengths=rowData(summed)$length,
                  sizeFactors = summed$sum,
                  verbose = TRUE)
y$offset <- cqn.summed$glm.offset

#estimate dispersion parameter
y <- estimateGLMCommonDisp(y, design = design)

#run DE analysis and look at results
efit <- glmFit(y, design = design)
elrt <- glmLRT(efit, coef = 2)
res_Pyramidal<-topTags(elrt,n=Inf)
FDR_filtered_Pyramidal<-as.data.frame(res_Pyramidal)
FDR_filtered_Pyramidal<-FDR_filtered_Pyramidal[FDR_filtered_Pyramidal$FDR<=0.05,]

#Get DEGs
FDR_filtered_Pyramidal<-as.data.frame(res_Pyramidal)
FDR_filtered_Pyramidal<-FDR_filtered_Pyramidal[FDR_filtered_Pyramidal$FDR<=0.05,]
Pyramidal_DEGs<-rownames(FDR_filtered_Pyramidal)

#Get DE IEGs
Pyramidal_DE_IEGs<-intersect(Pyramidal_DEGs,arg_list)

#MA plot 
pdf('Pyramidal_MA_plot.pdf')
plotSmear(elrt, de.tags = Pyramidal_DEGs)
dev.off()



de.specific <- pseudoBulkDGE(summed,
                                  label=summed$annotation,
                                  design=~factor(condition),
                                  coef=2,
                                  condition=summed$condition
)

FDR_CA1<-as.data.frame(de.specific[[1]])
FDR_CA1<-FDR_CA1[FDR_CA1$FDR<=0.05,]
FDR_CA1<-FDR_CA1[!is.na(FDR_CA1$FDR),]
CA1_DEGs<-rownames(FDR_CA1)

