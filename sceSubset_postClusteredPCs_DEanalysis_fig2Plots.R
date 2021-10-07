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
library(edgeR)
library(org.Mm.eg.db)
library(clusterProfiler)
load("/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50_clusteredPCs.rda")

##UMAP by condition for figure 2
pdf('20210913_UMAP_bySample.pdf',h=7,w=7)
plotUMAP(sce.subset,colour_by='sample_name',text_by='annotation',text_size=3,add_legend=F) 
dev.off()

##Differential expression analysis, broad===================================
###Set up data==================================================================
#Subset for DG only and aggregateAcrossCells to pseudobulk
summed <- aggregateAcrossCells(sce.subset,
                               ids=colData(sce.subset)
                               [c('condition','sample_name','annotation')])

#Use perCellQCMetrics to get library sizes for each sample 
#(didn't carry over after pseudobulk, may be an easier way but this worked
#so I ran with it)
location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(summed)$ID, 
                   column="SEQNAME", keytype="GENEID")
stats <- perCellQCMetrics(summed)
colnames(stats)
colData(summed)$sum<-stats$sum

#filter out small sample sizes
summed <- summed[,summed$ncells >= 20]

#Filter out lowly expressed genes 
summed<-summed[filterByExpr(summed, group=summed$sample_name),]

#make DGElist
y <- DGEList(counts(summed), samples=colData(summed),
             genes=rowData(summed))
y$samples$condition<-factor(y$samples$condition)

#make design matrix
design <- model.matrix(~condition+annotation, y$samples)

#normalization 
y <- calcNormFactors(y)

#estimate dispersion parameter
y <- estimateDisp(y, design = design)

#run DE analysis and look at results
efit <- glmQLFit(y, design = design)
elrt <- glmQLFTest(efit, coef = 2)
res<-topTags(elrt,n=Inf)
x<-as.data.frame(res)
x<-x[x$FDR<=0.05,]

##MA plot for figure 2
ids <- c("Nptx2" , 'Bdnf' , 'Egr3')
gene.labels <- x[x$Symbol %in% ids,]
pdf('202109214overall_MA_plot.pdf',height=5,width=5)
plotSmear(elrt, de.tags = rownames(x))
text(x=gene.labels$logCPM, y=gene.labels$logFC,
     labels=rownames(gene.labels), cex=1, pos=1)

dev.off()

#grouped heatmap for figure 2
pdf('heatmap_IEGs_bySample.pdf')
features=c('Arc','Bdnf','Nptx2','Jun','Junb','Nr4a1','Baz1a')
plotGroupedHeatmap(sce.subset,features=features,cluster_rows=F,cluster_cols=F,group='sample_name',center=F)
dev.off()

#GO analysis for figure 2
enrich_go <- enrichGO(gene = x$ID[x$FDR  < 0.05],
                      OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "BP",
                      pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)

## Visualize enrichment results
pdf('barplot_GOanalysis_30cats.pdf')
barplot(enrich_go, font.size = 7,showCategory=30)
dev.off()

save(sce.subset, hdg, pc.choice.hpc,marks,file="/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50_clusteredPCs.rda")
save(summed,y,x,res,enrich_go,file='/users/enelson/mouseHPC_broad_differentialExpressionAnalysis_results.rda')
