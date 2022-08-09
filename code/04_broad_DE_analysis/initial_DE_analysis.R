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
load("/users/enelson/sceSubset.rda")


##Differential expression analysis, broad===================================
###Set up data==================================================================
sce.subset$cellType_DE<-as.character(sce.subset$annotation)
sce.subset$cellType_DE<-factor(ifelse(sce.subset$cellType_DE %in% c("DG.1","DG.2"),
                                      "DG",sce.subset$cellType_DE))
#aggregateAcrossCells to pseudobulk
summed <- aggregateAcrossCells(sce.subset,
                               ids=colData(sce.subset)
                               [c('Sample')])

#Filter out lowly expressed genes 


#make DGElist
y <- DGEList(counts(summed), samples=colData(summed),
             genes=rowData(summed))
y$samples$condition<-factor(y$samples$condition)
y<-y[filterByExpr(y, group=summed$Sample),]

#make design matrix
design <- model.matrix(~condition,y$samples)

#normalization 
y <- calcNormFactors(y)

#estimate dispersion parameter
y <- estimateDisp(y, design = design)

#run DE analysis and look at results
efit <- glmQLFit(y, design = design)
elrt <- glmQLFTest(efit, coef=2)
res<-topTags(elrt,n=Inf)
x<-as.data.frame(res)
x<-x[x$FDR<=0.05,]

##MA plot for figure 2
ids <- c("Nptx2" , 'Bdnf', 'Gadd45b', 'Arc')
gene.labels <- x[x$Symbol %in% ids,]
#pdf('overall_MA_plot.pdf',height=5,width=5)
plotSmear(elrt, de.tags = rownames(x))
text(x=gene.labels$logCPM, y=gene.labels$logFC,
     labels=rownames(gene.labels), cex=1, pos=1)

dev.off()

#grouped heatmap for figure 2
pdf('heatmap_IEGs_bySample.pdf',w=4,h=3)
features=c('Bdnf','Nptx2','Baz1a','Inhba','Mapk4','Arc')
plotGroupedHeatmap(sce.subset,features=features,cluster_rows=F,cluster_cols=F,group='Sample',center=T,fontsize=11)
dev.off()

##volcano plot

data <- data %>% 
  mutate(
    Expression = case_when(logFC >= 0 & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= 0 & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
p2 <- ggplot(data, aes(logFC, -log(PValue,10))) +
  geom_point(aes(color = Expression), size = 1/2) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"PValue")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

p2
pdf('volcano.pdf',h=3,w=4.25)
p2 + geom_label_repel(data = top_genes,
                      mapping = aes(logFC, -log(PValue,10),
                      label=rownames(top_genes))) + 
                      theme(legend.position='none')


#GO analysis for figure 2
enrich_go <- enrichGO(gene = x$ID,
                      OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "BP",
                      pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)

enrich_up <- enrichGO(gene = x$ID[x$logFC > 0],
                      OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "BP",
                      pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)

#enrich_down2 <- enrichGO(gene = down$ID,
#                      OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "MF",
#                      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

## Visualize enrichment results
pdf('barplot_GOanalysis_10cats.pdf',h=6,w=4.25)
barplot(enrich_up, font.size = 11,showCategory=10)
dev.off()

hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
up_mapping <- getBM(
  attributes = c('entrezgene_id'), 
  filters = 'ensembl_gene_id',
  values = up$ID,
  mart = hsmart
)

save(sce.subset, hdg, pc.choice.hpc,marks,file="/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50_clusteredPCs.rda")
save(summed,y,x,res,enrich_go,file='/users/enelson/mouseHPC_broad_differentialExpressionAnalysis_results.rda')

sce.subset$cellType<-as.character(sce.subset$annotation)
sce.subset$cellType<-ifelse(sce.subset$cellType %in% c('RHP.L5/6/6b.1','RHP.L5/6/6b.2','RHP.L5/6/6b.3','RHP.L5/6/6b.4'),'RHP.L5/6/6b',
ifelse(sce.subset$cellType %in% c('RHP.L2/3.1','RHP.L2/3.2','RHP.L2/3.3','RHP.L2/3.4','RHP.L2/3.5','RHP.L2/3.6'),'RHP.L2/3',sce.subset$cellType))
sce.subset$cellType<-as.factor(sce.subset$cellType)

