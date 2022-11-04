library(SingleCellExperiment)
library(EnsDb.Mmusculus.v79)
library(scater)
library(scran)
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
sce.subset$cellType_DE<-factor(ifelse(sce.subset$cellType_DE %in% c("GC.1","GC.2"),
                                      "GC",sce.subset$cellType_DE))
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
data<-as.data.frame(res)
data <- data %>% 
  mutate(
    Expression = case_when(logFC >= 0 & FDR <= 0.05 ~ "Upregulated",
                           logFC <= 0 & FDR <= 0.05 ~ "Downregulated",
                          # logFC<1 & logFC>0 & FDR <= 0.05 ~ "Upregulated, FDR only",
                          # logFC<0 & logFC>-1 & FDR <= 0.05 ~ "Downregulated, FDR only",
                           TRUE ~ "Unchanged")
  )

top_genes <- bind_rows(
  data %>% 
    filter(Expression == 'Upregulated') %>% 
    arrange(FDR, desc(abs(logFC))) %>% 
    head(5),
  data %>% 
    filter(Expression == 'Downregulated') %>% 
    arrange(FDR, desc(abs(logFC))) %>% 
    head(5)
)

p2 <- ggplot(data, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = Expression), size = 1/2) +
  xlab(expression("log"[2]*"fold change")) + 
  ylab(expression("-log"[10]*"adjusted p-value")) +
  scale_color_manual(values = c(#"darkmagenta",
                                "indianred2", 
                                "gray50", 
                               # "seagreen", 
                                "cornflowerblue")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position='right') +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
 # geom_vline(xintercept = 1, linetype = "dashed")+
#  geom_vline(xintercept = -1, linetype = "dashed")+
  geom_text_repel(data = top_genes,
                  mapping = aes(logFC, -log(FDR,10),
                  color=Expression,
                  label=rownames(top_genes)))

pdf('volcano.pdf',h=7,w=8.5)
p2
dev.off()

#GO analysis for figure 2
enrich_go <- enrichGO(gene = x$ID,
                      OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "BP",
                      pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)


## Visualize enrichment results
pdf('barplot_GOanalysis_10cats.pdf',h=6,w=4.25)
barplot(enrich_up, font.size = 11,showCategory=10)
dev.off()

save(sce.subset, hdg, pc.choice.hpc,marks,file="/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50_clusteredPCs.rda")
save(summed,y,x,res,enrich_go,file='/users/enelson/mouseHPC_broad_differentialExpressionAnalysis_results.rda')


