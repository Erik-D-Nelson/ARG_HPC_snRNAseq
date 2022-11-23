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
load("processed_data/sce_subset.rda")


##Differential expression analysis, broad===================================
###Set up data==================================================================
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


#grouped heatmaps for figure 2
pdf('plots/fig2/heatmap_ARGs_bySample.pdf',w=4,h=3)
features=list("rPRG"<-c('Arc' , 'Egr3', 'Gadd45b', 'Junb'),
           "dPRG"<-c('Baz1a','Brinp1','Grasp','Bdnf'),
           "SRG"<-c('Nptx2','Mfap3l','Kcna1','Kcnj2'))
for(i in 1:length(features)){
plotGroupedHeatmap(sce.subset,features=features[[i]],cluster_rows=F,cluster_cols=F,
                   group='Sample',center=T,fontsize=11)}
dev.off()

##volcano plot
data<-as.data.frame(res)
data <- data %>% 
  mutate(
    Expression = case_when(logFC >= 0 & FDR <= 0.05 ~ "Upregulated",
                           logFC <= 0 & FDR <= 0.05 ~ "Downregulated",
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
  scale_color_manual(values = c("indianred2", 
                                "gray50", 
                                "cornflowerblue")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position='right') +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_text_repel(data = top_genes,
                  mapping = aes(logFC, -log(FDR,10),
                  color=Expression,
                  label=rownames(top_genes)))

pdf('plots/fig2/volcano.pdf',h=7,w=8.5)
p2
dev.off()

#GO analysis for figure 2
x<-x[x$logFC > 0,]
enrich_go <- enrichGO(gene = x$ID,
                      OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "BP",
                      pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)


## Visualize enrichment results
pdf('plots/fig2/barplot_GOanalysis_10cats.pdf',h=6,w=4.25)
barplot(enrich_go, font.size = 11,showCategory=10)
dev.off()

save(summed,x,res,enrich_go,file='processed_data/initial_de_analysis.rda')

##Pvalue histogram
dist<-ggplot(as.data.frame(res),aes(x=PValue))+ geom_histogram(bins=30)

pdf('plots/figS3/figS3.pdf',h=15,w=15)
dist
dev.off()
