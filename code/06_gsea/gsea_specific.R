library(scater)
library(scran)
library(SingleCellExperiment)
library(edgeR)
library(org.Mm.eg.db)
library(clusterProfiler)
library(AnnotationDbi)

###First, load the cell type specific stats and get t-stats
load('processed_data/specific_de_results.rda')

for(i in 1:length(de.dge)){
de.dge[[i]]$entrez = mapIds(org.Mm.eg.db,
                    keys=rownames(de.dge[[i]]), 
                    column="ENTREZID",
                    keytype="SYMBOL",
                    multiVals="first")
}
#now, retain only those genes with entrez IDs

t<-list()
for(i in 1:length(de.dge)){
de.dge[[i]]<-de.dge[[i]][!is.na(de.dge[[i]]$entrez),]
de.dge[[i]]<-de.dge[[i]][!is.na(de.dge[[i]]$logFC),]
t[[i]]<-sign(de.dge[[i]]$logFC) * sqrt(de.dge[[i]]$F)
names(t[[i]])<-de.dge[[i]]$entrez
t[[i]]<-t[[i]][order(t[[i]],decreasing=T)]
}
names(t)<-names(de.dge)
###We have our ordered gene list now! Let's start by running gsea for GO terms
#To start, let's do gene sets size 25-500, but we can increase minSetSize later
#if needed
go_CC<-list()
for(i in 1:length(t)){
go_CC <- gseGO(geneList   = t[[i]],
               OrgDb        = org.Mm.eg.db,
               ont          = "CC",
               minGSSize    = 25,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               eps          = 0,
               verbose      = FALSE)
}
print(go_CC[[3]])

go_BP<-list()
for(i in 1:length(t)){
  go_BP <- gseGO(geneList   = t[[i]],
                 OrgDb        = org.Mm.eg.db,
                 ont          = "BP",
                 minGSSize    = 25,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 eps          = 0,
                 verbose      = FALSE)
}
print(go_BP[[3]])

go_MF<-list()
for(i in 1:length(t)){
  go_MF <- gseGO(geneList   = t[[i]],
                 OrgDb        = org.Mm.eg.db,
                 ont          = "MF",
                 minGSSize    = 25,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 eps          = 0,
                 verbose      = FALSE)
}
print(go_MF[[3]])

##We now have GSEA results for all 3 GO ontology groups
#Let's save and start on exploring the results and making plots.
save(go_CC,go_BP,go_MF,t,file='processed_data/specific_gsea_results.rda')
##
