#we want to do 3 things with this script
##1) Get list of gene types retained when subsetting for entrez ID
###     -most of these should be noncoding guys
##2) Set up ranked gene lists for each pattern
##3) Run GSEA for all 3 GO ontology types

##first, get entrez genes
entrez = mapIds(org.Mm.eg.db,
        keys=rownames(sce.subset), 
        column="ENTREZID",
        keytype="SYMBOL",
        multiVals="first")

entrezSymbol= mapIds(org.Mm.eg.db,
                     keys=rowData(sce.subset)$Symbol, 
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")

entrezEnsembl= mapIds(org.Mm.eg.db,
                     keys=rowData(sce.subset)$ID, 
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")


entrezCombined<-entrez
entrezCombined<-ifelse(is.na(entrezCombined) & 
                               !is.na(entrezSymbol),
                       entrezSymbol,entrezCombined)
entrezCombined<-ifelse(is.na(entrezCombined) & 
                         !is.na(entrezEnsembl),
                       entrezEnsembl,entrezCombined)
                             

##print table showing how many NAs
table(is.na(entrez))
#FALSE  TRUE
#17812  2017

##What kinds of genes are these?



#now, retain only those genes with entrez IDs
sce.subset<-sce.subset[!is.na(rowData(sce.subset)$Entrez),]
index<-data.frame('Symbol'=rownames(sce.subset),
                  'Entrez'=rowData(sce.subset)$Entrez)
scores<-scores[rownames(scores) %in% index$Symbol,]
rownames(scores)<-index$Entrez
t<-list()
for(i in 1:ncol(scores)){
  t[[i]]<-scores[,i]
  t[[i]]<-t[[i]][order(t[[i]],decreasing=T)]
}
names(t)<-colnames(scores)
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
