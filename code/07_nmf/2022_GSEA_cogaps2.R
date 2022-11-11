library(clusterProfiler)
library(org.Mm.eg.db)
load("/users/enelson/ranked_genes_cogaps.rda")
###We have our ordered gene list now! Let's start by running gsea for GO terms
#To start, let's do gene sets size 25-500, but we can increase minSetSize later
#if needed
##go_CC<-list()
##for(i in 1:length(t)){
##  go_CC[[i]] <- gseGO(geneList   = t[[i]],
##                 OrgDb        = org.Mm.eg.db,
 ##                ont          = "CC",
 ##                minGSSize    = 50,
 ##                maxGSSize    = 500,
 ##                pvalueCutoff = 0.05,
 ##                eps          = 0,
  ##               nPermSimple  = 5000,
  ##               verbose      = FALSE)
##}
##print(go_CC[[55]])
##save(go_CC,file='CC_cogaps.rda')
##rm(go_CC)

go_BP<-list()
print("Running BP GO Analysis...")
for(i in 1:length(t)){
  go_BP[[i]] <- gseGO(geneList   = t[[i]],
                 OrgDb        = org.Mm.eg.db,
                 ont          = "BP",
                 minGSSize    = 50,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 eps          = 0,
                 nPermSimple  = 2500,
                 verbose      = FALSE)
}
print(go_BP[[55]])
print(go_BP[[1]])
save(go_BP,file='BP_cogaps.rda')
rm(go_BP)

go_MF<-list()
for(i in 1:length(t)){
  go_MF[[i]] <- gseGO(geneList   = t[[i]],
                 OrgDb        = org.Mm.eg.db,
                 ont          = "MF",
                 minGSSize    = 50,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 eps          = 0,
                 nPermSimple  = 2500,
                 verbose      = FALSE)
}
print(go_MF[[55]])
save(go_MF,file='MF_cogaps.rda')
rm(go_MF)


