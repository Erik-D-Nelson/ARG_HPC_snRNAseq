################################################################################
### LIBD pilot 10x snRNA-seq: HPC samples
###   **Region-specific analyses**
###   - R-batch job for detxn of optimal PC space with 'sce.subset' object
################################################################################

library(SingleCellExperiment)
library(scRNAseq)
library(scater)
library(scran)
library(jaffelab)

# ===


load("/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50.rda",verbose=T)

## ** [corrected] PCA already done (interactively, using `fastMNN`) - took top 100 PCs

# getClusteredPCs() to identify working PC space
pc.choice.hpc <- getClusteredPCs(reducedDim(sce.subset, "PCA"))

# How many PCs should use in this space?
metadata(pc.choice.hpc)$chosen


## Plot n Clusters vs. d PCs
pdf("20210913_sceSubset_getClusteredPCs.pdf")
plot(pc.choice.hpc$n.pcs, pc.choice.hpc$n.clusters,
     main="PC Choice via getClusteredPCs()",
     abline(v=metadata(pc.choice.hpc)$chosen, col="red", lty="dashed", lwd=0.8)) 
dev.off()


# Save
save(sce.subset, hdg, pc.choice.hpc,
     file="/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50_clusteredPCs.rda")

rm(list=ls())
sessionInfo()