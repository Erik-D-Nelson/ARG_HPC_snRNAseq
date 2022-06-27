library(SingleCellExperiment)
library(scater)
library(scran)
library(schex)
load(file='sceSubset.rda')

##Make new sce with hexbins using make_hexbin
hex <- make_hexbin(sce.subset, nbins = 100,
                        dimension_reduction = "UMAP", use_dims=c(1,2))

##Make cluster labels for plot
label_df <- make_hexbin_label(hex, col="cellType")

#make plot 
pp <- plot_hexbin_meta(hex, col="Pattern_20", action="median",
                       xlab='UMAP1',ylab='UMAP2') + ggtitle("Pattern 20") 
pdf("UMAP_hex_condition.pdf",h=5,w=6)
pp + ggrepel::geom_label_repel(data = label_df, aes(x=x, y=y,label=label), 
                               colour="black",  label.size = NA, fill = NA)+
  theme(legend.position='right')#,text = element_text(size = 20))
dev.off()



##Plot a feature
gene_id <-"Nr2f2"
gg <- schex::plot_hexbin_feature(hex, type="logcounts", feature=gene_id, 
                                 action="median", xlab="UMAP1", ylab="UMAP2", 
                                 title=paste0("Median of ", gene_id))
#gg + theme_void()


#sce.DG$annotation<-ifelse(sce.DG$cluster %in% c(2,8,9),'DG.sham.1',
 #                         ifelse(sce.DG$cluster==1,'DG.seizure.1',
  #                               ifelse(sce.DG$cluster==5,'DG.sham.2',
   #                                     ifelse(sce.DG$cluster==7,'DG.seizure.2',
    #                                           ifelse(sce.DG$cluster %in% c(3,6),'DG.sham.3',
     #                                                 'DG.seizure.3')))))