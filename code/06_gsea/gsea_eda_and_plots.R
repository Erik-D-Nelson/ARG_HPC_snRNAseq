##New GSEA figs!
library(scater)
library(scran)
library(SingleCellExperiment)
library(clusterProfiler)
library(org.Mm.eg.db)
library(rrvgo)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(GOSemSim)

##We need to build some dotplots for different GO clusters
load('processed_data/specific_gsea_results.rda')
load('processed_data/de_dge.rda')

##get Count for each list element

bp_Count<-list()
for(i in 1:length(go_BP)){
  ifelse(nrow(go_BP[[i]]) > 0,
    bp_Count[[i]]<-go_BP[[i]]$core_enrichment,0)
}

for(i in 1:length(bp_Count)){    
  for(x in 1:length(bp_Count[[i]])){
    if(i %in% c(1:15,17)){
      bp_Count[[i]][x]<-strsplit(bp_Count[[i]][[x]],'/')
    }
  }
}

bp<-lapply(go_BP,as.data.frame)
for(i in c(1:15,17)){
  bp[[i]]$Count<-lengths(bp_Count[[i]])
}

mf_Count<-list()
for(i in 1:length(go_MF)){
  ifelse(nrow(go_MF[[i]]) > 0,
         mf_Count[[i]]<-go_MF[[i]]$core_enrichment,0)
}

for(i in 1:length(mf_Count)){
  if(length(mf_Count[[i]])>0){
    for(x in 1:length(mf_Count[[i]])){
      mf_Count[[i]][x]<-strsplit(mf_Count[[i]][[x]],'/')
    }
  }
}

mf<-lapply(go_MF,as.data.frame)
for(i in 1:length(mf)){
if(length(mf_Count[[i]])>0){
  mf[[i]]$Count<-lengths(mf_Count[[i]])
}
}
cc_Count<-list()
for(i in 1:length(go_CC)){
  ifelse(nrow(go_CC[[i]]) > 0,
         cc_Count[[i]]<-go_CC[[i]]$core_enrichment,0)
}

for(i in 1:length(cc_Count)){
  if(length(cc_Count[[i]])>0){
    for(x in 1:length(cc_Count[[i]])){
      cc_Count[[i]][x]<-strsplit(cc_Count[[i]][[x]],'/')
    }
  }
}

cc<-lapply(go_CC,as.data.frame)
for(i in 1:length(cc)){
if(length(cc_Count[[i]])>0){
  cc[[i]]$Count<-lengths(cc_Count[[i]])
}
}

##now we have Counts!!! lets get all of the GO terms across all cell types for clustering
bp_total<-list()
for(i in 1:length(bp)){
  bp_total[[i]]<-rownames(bp[[i]])
}

bp_total<-unique(do.call(c, bp_total))

for(i in 1:length(bp)){
  bp[[i]]<-bp[[i]][match(bp_total,rownames(bp[[i]])),]
}

cc_total<-list()
for(i in 1:length(cc)){
  cc_total[[i]]<-rownames(cc[[i]])
}

cc_total<-unique(do.call(c, cc_total))

for(i in 1:length(cc)){
  cc[[i]]<-cc[[i]][match(cc_total,rownames(cc[[i]])),]
}

mf_total<-list()
for(i in 1:length(mf)){
  mf_total[[i]]<-rownames(mf[[i]])
}

mf_total<-unique(do.call(c, mf_total))

for(i in 1:length(mf)){
  mf[[i]]<-mf[[i]][match(mf_total,rownames(mf[[i]])),]
}

##now we need to get a named vector of average fdr-adjusted p-values for naming clusters in rrvgo
scores_fun<-function(x,y){
  for(i in 1:length(x)){
    x[[i]]<-x[[i]][match(y,rownames(x[[i]])),]
  }
  
  fdr<-list()
  for(i in 1:length(x)){
    fdr[[i]]<-x[[i]]$p.adjust
  }
  fdr_df<-do.call(cbind.data.frame, fdr)
  fdr_df[is.na(fdr_df)]<-1
  fdr_df= -log10(fdr_df)
  fdr_df$mean<-rowMeans(fdr_df)
  scores<-fdr_df$mean
  names(scores)<-y
  return(scores)
}

results<-list(bp,cc,mf)
terms<-list(bp_total,cc_total,mf_total)

scores<-list()
for(i in 1:length(results)){
  scores[[i]]<-scores_fun(results[[i]],terms[[i]])
}
# ##we need to cluster all terms from all ontologies


bp_data<-godata(
  OrgDb = 'org.Mm.eg.db',
  keytype = "ENTREZID",
  'BP',
  computeIC = TRUE
)

cc_data<-godata(
  OrgDb = 'org.Mm.eg.db',
  keytype = "ENTREZID",
  'CC',
  computeIC = TRUE
)

mf_data<-godata(
  OrgDb = 'org.Mm.eg.db',
  keytype = "ENTREZID",
  'MF',
  computeIC = TRUE
)

bp_simmatrix<-calculateSimMatrix(
  bp_total,
  'org.Mm.eg.db',
  semdata = bp_data,
  ont = c("BP"),
  method = "Wang"
)


#bp_scores<-bp_scores[names(bp_scores) %in% rownames(bp_simmatrix)]
bp_reducedTerms <- reduceSimMatrix(bp_simmatrix,
                                threshold=0.9,
                                scores=scores[[1]],
                                orgdb="org.Mm.eg.db")


cc_simmatrix<-calculateSimMatrix(
  cc_total,
  'org.Mm.eg.db',
  semdata = cc_data,
  ont = c("CC"),
  method = "Wang"
)


cc_reducedTerms <- reduceSimMatrix(cc_simmatrix,
                                   threshold=0.7,
                                   scores=scores[[2]],
                                   orgdb="org.Mm.eg.db")

mf_simmatrix<-calculateSimMatrix(
  mf_total,
  'org.Mm.eg.db',
  semdata = mf_data,
  ont = c("MF"),
  method = "Wang"
)


mf_reducedTerms <- reduceSimMatrix(mf_simmatrix,
                                   threshold=0.7,
                                  scores=scores[[3]],
                                   orgdb="org.Mm.eg.db")

parents_bp<-unique(bp_reducedTerms$parentTerm)
parents_cc<-unique(cc_reducedTerms$parentTerm)
parents_mf<-unique(mf_reducedTerms$parentTerm)




##Now let's figure out which of these clusters have the highest -log10FDR across all cell types
##This will help us decide what to plot
group_bp<-list()
for(i in 1:length(parents_bp)){
  group_bp[[i]]<-
    rownames(bp_reducedTerms)[bp_reducedTerms$parentTerm==parents_bp[[i]]]
}
names(group_bp)<-parents_bp

group_cc<-list()
for(i in 1:length(parents_cc)){
  group_cc[[i]]<-
    rownames(cc_reducedTerms)[cc_reducedTerms$parentTerm==parents_cc[[i]]]
}
names(group_cc)<-parents_cc

group_mf<-list()
for(i in 1:length(parents_mf)){
  group_mf[[i]]<-
    rownames(mf_reducedTerms)[mf_reducedTerms$parentTerm==parents_mf[[i]]]
}
names(group_mf)<-parents_mf

##Now let's go ahead calculate gene ratios and find which clusters have the highest mean fdr
#BP:
for(i in 1:length(bp)){
  bp[[i]]$setSize[is.na(bp[[i]]$setSize)]<-0
}
bp[[16]]$Count<-rep(NA,814)
for(i in 1:length(bp)){
  bp[[i]]$Count[is.na(bp[[i]]$Count)]<-0
  bp[[i]]$geneRatio<-ifelse(bp[[i]]$setSize > 0, 
                            bp[[i]]$Count/bp[[i]]$setSize,0)
}


fdr_bp<-list()
for(i in 1:length(bp)){
  fdr_bp[[i]]<-bp[[i]]$p.adjust
}

bp_fdr_df<-do.call(cbind.data.frame, fdr_bp)
bp_fdr_df[is.na(bp_fdr_df)]<-1
bp_fdr_df= -log10(bp_fdr_df)
colnames(bp_fdr_df)<-names(de.dge)
bp_fdr_df$mean<-rowMeans(bp_fdr_df)
rownames(bp_fdr_df)<-bp_total


bp_parent_df<-list()
bp_means<-list()
for(i in 1:length(group_bp)){
  bp_parent_df[[i]]<-bp_fdr_df[rownames(bp_fdr_df) %in% group_bp[[i]],]
  bp_means[[i]]<-mean(bp_parent_df[[i]]$mean)
}

bp_median<-list()
for(i in 1:length(group_bp)){
bp_median[[i]]<-median(bp_parent_df[[i]]$mean)
}

geneRatio_bp<-list()
for(i in 1:length(bp)){
  geneRatio_bp[[i]]<-bp[[i]]$geneRatio
}

bp_geneRatio_df<-do.call(cbind.data.frame, geneRatio_bp)
bp_geneRatio_df[is.na(bp_geneRatio_df)]<-0
colnames(bp_geneRatio_df)<-names(de.dge)
rownames(bp_geneRatio_df)<-bp_total




bp_parent_geneRatio<-list()
for(i in 1:length(group_bp)){
  bp_parent_geneRatio[[i]]<-bp_geneRatio_df[rownames(bp_geneRatio_df) %in% group_bp[[i]],]
}



names(bp_parent_df)<-names(group_bp)
names(bp_parent_geneRatio)<-names(group_bp)
names(bp_means)<-names(group_bp)
names(bp_median)<-names(group_bp)

#CC:
for(i in 1:length(cc)){
  cc[[i]]$setSize[is.na(cc[[i]]$setSize)]<-0
  cc[[i]]$Count[is.na(cc[[i]]$Count)]<-0
  cc[[i]]$geneRatio<-ifelse(cc[[i]]$setSize > 0, 
                            cc[[i]]$Count/cc[[i]]$setSize,0)
  
}

fdr_cc<-list()
for(i in 1:length(cc)){
  fdr_cc[[i]]<-cc[[i]]$p.adjust
}

cc_fdr_df<-do.call(cbind.data.frame, fdr_cc)
cc_fdr_df[is.na(cc_fdr_df)]<-1
cc_fdr_df= -log10(cc_fdr_df)
colnames(cc_fdr_df)<-names(de.dge)
cc_fdr_df$mean<-rowMeans(cc_fdr_df)
rownames(cc_fdr_df)<-cc_total


cc_parent_df<-list()
cc_means<-list()
for(i in 1:length(group_cc)){
  cc_parent_df[[i]]<-cc_fdr_df[rownames(cc_fdr_df) %in% group_cc[[i]],]
  cc_means[[i]]<-mean(cc_parent_df[[i]]$mean)
}

cc_median<-list()
for(i in 1:length(group_cc)){
  cc_median[[i]]<-median(cc_parent_df[[i]]$mean)
}

geneRatio_cc<-list()
for(i in 1:length(cc)){
  geneRatio_cc[[i]]<-cc[[i]]$geneRatio
}

cc_geneRatio_df<-do.call(cbind.data.frame, geneRatio_cc)
cc_geneRatio_df[is.na(cc_geneRatio_df)]<-0
colnames(cc_geneRatio_df)<-names(de.dge)
rownames(cc_geneRatio_df)<-cc_total




cc_parent_geneRatio<-list()
for(i in 1:length(group_cc)){
  cc_parent_geneRatio[[i]]<-cc_geneRatio_df[rownames(cc_geneRatio_df) %in% group_cc[[i]],]
}



names(cc_parent_df)<-names(group_cc)
names(cc_parent_geneRatio)<-names(group_cc)
names(cc_means)<-names(group_cc)
names(cc_median)<-names(group_cc)

#MF:
for(i in 1:length(mf)){
  mf[[i]]$setSize[is.na(mf[[i]]$setSize)]<-0
}
mf[[2]]$Count<-rep(NA,97)
mf[[4]]$Count<-rep(NA,97)
mf[[10]]$Count<-rep(NA,97)
mf[[11]]$Count<-rep(NA,97)
mf[[15]]$Count<-rep(NA,97)

for(i in 1:length(mf)){
  mf[[i]]$Count[is.na(mf[[i]]$Count)]<-0
  mf[[i]]$geneRatio<-ifelse(mf[[i]]$setSize > 0, 
                            mf[[i]]$Count/mf[[i]]$setSize,0)
}



fdr_mf<-list()
for(i in 1:length(mf)){
  fdr_mf[[i]]<-mf[[i]]$p.adjust
}

mf_fdr_df<-do.call(cbind.data.frame, fdr_mf)
mf_fdr_df[is.na(mf_fdr_df)]<-1
mf_fdr_df= -log10(mf_fdr_df)
colnames(mf_fdr_df)<-names(de.dge)
mf_fdr_df$mean<-rowMeans(mf_fdr_df)
rownames(mf_fdr_df)<-mf_total


mf_parent_df<-list()
mf_means<-list()
for(i in 1:length(group_mf)){
  mf_parent_df[[i]]<-mf_fdr_df[rownames(mf_fdr_df) %in% group_mf[[i]],]
  mf_means[[i]]<-mean(mf_parent_df[[i]]$mean)
}

mf_median<-list()
for(i in 1:length(group_mf)){
  mf_median[[i]]<-median(mf_parent_df[[i]]$mean)
}

geneRatio_mf<-list()
for(i in 1:length(mf)){
  geneRatio_mf[[i]]<-mf[[i]]$geneRatio
}

mf_geneRatio_df<-do.call(cbind.data.frame, geneRatio_mf)
mf_geneRatio_df[is.na(mf_geneRatio_df)]<-0
colnames(mf_geneRatio_df)<-names(de.dge)
rownames(mf_geneRatio_df)<-mf_total




mf_parent_geneRatio<-list()
for(i in 1:length(group_mf)){
  mf_parent_geneRatio[[i]]<-mf_geneRatio_df[rownames(mf_geneRatio_df) %in% group_mf[[i]],]
}



names(mf_parent_df)<-names(group_mf)
names(mf_parent_geneRatio)<-names(group_mf)
names(mf_means)<-names(group_mf)
names(mf_median)<-names(group_mf)


##heatmaps of all clusters

bp_heats<-list()
for(i in 1:length(bp_parent_df)){
  bp_heats[[i]]<-pheatmap(bp_parent_df[[i]],cluster_rows = F,cluster_cols = F,
                          scale='row')
}

cc_heats<-list()
for(i in 1:length(cc_parent_df)){
  cc_heats[[i]]<-pheatmap(cc_parent_df[[i]],cluster_rows = F,cluster_cols = F,
                          scale='row')
}


mf_heats<-list()
for(i in 1:length(mf_parent_df)){
  mf_heats[[i]]<-pheatmap(mf_parent_df[[i]],cluster_rows = F,cluster_cols = F,
                          scale='row')
}


###OK now lets make some plots
##first, write a function to return balloon plots
balloon_fun<-function(x,y){
  ###set up terms for color/shading in balloon plot
  y$mean<-x$mean
  x<-x[order(x$mean,decreasing=T),]
  y<-y[order(y$mean,decreasing=T),]
  rownames(x)<-go2term(rownames(x))$Term[match(rownames(x),
                                               go2term(rownames(x))$go_id)]
  x<-x[,1:17]
  x<-x[,c(11,3,15,1,8,16,10,4,17,6,7,9,5,12,14,2,13)]
  
  ###set up gene ratio for size of circles in balloon plot
  rownames(y)<-go2term(rownames(y))$Term[match(rownames(y),
                                               go2term(rownames(y))$go_id)]
  y<-y[,1:17]
  x<-x[,c(11,3,15,1,8,16,10,4,17,6,7,9,5,12,14,2,13)]
  
  
  x<-x[c(1:4),]
  y<-y[c(1:4),]
  x$ID<-factor(rownames(x))
  x_melt<-melt(x)

  
  
  pp<-ggplot(plot_go, aes(x =variable, y = ID)) + 
    geom_point(aes(colour=value,size=size)) +
    theme(panel.background=element_blank(),
          panel.border = element_rect(colour = "black", 
                                      fill=NA), 
          axis.text.x = element_text(angle = 45,
                                     vjust=.5)) + 
    scale_colour_gradient(low='white',
                          high='black',
                          limits=c(0,17)) +
    scale_size_area(limits=c(0,0.5))+
    labs(size='geneRatio',
         color='-log10 adjusted p-value',
         x='Annotated Cell Type',
         y='GO Term')
  outs<-list(x,x_melt,y,pp)
  return(outs)}
##now, make plots
##BP: apoptosis, growth factor signaling
signal<-balloon_fun(x=bp_parent_df[[2]],y=bp_parent_geneRatio[[2]])
growth<-balloon_fun(x=bp_parent_df[[11]],y=bp_parent_geneRatio[[11]])

##CC: synapse, secretory vesicle
synapse<-balloon_fun(x=cc_parent_df[[1]],y=cc_parent_geneRatio[[1]])
#vesicle<-balloon_fun(x=cc_parent_df[[5]],y=cc_parent_geneRatio[[5]])

##MF: phosphorylation and ubiquitination
ubiq<-balloon_fun(x=mf_parent_df[[8]],y=mf_parent_geneRatio[[8]])
kinase<-balloon_fun(x=mf_parent_df[[2]],y=mf_parent_geneRatio[[2]])


pdf('plots/fig4/synapse.pdf',w=6.5,h=4)
synapse[[4]]+ggtitle('')
dev.off()

pdf('plots/fig4/growth.pdf',w=6.5,h=4)
growth[[4]]+ggtitle('')
dev.off()

#pdf('plots/fig4/signal.pdf',w=6.5,h=4)
#signal[[4]]+ggtitle('')
#dev.off()

pdf('plots/fig4/kinase.pdf',w=6.5,h=4)
kinase[[4]]+ggtitle('')
dev.off()

pdf('plots/fig4/ubiq.pdf',w=6.5,h=4)
ubiq[[4]]+ggtitle('')
dev.off()

save(bp_reducedTerms,cc_reducedTerms,mf_reducedTerms,
     bp_simmatrix,cc_simmatrix,mf_simmatrix, file='go_clustering_results.rda')

##Fig S6
bp_plot<-list()
bp_ratios<-list()
for(i in order(unlist(bp_means),decreasing=T)){
  bp_plot[[i]]<-colMeans(bp_parent_df[[i]])
  bp_ratios[[i]]<-colMeans(bp_parent_geneRatio[[i]])
}

bp_plot<-do.call(rbind.data.frame,bp_plot)
bp_ratios<-do.call(rbind.data.frame,bp_ratios)

colnames(bp_plot)<-c("CA1", "CA2", "CA3.1", "CA3.2", "CA4",
                     "GC", "GABA.1", "GABA.2", "GABA.3", "GABA.4",
                     "GABA.5", "PS.1", "PS.2", "L2/3", "L5/Po",
                     "L6/6b", "Sub")
bp_plot<-bp_plot[,c(1:17)]
rownames(bp_plot)<-names(bp_parent_df)[c(order(unlist(bp_means),decreasing=T))]
colnames(bp_ratios)<-colnames(bp_plot)
bp_plot<-bp_plot[,c(6,12,3,1,9,13,11,4,17,7,8,10,5,14,16,2,15)]
bp_ratios<-bp_ratios[,c(6,12,3,1,9,13,11,4,17,7,8,10,5,14,16,2,15)]

bp_plot$ID<-factor(rownames(bp_plot),levels=rownames(bp_plot))
bp_plot_melt<-melt(bp_plot)
bp_plot_melt$size<-unlist(bp_ratios)

pdf('plots/figS6/pp_bp.pdf',h=7,w=15)
pp_bp<-ggplot(bp_plot_melt, aes(x =variable, y = ID)) + 
  geom_point(aes(colour=value,size=size)) +
  theme(panel.background=element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA), 
        axis.text.x = element_text(angle = 45,
                                   vjust=.5)) + 
  scale_y_discrete(limits=rev(levels(bp_plot$ID))) +
  scale_colour_gradient(low='white',
                        high='black') +
  scale_size_area()+
  labs(size='count',
       color='mean of -log10 adjusted P values',
       bp_plot='Annotated Cell Type',
       bp_ratios='GO term cluster')
dev.off()

mf_plot<-list()
mf_ratios<-list()
for(i in order(unlist(mf_means),decreasing=T)){
  mf_plot[[i]]<-colMeans(mf_parent_df[[i]])
  mf_ratios[[i]]<-colMeans(mf_parent_geneRatio[[i]])
}

mf_plot<-do.call(rbind.data.frame,mf_plot)
mf_ratios<-do.call(rbind.data.frame,mf_ratios)

colnames(mf_plot)<-c("CA1", "CA2", "CA3.1", "CA3.2", "CA4",
                     "GC", "GABA.1", "GABA.2", "GABA.3", "GABA.4",
                     "GABA.5", "PS.1", "PS.2", "L2/3", "L5/Po",
                     "L6/6b", "Sub")
mf_plot<-mf_plot[,c(1:17)]
rownames(mf_plot)<-names(mf_parent_df)[c(order(unlist(mf_means),decreasing=T))]
colnames(mf_ratios)<-colnames(mf_plot)
mf_plot<-mf_plot[,c(6,12,3,1,9,13,11,4,17,7,8,10,5,14,16,2,15)]
mf_ratios<-mf_ratios[,c(6,12,3,1,9,13,11,4,17,7,8,10,5,14,16,2,15)]

mf_plot$ID<-factor(rownames(mf_plot),levels=rownames(mf_plot))
mf_plot_melt<-melt(mf_plot)
mf_plot_melt$size<-unlist(mf_ratios)

pdf('plots/figS6/pp_mf.pdf',h=7,w=15)
pp_mf<-ggplot(mf_plot_melt, aes(x =variable, y = ID)) + 
  geom_point(aes(colour=value,size=size)) +
  theme(panel.background=element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA), 
        axis.text.x = element_text(angle = 45,
                                   vjust=.5)) + 
  scale_y_discrete(limits=rev(levels(mf_plot$ID))) +
  scale_colour_gradient(low='white',
                        high='black') +
  scale_size_area()+
  labs(size='count',
       color='mean of -log10 adjusted P values',
       mf_plot='Annotated Cell Type',
       mf_ratios='GO term cluster')
dev.off()


cc_plot<-list()
cc_ratios<-list()
for(i in order(unlist(cc_means),decreasing=T)){
  cc_plot[[i]]<-colMeans(cc_parent_df[[i]])
  cc_ratios[[i]]<-colMeans(cc_parent_geneRatio[[i]])
}

cc_plot<-do.call(rbind.data.frame,cc_plot)
cc_ratios<-do.call(rbind.data.frame,cc_ratios)

colnames(cc_plot)<-c("CA1", "CA2", "CA3.1", "CA3.2", "CA4",
                     "GC", "GABA.1", "GABA.2", "GABA.3", "GABA.4",
                     "GABA.5", "PS.1", "PS.2", "L2/3", "L5/Po",
                     "L6/6b", "Sub")
cc_plot<-cc_plot[,c(1:17)]
rownames(cc_plot)<-names(cc_parent_df)[c(order(unlist(cc_means),decreasing=T))]
colnames(cc_ratios)<-colnames(cc_plot)
cc_plot<-cc_plot[,c(6,12,3,1,9,13,11,4,17,7,8,10,5,14,16,2,15)]
cc_ratios<-cc_ratios[,c(6,12,3,1,9,13,11,4,17,7,8,10,5,14,16,2,15)]

cc_plot$ID<-factor(rownames(cc_plot),levels=rownames(cc_plot))
cc_plot_melt<-melt(cc_plot)
cc_plot_melt$size<-unlist(cc_ratios)

pdf('plots/figS6/pp_cc.pdf',h=7,w=15)
pp_cc<-ggplot(cc_plot_melt, aes(x =variable, y = ID)) + 
  geom_point(aes(colour=value,size=size)) +
  theme(panel.background=element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA), 
        axis.text.x = element_text(angle = 45,
                                   vjust=.5)) + 
  scale_y_discrete(limits=rev(levels(cc_plot$ID))) +
  scale_colour_gradient(low='white',
                        high='black') +
  scale_size_area()+
  labs(size='count',
       color='mean of -log10 adjusted P values',
       cc_plot='Annotated Cell Type',
       cc_ratios='GO term cluster')
dev.off()

pdf('plots/figS5/mf_simmatrix.pdf',h=8,w=15)
heatmapPlot(mf_simmatrix,
            mf_reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=14,show_rownames=F,show_colnames=F)
dev.off()

pdf('plots/figS5/bp_simmatrix.pdf',h=8,w=16)
heatmapPlot(bp_simmatrix,
            bp_reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=14,show_rownames=F,show_colnames=F)
dev.off()

pdf('plots/figS5/cc_simmatrix.pdf',h=8,w=12)
heatmapPlot(cc_simmatrix,
            cc_reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=14,show_rownames=F,show_colnames=F)
dev.off()