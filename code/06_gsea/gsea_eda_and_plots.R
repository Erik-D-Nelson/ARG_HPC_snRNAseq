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
load('specific_gsea_results.rda')
load('de.dge.rda')

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

cc_total<-list()
for(i in 1:length(cc)){
  cc_total[[i]]<-rownames(cc[[i]])
}

cc_total<-unique(do.call(c, cc_total))

mf_total<-list()
for(i in 1:length(mf)){
  mf_total[[i]]<-rownames(mf[[i]])
}

mf_total<-unique(do.call(c, mf_total))

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

##we need to cluster all terms from all ontologies


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
                                   threshold=0.8,
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

##Now let's go ahead and find which clusters have the highest mean fdr
#BP:
for(i in 1:length(bp)){
  bp[[i]]<-bp[[i]][match(bp_total,rownames(bp[[i]])),]
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

count_bp<-list()
for(i in 1:length(bp)){
  count_bp[[i]]<-bp[[i]]$Count
}
count_bp[[16]]<-rep(NA,814)

bp_count_df<-do.call(cbind.data.frame, count_bp)
bp_count_df[is.na(bp_count_df)]<-0
colnames(bp_count_df)<-names(de.dge)
rownames(bp_count_df)<-bp_total


bp_parent_count<-list()
for(i in 1:length(group_bp)){
  bp_parent_count[[i]]<-bp_count_df[rownames(bp_count_df) %in% group_bp[[i]],]
}



names(bp_parent_df)<-names(group_bp)
names(bp_parent_count)<-names(group_bp)
names(bp_means)<-names(group_bp)
names(bp_median)<-names(group_bp)

#CC:
for(i in 1:length(cc)){
  cc[[i]]<-cc[[i]][match(cc_total,rownames(cc[[i]])),]
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

count_cc<-list()
for(i in 1:length(cc)){
  count_cc[[i]]<-cc[[i]]$Count
}

cc_count_df<-do.call(cbind.data.frame, count_cc)
cc_count_df[is.na(cc_count_df)]<-0
colnames(cc_count_df)<-names(de.dge)
rownames(cc_count_df)<-cc_total


cc_parent_count<-list()
for(i in 1:length(group_cc)){
  cc_parent_count[[i]]<-cc_count_df[rownames(cc_count_df) %in% group_cc[[i]],]
}

names(cc_parent_df)<-names(group_cc)
names(cc_parent_count)<-names(group_cc)
names(cc_means)<-names(group_cc)
names(cc_median)<-names(group_cc)

#MF:
for(i in 1:length(mf)){
  mf[[i]]<-mf[[i]][match(mf_total,rownames(mf[[i]])),]
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

count_mf<-list()
for(i in 1:length(mf)){
  count_mf[[i]]<-mf[[i]]$Count
}
count_mf[[2]]<-rep(NA,97)
count_mf[[4]]<-rep(NA,97)
count_mf[[10]]<-rep(NA,97)
count_mf[[11]]<-rep(NA,97)
count_mf[[15]]<-rep(NA,97)

mf_count_df<-do.call(cbind.data.frame, count_mf)
mf_count_df[is.na(mf_count_df)]<-0
colnames(mf_count_df)<-names(de.dge)
rownames(mf_count_df)<-mf_total


mf_parent_count<-list()
for(i in 1:length(group_mf)){
  mf_parent_count[[i]]<-mf_count_df[rownames(mf_count_df) %in% group_mf[[i]],]
}



names(mf_parent_df)<-names(group_mf)
names(mf_parent_count)<-names(group_mf)
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


###OK now lets make some plotskis
##CC we'll plot synapse and secretory vesicle
synapse<-cc_parent_df[[1]]
synapse<-synapse[order(synapse$mean,decreasing=T),]
rownames(synapse)<-go2term(rownames(synapse))$Term
synapse<-synapse[,1:17]
colnames(synapse)<-c("CA1", "CA2", "CA3.1", "CA3.2", "CA4",
  "DG", "GABA.1", "GABA.2", "GABA.3", "GABA.4",
  "GABA.5", "PS.1", "PS.2", "L2/3", "L5/Po",
  "L6/6b", "Sub")
synapse<-synapse[,c(6,3,12,1,9,13,11,4,17,7,8,10,5,14,16,2,15)]

synapse_count<-cc_parent_count[[1]]
rownames(synapse_count)<-go2term(rownames(synapse_count))$Term
synapse_count<-synapse_count[match(rownames(synapse),rownames(synapse_count)),]
synapse_count<-synapse_count[,1:17]
colnames(synapse_count)<-c("CA1", "CA2", "CA3.1", "CA3.2", "CA4",
                           "DG", "GABA.1", "GABA.2", "GABA.3", "GABA.4",
                           "GABA.5", "PS.1", "PS.2", "L2/3", "L5/Po",
                           "L6/6b", "Sub")
synapse_count<-synapse_count[,c(6,3,12,1,9,13,11,4,17,7,8,10,5,14,16,2,15)]


synapse<-synapse[c(1:5),]
synapse_count<-synapse_count[c(1:5),]
synapse$ID<-factor(rownames(synapse))
synapse_melt<-melt(synapse)


ggplot(synapse_melt, aes(x =variable, y = ID)) + 
  geom_point(aes(colour=value,size=unlist(synapse_count))) +
  theme(panel.background=element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA), 
        axis.text.x = element_text(angle = 45,
                                   vjust=.5)) + 
  scale_colour_gradient(low='white',
                        high='black') +
  scale_size_area(breaks=c(0,50,100,150,200))+
  labs(size='Count',
       color='-log10 adjusted P value',
       x='Annotation',
       y='GO Term')
dev.off()


ggplot(ion_melt, aes(x =variable, y = ID)) + 
  geom_point(aes(colour=value,size=unlist(key_ion))) +
  theme(panel.background=element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA), 
        axis.text.x = element_text(angle = 45,
                                   vjust=.5),
        text = element_text(size = 12)) + 
  scale_colour_gradient(low='white',
                        high='black') +
  scale_size_area(breaks=c(0,50,100,150,200))+
  labs(size='Count',
       color='-log10 adjusted P value',
       x='Annotation',
       y='GO Term') + 
  ggtitle('Transcription factor binding')



  reducedTerms_CC$term[reducedTerms_cc$parentTerm=='glutamatergic synapse']



synapse<-reducedTerms_CC$term[reducedTerms_CC$parentTerm=='glutamatergic synapse']

dna<-c('DNA-binding transcription activator activity',
       'transcription coregulator activity',
       'transcription corepressor activity',
       'transcription coactivator activity')

DNA-binding transcription activator activity, RNA polymerase II-specific

signal<-bp_reducedTerms$term[bp_reducedTerms$parentTerm=="G protein-coupled receptor signaling pathway"]
  
  c('modulation of chemical synaptic transmission',
          'transmembrane receptor protein tyrosine kinase signaling pathway',
          'regulation of synaptic plasticity',
          "G protein-coupled receptor signaling pathway")

tgfb<-c('cellular response to chemical stress',
        'cellular response to oxidative stress',
        'response to cytokine',
        'transforming growth factor beta receptor signaling pathway')
cc_plot<-list()
for(i in 1:length(cc)){

  
  cc_plot[[i]]<-cc[[i]][cc[[i]]$Description %in% c(rownames(synapse)),]
}



cc_plot_final<-cc_plot
for(i in 1:length(cc_plot_final)){
  cc_plot_final[[i]]<-cc_plot_final[[i]][match(synapse,
                                               cc_plot_final[[i]]$Description),]
  cc_plot_final[[i]]$Description<-factor(cc_plot_final[[i]]$Description)
}


key<-lapply(cc_plot_final,pull,Count)
keydf<-data.frame('ID'=synapse)
keydf<-key
keydf[is.na(keydf)]<-0
keydf$ID<-NULL
key_signal<-keydf

p<-lapply(bp_plot_final,pull,p.adjust)
for(i in 1:length(p)){
  p[[i]][is.na(p[[i]])]<-1
}

bp_signal<-data.frame('ID'=signal)

bp_signal[,2:18]<-p
colnames(bp_signal)[2:18]<-names(de.dge)



colnames(bp_signal)[2:18]<-recode(colnames(bp_signal)[2:18], 'ProS.1'='PS.1','ProS.2'='PS.2',
                                  'RHP.L2/3'='L2/3','RHP.L5/Po'='L5/Po',
                                  'RHP.L6/6b'='L6/6b','GABA.CGE.1'='GABA.1',
                                  'GABA.CGE.2'='GABA.2','GABA.MGE.1'='GABA.3',
                                  'GABA.MGE.2'='GABA.4','GABA.MGE.3'='GABA.5')
bp_signal<-bp_signal[,c(1,7,2:6,13,14,18,15,16,17,8:12)]
key_signal<-key_signal[,c(6,1:5,12,13,17,14,15,16,7:11)]


bp_signal[,2:18]=-log10(bp_signal[,2:18])
signal_melt<-melt(bp_signal)
bp_signal<-recode(bp_signal,variable=Annotation)


pdf('signaling.pdf',h=6,w=8)
ggplot(signal_melt, aes(x =variable, y = ID)) + 
  geom_point(aes(colour=value,size=unlist(value))) +
  theme(panel.background=element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA), 
        axis.text.x = element_text(angle = 45,
                                   vjust=.5)) + 
  scale_colour_gradient(low='white',
                        high='red') +
  labs(size='Count',
       color='-log10 adjusted P value')
dev.off()

