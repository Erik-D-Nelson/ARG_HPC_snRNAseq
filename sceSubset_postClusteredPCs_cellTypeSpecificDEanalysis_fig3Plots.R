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
library(UpSetR)



for(i in 1:length(x)){
  x[[i]]<-x[[i]][!is.na(x[[i]]$FDR),]
}

for(i in 1:length(x)){
  x[[i]]<-x[[i]][x[[i]]$logFC > 0,]
}

for(i in 1:length(x)){
  x[[i]]<-x[[i]][order(x[[i]]$FDR),]
}

df<-list()
for(i in 1:length(x)){
  x[[i]]<-x[[i]][rownames(x[[i]]) %in% chosen,]
  df[[i]]<-x[[i]]$FDR
  names(df[[i]])<-rownames(x[[i]])
}

for(i in 1:length(df)){
df[[i]]<-df[[i]][match(names(df[[1]]),names(df[[i]]))]
}

CA1<-intersect(rownames(x[[1]]),retrieved$SYMBOL)[1:5]


load("/users/enelson/20210907_sceSubset_postClustering_moreHDG_k20&50_clusteredPCs.rda")
load('/users/enelson/mouseHPC_broad_differentialExpressionAnalysis_results.rda')
load('argList.rda') #this was compiled from an excel sheet from a paper from Jesse Gray's lab, linked later in script

###Cell type-specific differential expression analysis==========================
summed$cellType_DE<-as.character(summed$cellType)
summed$cellType_DE<-factor(ifelse(summed$cellType_DE %in% c('DG.1','DG.2'),'DG',summed$cellType_DE))

de.dge2 <- pseudoBulkDGE(summed,
                        label=summed$cellType_DE,
                        design=~condition,
                        coef=2,
                        condition=summed$condition,
)
is.de <- decideTestsPerLabel(de.dge2, threshold=0.05)
summarizeTestsPerLabel(is.de)

a <- plot.de %>%
  as_tibble(rownames = "annotations") %>%
  gather(Gene, Member, -annotations) %>%
  filter(Member) %>%
  select(- Member)
a<-a %>%
  group_by(Gene) %>%
  summarize(annotation = list(annotations))
P<-a %>% ggplot(aes(x = annotation)) +
      geom_bar(stat='count') +
      stat_count(geom = "text", colour = "white", size = 3,
        aes(label = ..count..),position=position_stack(vjust=0.5))+
      scale_x_upset() +
      scale_y_continuous(trans='log10') +
      ylab(expression('log'[10]*'intersection size')) +
      theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(size = 14))+
      ggtitle('Upregulated DE genes')
    


names(de.dge)<-c('CA1','CA2','CA3.1','CA3.2','CA4',
                 'DG','GABA.1','GABA.2','GABA.3','GABA.4',
                 'GABA.5','PS.1','PS.2','L2/3','L5/Po','L6/6b',
                 'Sub')
is.de <- decideTestsPerLabel(de.dge, threshold=0.05)
summarizeTestsPerLabel(is.de)

BPtotal<-BP_filter[c(1,3,4,6,9,11,12)]
BPp<-list()
for(i in 1:length(BPtotal)){
  BPp[[i]]<-BPtotal[[i]]$p.adjust
  names(BPp[[i]])<-rownames(BPtotal[[i]])
}

for(i in 1:length(BPp)){
BPp[[i]]<-BPp[[i]][match(frame,names(BPp[[i]]))]
}

goData<-godata(OrgDb = 'org.Mm.eg.db',
                 keytype = "ENTREZID",
                 ont='BP',
                 keytype = "ENTREZID",
                 ont="BP",
                 computeIC = TRUE,
                 processTCSS = TRUE,
                 cutoff = 0.7
)

x<-simplify(
  go_BP[[6]],
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = goData
)


ranks<-list()
for(i in 1:ncol(frame)){
ranks[[i]]<-ifelse(!is.na(frame[,i]),rank(-frame[,i]),NA)}
y<-ifelse(!is.na(frame[,4]),frame[,4],rowMeans(y,na.rm=T))

for(i in 1:length(BPtotal)){
BPtotal[[i]]<-(BPtotal[[i]][match(rownames(BPtotal[[4]]),rownames(BPtotal[[i]])),])
}
x<-as.data.frame(summarizeTestsPerLabel(de.dge))
colnames(x)<-c('Downregulated',0,'Upregulated',NA)
x$sum<-x$Upregulated + x$Downregulated
x<-x[order(x$sum,decreasing=T),]



test = data.frame(
  group = c(rep("Upregulated", nrow(x)), rep("Downregulated", nrow(x))),
  annotation = rownames(x),
  DEGs = c(x$Upregulated, -x$Downregulated),
  stringsAsFactors = F)

test$annotation<-factor(test$annotation,levels=rownames(x))




pdf('test_barplot.pdf',h=5,w=6)
ggplot(test, aes(x = annotation, y = DEGs, fill = group))+
  geom_col()  +
  scale_y_continuous(limits=c(-5000,4630),breaks=seq(-5000,4630,by=50)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        text=element_text(size = 14),
        legend.position='bottom') +
  geom_text(aes(label = abs(DEGs)),vjust=c(rep(-.25,17),rep(1.2,17)),size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  scale_y_break(c(300,4550)) + scale_y_break(c(-50,-4950))  
dev.off()


top_genes <- bind_rows(
  data %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(FDR, desc(abs(logFC))),
  data %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(FDR, desc(abs(logFC)))
)
top_genes<-data[rownames(data) %in% features,]

data<-as.data.frame(de.dge[[6]])
data <- data %>% 
  mutate(
    Expression = case_when(logFC >= 0 & FDR <= 0.05 ~ "Upregulated",
                           logFC <= 0 & FDR <= 0.05 ~ "Downregulated",
                           TRUE ~ "Unchanged")
  )
top_genes<-data[rownames(data) %in% features,]

pdf('dg_volcano.pdf',h=2.5,w=4)
ggplot(data, aes(logFC, -log(PValue,10))) +
  geom_point(aes(color = Expression), size = 1/2) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"PValue")) +
  scale_color_manual(values = c("dodgerblue3", "black", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position='none') +
  ggtitle('DG')+
geom_label_repel(data = top_genes,
          mapping = aes(logFC, -log(PValue,10),
          label=rownames(top_genes),
          color=Expression)) + 
          theme(legend.position='none') 
dev.off()

  #  geom_bar(stat="identity", position="identity") +
  # ylim(min, max(x$goals)) +
  #scale_y_continuous(breaks = seq(min(test$y), max(test$y), 1)) +

de.specific <- pseudoBulkSpecific(summed,
                        label=summed$cellType_DE,
                        design=~condition,
                        coef=2,
                        condition=summed$condition
)
is.de.specific <- decideTestsPerLabel(de.specific, threshold=0.05)
summarizeTestsPerLabel(is.de.specific)

##Pull out DEGs for 3 labels with most DE genes

is.de<-as.data.frame(is.de)
DEG.all<-list()
de.dge<-de.dge.og
for(i in 1:length(de.dge)){
    de.dge[[i]]$F<-NULL
    de.dge[[i]]$logCPM<-NULL}
  DEG.all[[i]]<-rownames(de.dge[[i]][de.dge[[i]]$FDR<0.05,])
}


dist<-list()
for(i in 1:length(de.dge)){  
dist[[i]]<-ggplot(as.data.frame(de.dge[[i]]),aes(x=PValue))+
  geom_histogram(bins=50)+
  ggtitle(paste0('Gene PValue Histogram: ',names(de.dge)[i]))}

for(i in 1:length(de.dge)){ 
  plots[[i]]<-plot(hist[[i]])}

DEG.all[[i]]<-is.de[[i]][is.de$DG %in% c(1),]
is.de<-as.data.frame(is.de)
DEG<-list()
for(i in 1:length(colnames(is.de))){
  DEG[[i]]<-rownames(is.de)[is.de[,i] %in% c(1)]
}

is.de<-as.data.frame(is.de)
DEG.down<-list()
for(i in 1:length(colnames(is.de))){
  DEG.down[[i]]<-rownames(is.de)[is.de[,i] %in% c(-1)]
}

DEG.DG<-rownames(is.de)[is.de$DG %in% c(1)]
DEG.CA1<-rownames(is.de)[is.de$CA1 %in% c(1)]
DEG.CA3.1<-rownames(is.de)[is.de$CA3.1 %in% c(1)]
DEG.PS.1<-rownames(is.de)[is.de$ProS.1 %in% c(1)]
DEG.all<-list(DEG.DG,DEG.CA1,DEG.CA3.1,DEG.PS.1)
DEG.up<-unique(do.call(c,DEG.up))

DEG.DG.down<-rownames(is.de)[is.de$DG %in% c(-1)]
breaks2<-seq(-1.5,3,by=0.045) 

for(i in 1:length(de.dge.synapse)){
  de.dge.synapse[[i]]<-de.dge.synapse[[i]][de.dge.synapse[[i]]$FDR < 0.05,]
  de.dge.synapse[[i]]<-de.dge.synapse[[i]][rownames(de.dge.synapse[[i]]) %in% location,]}
de.dge[[i]][is.na(de.dge[[i]]$logFC),]<-0
}

map<-data.frame(row.names=rownames(de.dge[[1]]))
for(i in 1:length(de.dge)){
map[i]<-de.dge[[i]]$logFC
}
colnames(map)<-names(de.dge)



<-data.frame('CA1'=de.dge[['CA1']]$logFC,
                'CA3.1'=de.dge[['CA3.1']]$logFC,
                'DG'=de.dge[['DG']]$logFC,
                'PS.1'=de.dge[['ProS.1']]$logFC,
                row.names=rownames(de.dge[['CA1']]))
map <- map[order(map$DG,decreasing=T),]



##Pull out DEGs that are also ARGs (as defined by Jesse Gray's group's paper
##Link is here: https://www.cell.com/neuron/fulltext/S0896-6273(18)30285-X?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS089662731830285X%3Fshowall%3Dtrue)
DEG.ARG.DG<-DEG.DG[DEG.DG %in% arg.list]
DEG.ARG.CA1<-DEG.CA1[DEG.CA1 %in% arg.list]
DEG.ARG.CA3.1<-DEG.CA3.1[DEG.CA3.1 %in% arg.list]
DEG.ARG.PS.1<-DEG.PS.1[DEG.PS.1 %in% arg.list]

DEG.SRG<-list()
DEG.rPRG<-list()
DEG.dPRG<-list()
for(i in 1:length(DEG.arg.all)){
  DEG.SRG[[i]]<-intersect(DEG.arg.all[[i]],SRG)
  DEG.dPRG[[i]]<-intersect(DEG.arg.all[[i]],dPRG)
  DEG.rPRG[[i]]<-intersect(DEG.arg.all[[i]],rPRG)
}
names(DEG.SRG)<-names(DEG.arg.all)
names(DEG.dPRG)<-names(DEG.SRG)
names(DEG.rPRG)<-names(DEG.SRG)


##make lists of lists
#DEG.all<-list(DEG.DG,DEG.CA1,DEG.CA3.1,DEG.PS.1)
names(DEG.all)<-c('DG','CA1','CA3.1','PS.1')
DEG.arg.all<-list(DEG.ARG.DG,DEG.ARG.CA1,DEG.ARG.CA3.1,DEG.ARG.PS.1)

##make upset lists
#upset.all<-fromList(DEG.all)
upset.arg<-fromList(DEG.arg.all)

##upset plots for fig 3
#jpeg('upsetPlot_allDEGs.jpg')
upset(upset.all,yaxt='n')

dev.off()

jpeg('upsetPlot_ARGSonly.jpg')
upset(upset.arg)
dev.off()




##make new SCE object for CA3.D,CA1.V,and DG for the violin plots
sce.fig3<-sce.subset[,sce.subset$cellType %in% c('CA1','CA3.1','ProS.1','DG.1','DG.2','CA2')]

##set up labels for plots
sce.fig3$cellType<-droplevels(sce.fig3$cellType)
sce.fig3$annotation<-as.character(sce.fig3$cellType)
sce.fig3$annotation<-factor(ifelse(sce.fig3$condition=='control',
         paste0(sce.fig3$annotation,'.sham'),paste0(sce.fig3$annotation,'.seizure')))
         ifelse(sce.fig3$cellType=='DG','DG', 
                ifelse(sce.fig3$cellType=='DG.1','DG.seizure',sce.fig3$cellType)
)))
'Dner','Mfap3l','Arc','Fam126b','Mir670hg','Pde4b'
##Violin plots
pdf('2022_seizure_violin.pdf',w=4,h=8)
plotExpression(sce.subset,features=features,
               x="annotation", colour_by="cellType", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3)+
  theme(axis.text.x = element_text(angle = 45),
        text=element_text(size = 14))
#  ggtitle('Upregulated in CA1.2 & DG only')
dev.off()

cur.specific <- de.specific[["CA2"]]
cur.specific <- cur.specific[order(cur.specific$PValue),]
cur.specific

sce.fig3$annotation<-sce.fig3$cellType
sce.fig3$cellType<-ifelse(sce.fig3$annotation %in% c(DG))

pdf('2022_seizure_violin_all.pdf',h=2.5,w=4)
plotExpression(sce.subset,features='Mir670hg',
               x="active", colour_by="cellType", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3)+
  theme(axis.text.x = element_text(angle = 45,vjust=0.6),
        text=element_text(size = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  ggtitle('Upregulated in all')
dev.off()

pdf('20210922_seizure_violin_DGOnly.pdf',h=2.5,w=4)
plotExpression(sce.fig3,features=c("Inhba"),
               x="annotation", colour_by="cellType", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3) + 
  theme(axis.text.x = element_text(angle = 45,vjust=0.6),
                                    text=element_text(size = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  ggtitle('Upregulated in DG')
dev.off()

pdf('20210922_seizure_violin_CA_only.pdf',h=2.5,w=4)
plotExpression(sce.fig3,features=c("Trpc4"),
               x="annotation", colour_by="cellType", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3)+
  theme(axis.text.x = element_text(angle = 45,vjust=0.6),
        text=element_text(size = 10))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  ggtitle('Not upregulated in DG')
dev.off()

save(de.specific,is.de,DEG.all,DEG.arg.all,sce.fig3,file='HPC_fig3Objects.rda')

enrich<-list()
for(i in 1:length(DEG.up)){
enrich[[i]] <- enrichGO(gene = DEG.up[[i]],
                      OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "all",
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                      universe=rownames(de.dge[[i]]))
}
names(enrich)<-names(DEG.up)

BP_filter<-list()
for(i in 1:length(go_BP)){
BP_filter[[i]]<-go_BP[[i]][go_BP[[i]]$p.adjust < 0.001,]
}

enrich_CA3 <- enrichGO(gene = DEG.CA3.1,
                          OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "all",
                          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

enrich_PS.1 <- enrichGO(gene = DEG.PS.1,
                          OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "all",
                          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

enrich_DG <- enrichGO(gene = DEG.DG,
                         OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "all",
                         pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

enrich_DG_CC <- enrichGO(gene = DEG.DG,
                          OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "CC",
                          pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)

enrich_PS_BP <- enrichGO(gene = DEG.PS.1,
                         OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                         pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)

enrich_PS_CC <- enrichGO(gene = DEG.PS.1,
                         OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "CC",
                         pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)

t<-list()
for(i in 1:length(de.dge)){
  t[[i]]<-sign(de.dge[[i]]$logFC) * sqrt(de.dge[[i]]$F)
  names(t[[i]])<-rownames(de.dge[[i]])
  t[[i]]<-t[[i]][!is.na(t[[i]])]
  t[[i]]<-t[[i]][order(t[[i]],decreasing=T)]

}

for(i in 1:length(pvals.test)){
  for(x in 1:length(pvals.test[[i]])){
    names(pvals.test[[i]][[x]])<-enrichList[[i]][[x]]$Description
  }
}

for(i in 1:length(pvals.test)){
  for(x in 1:length(pvals.test[[i]])){
  pvals.test[[i]][[x]]<-pvals.test[[i]][[x]][match(clusterTerms[[i]],
                                                   names(pvals.test[[i]][[x]]))]
  }
}

for(i in 1:length(pvals.test)){
  for(x in 1:length(pvals.test[[i]])){
    names(pvals.test[[i]][[x]])<-clusterTerms[[i]]
  }
}

map<-list()
for(i in 1:length(pvals.test)){
  map[[i]]<-as.data.frame(do.call(cbind,pvals.test[[i]]))
  map[[i]][is.na(map[[i]])]<-1
  map[[i]] = -log10(map[[i]])
}

pdf('go_cluster_heatmaps.pdf',w=16,h=8)
for(i in 1:length(map)){
  ifelse(nrow(map[[i]])>1,pheatmap(map[[i]]),
         pheatmap(map[[i]],cluster_rows=F))}

for(i in 1:length(map)){
  for(x in 1:length(map[[i]])){
    colnames(map[[i]][[x]])<-names(de.dge
  }
}

ggplot(test, aes(x = annotation, y = DEGs, fill = group))+
  geom_col()  +
  scale_y_continuous(limits=c(-5000,4630),breaks=seq(-5000,4630,by=50)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        text=element_text(size = 14),
        legend.position='bottom') +
  geom_text(aes(label = abs(DEGs)),vjust=c(rep(-.25,17),rep(1.2,17)),size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  scale_y_break(c(300,4550)) + scale_y_break(c(-50,-4950))  


for(i in 1:length(map)){
colnames(map[[i]])[c(7:17)]<-c('GABA.1','GABA.2','GABA.3','GABA.4','GABA.5',
                               'PS.1','PS.2','L2/3','L5/Po','L6/6b','Sub')}


x$sum<-x$Upregulated + x$Downregulated
x<-x[order(x$sum,decreasing=T),]

go_overview = data.frame(
  'GO' = c(rep("BP", 17), rep("CC", 17), rep("MF", 17)),
  'annotation' = rep(names(bp_lengths),3),
  'enriched GO terms' = c(bp_lengths, cc_lengths, mf_lengths),
  stringsAsFactors = F)

pdf('go_barplot.pdf',h=5,w=6)
ggplot(go_overview, aes(x = annotation, y = enriched.GO.terms, fill = GO))+
  geom_col(position='stack')  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        text=element_text(size = 14),
          legend.position = c(.95, .95),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6)) +
#  geom_text(aes(label = abs(DEGs)),vjust=c(rep(-.25,17),rep(1.2,17)),size=3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
dev.off()

tiff('simMatrix.tiff',h=5,w=6)
pheatmap(simMatrix,
            annotation_row=redTerms,
            fontsize=8,
            show_rownames=F,show_colnames=F,
            legend=T,
            annotation_legend=T,
            annotation_names_row=F,
            main='BP Term Similarity')
dev.off()


clust1<-c(
  'modulation of chemical synaptic transmission',
  "regulation of trans-synaptic signaling",
  'regulation of synaptic plasticity',
  'positive regulation of synaptic transmission',
  "G protein-coupled receptor signaling pathway")

clust5<-c('neurotransmitter transport',
          'vesicle−mediated transport in synapse',
          'synaptic vesicle cycle',
          'signal release',
          'neurotransmitter secretion')

clust13<-c('synapse organization',
           'dendrite development',
           'regulation of neuron projection development',
           'synapse assembly',
           'cell junction assembly')

clust4<-c('cellular response to growth factor stimulus',
          'response to organic cyclic compound',
          'transforming growth factor beta receptor signaling pathway',
          'response to hormone',
          'response to cytokine')