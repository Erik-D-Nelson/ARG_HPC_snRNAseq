library(SingleCellExperiment)
library(EnsDb.Mmusculus.v79)
library(scater)
library(scran)
library(ggplot2)
library(dplyr)
library(magrittr)
library(bluster)
library(pheatmap)
library(edgeR)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggupset)
library(edgeR)


###Cell type-specific differential expression analysis==========================
##Set up data
sce.subset$cellType_DE<-as.character(sce.subset$cellType)
sce.subset$cellType_DE<-factor(ifelse(sce.subset$cellType_DE %in% c('GC.1','GC.2'),'GC',sce.subset$cellType_DE))
summed <- aggregateAcrossCells(sce.subset,
                               ids=colData(sce.subset)
                               [c('Sample','cellType_DE')])
#do DE analysis
de.dge <- pseudoBulkDGE(summed,
                         label=summed$cellType_DE,
                         design=~condition,
                         coef=2,
                         condition=summed$condition,
)
is.de <- decideTestsPerLabel(de.dge, threshold=0.05)
summarizeTestsPerLabel(is.de)

##Set up data for double-sided barplot for fig 1
x<-as.data.frame(summarizeTestsPerLabel(de.dge))
colnames(x)<-c('Downregulated',0,'Upregulated',NA)
x$sum<-x$Upregulated + x$Downregulated
x<-x[order(x$sum,decreasing=F),]

bar_df = data.frame(
  group = c(rep("Upregulated", nrow(x)), rep("Downregulated", nrow(x))),
  annotation = rownames(x),
  DEGs = c(x$Upregulated, -x$Downregulated),
  stringsAsFactors = F)

bar_df$annotation<-factor(bar_df$annotation,levels=rownames(x))

##make barplot
pdf('DEG_barplot.pdf',h=3.5,w=6)
ggplot(bar_df, aes(x = annotation, y = DEGs, fill = group))+
  geom_col()  +
  scale_fill_manual(values=c('cornflowerblue','indianred2'))+
  scale_y_continuous(limits=c(-5000,4630),
                     breaks=c(-10000,-1000,-100,-10,0,10,100,1000,10000),
                     trans='pseudo_log') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        legend.position='right') +
  xlab('Annotated cell type') + ylab('DEGs (log10 scale)')+
  geom_text(aes(label = abs(DEGs)),vjust=c(rep(-.25,17),rep(1.2,17)),size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

##set up data for upset barplot
plot.de<-as.data.frame(is.de)
plot.de<-plot.de[,c(1,3,11,15)]
for(i in 1:ncol(plot.de)){
 plot.de[,i]<-ifelse(plot.de[,i]==1,T,F)
}

a <- plot.de %>%
  as_tibble(rownames = "gene") %>%
  gather(annotation, Member, -gene) %>%
  filter(Member) %>%
  select(- Member)
a<-a %>%
  group_by(gene) %>%
  summarize(annotation = list(annotation))

##make upset barplot
pdf('fig3_upset_plot.pdf',h=3.5,w=3.6)
a %>% ggplot(aes(x = annotation)) +
      geom_bar(stat='count',fill='cornflowerblue') +
      stat_count(geom = "text", colour = "black", size = 3,
        aes(label = ..count..),position=position_stack(vjust=0.5))+
      scale_x_upset(reverse=T) +
      scale_y_continuous(trans='log10') +
      ylab('Upregulated DEGs (log10 scale)') + xlab('Annotated Cell Type')+
      theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(colour='black'),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 13))#+
      #ggtitle('Upregulated DE genes')
dev.off()

##volcano plots

####violin plots
##make new SCE object for CA3.1,CA1,PS.1, and GC for the violin plots
sce.fig3<-sce.subset[,sce.subset$cellType %in% c('CA1','CA3.1','PS.1','GC.1','GC.2')]

##set up labels for plots
sce.fig3$cellType<-droplevels(sce.fig3$cellType)
sce.fig3$annotation<-as.character(sce.fig3$cellType)
sce.fig3$annotation<-factor(ifelse(sce.fig3$condition=='control',
                                   paste0(sce.fig3$annotation,'.sham'),paste0(sce.fig3$annotation,'.seizure')))
ifelse(sce.fig3$cellType=='DG','DG', 
       ifelse(sce.fig3$cellType=='DG.1','DG.seizure',sce.fig3$cellType)
)))



names(de.dge)<-c('CA1','CA2','CA3.1','CA3.2','CA4',
                 'DG','GABA.1','GABA.2','GABA.3','GABA.4',
                 'GABA.5','PS.1','PS.2','L2/3','L5/Po','L6/6b',
                 'Sub')
is.de <- decideTestsPerLabel(de.dge, threshold=0.05)
summarizeTestsPerLabel(is.de)



goData<-godata(OrgDb = 'org.Mm.eg.db',
                 keytype = "ENTREZID",
                 ont='BP',
                 keytype = "ENTREZID",
                 ont="BP",
                 computeIC = TRUE,
                 processTCSS = TRUE,
                 cutoff = 0.7
)



ranks<-list()
for(i in 1:ncol(frame)){
ranks[[i]]<-ifelse(!is.na(frame[,i]),rank(-frame[,i]),NA)}
y<-ifelse(!is.na(frame[,4]),frame[,4],rowMeans(y,na.rm=T))

for(i in 1:length(BPtotal)){
BPtotal[[i]]<-(BPtotal[[i]][match(rownames(BPtotal[[4]]),rownames(BPtotal[[i]])),])
}



test = data.frame(
  group = c(rep("Upregulated", nrow(x)), rep("Downregulated", nrow(x))),
  annotation = rownames(x),
  DEGs = c(x$Upregulated, -x$Downregulated),
  stringsAsFactors = F)

test$annotation<-factor(test$annotation,levels=rownames(x))





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
  scale_color_manual(values = c("cornflowerblue", "gray50", "indianred")) +
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



##Pull out DEGs for 3 labels with most DE genes

is.de<-as.data.frame(is.de)
DEG.all<-list()
de.dge.og<-de.dge
for(i in 1:length(de.dge)){
    de.dge[[i]]$F<-NULL
    de.dge[[i]]$logCPM<-NULL
    de.dge[[i]]$FDR[is.na(de.dge[[i]]$FDR)]<-1
  DEG.all[[i]]<-rownames(de.dge[[i]][de.dge[[i]]$FDR<0.05,])}


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
  names(DEG[[i]]<-rownames(is.de)[is.de[,i] %in% c(1)])
}

is.de<-as.data.frame(is.de)
DEG.down<-list()
for(i in 1:length(colnames(is.de))){
  DEG.down[[i]]<-rownames(is.de)[is.de[,i] %in% c(-1)]
}

<-data.frame('CA1'=de.dge[['CA1']]$logFC,
                'CA3.1'=de.dge[['CA3.1']]$logFC,
                'DG'=de.dge[['DG']]$logFC,
                'PS.1'=de.dge[['ProS.1']]$logFC,
                row.names=rownames(de.dge[['CA1']]))
map <- map[order(map$DG,decreasing=T),]


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
for(i in 1:length(bp)){
  map<-as.data.frame(do.call(cbind,bp))
map[is.na(map)]<-1
map = -log10(map)
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

clust2<-c('transmembrane receptor protein tyrosine kinase signaling pathway',
          'MAPK cascade',
          'regulation of MAPK cascade',
          'small GTPase mediated signal transduction',
          'ERK1 and ERK2 cascade')


plotExpression(sce.fig3,features=c('Mir670hg','Fam126b','Ppm1h',
                                  'Tll1', 'Baz1a','Inhba',
                                  'Kdm2b','Thsd4',''),
               x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3) + 
  theme(axis.text.x = element_text(angle = 45,vjust=0.6),
        text=element_text(size = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(co  theme(axis.text.x = element_text(angle = 45,vjust=0.6),
lour = "black"))

