library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(dplyr)
library(magrittr)
library(pheatmap)
library(edgeR)
library(ggupset)
library(edgeR)
library(ggrepel)

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

###volcano plots
##set up the data
features<-c('Mir670hg','Ppm1h',
  'Baz1a','Tll1','Inhba',
  'Kdm2b')

volcano_de<-de.dge[c(1,3,11,15)]
names(volcano_de)<-c('CA1','CA3.1','GC','PS.1')
volcano_plots<-list()

##set up the plots
for(i in 1:length(volcano_de)){
data<-as.data.frame(volcano_de[[i]])
data <- data %>% 
  mutate(
    Expression = case_when(logFC >= 0 & FDR <= 0.05 ~ "Upregulated",
                           logFC <= 0 & FDR <= 0.05 ~ "Downregulated",
                           TRUE ~ "Unchanged")
  )
genes<-data[rownames(data) %in% features,]

volcano_plots[[i]]<-ggplot(data, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = Expression), size = 1/2) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_color_manual(values = c("indianred", "gray50", "cornflowerblue")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(colour='black'),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 13),
        legend.position='none') +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  ggtitle(names(volcano_de)[[i]])+
  geom_text_repel(data = genes,
                   mapping = aes(logFC, -log(FDR,10),
                                 label=rownames(genes),
                                 color=Expression),hjust=25) + 
  theme(legend.position='none') 
}
pdf('ca1_volcano.pdf',h=2.1,w=3.25)
volcano_plots[[1]]
dev.off()
pdf('ca3_volcano.pdf',h=2.1,w=3.25)
volcano_plots[[2]]
dev.off()
pdf('gc_volcano.pdf',h=2.1,w=3.25)
volcano_plots[[3]]
dev.off()
pdf('ps_volcano.pdf',h=2.1,w=3.25)
volcano_plots[[4]]
dev.off()


##make violin plots
##make new SCE object for CA3.D,CA1.V,and DG for the violin plots
sce.fig3<-sce.subset[,sce.subset$cellType_DE %in% c('CA1','CA3.1','PS.1','GC')]

##set up labels for plots
sce.fig3$cellType_DE<-droplevels(sce.fig3$cellType_DE)
sce.fig3$cellType:condition<-as.character(sce.fig3$cellType_DE)
sce.fig3$cellType:condition<-factor(ifelse(sce.fig3$condition=='Sham',
                                   paste0(sce.fig3$cellType:condition,'.Sham'),
                                   paste0(sce.fig3$cellType:condition,'.ECS')))
##Violin plots
features=c('Mir670hg','Ppm1h','Baz1a','Tll1','Inhba','Kdm2b')
pdf('fig3_violins.pdf',w=4,h=7)
plotExpression(sce.fig3,features=features,
               x="cellType:condition", colour_by="cellType_DE", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3)+
  theme(axis.text.x = element_text(angle = 45),
        text=element_text(size = 13))+
  labs(x='Cell Type:Condition',y='log2 normalized counts')
dev.off()



save(de.dge,file='processed_data/specific_de_results.rda')
