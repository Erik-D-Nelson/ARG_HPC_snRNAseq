#cd ~/arg-hpc/
####Code for figure 5: scCoGAPS
library(CoGAPS)
library(scater)
library(scran)
library(SingleCellExperiment)
library(schex)
library(here)
library(ggpubr)
library(reshape2)

##scCoGAPS has already finished running at this point. Let's load the results
load("2022_scCoGAPS.rda")
load("onehot_cellType.rda")
load("onehot_condition.rda")
ls()
#[1] "params"  "results"
results
#class: CogapsResult 
#dim: 15990 57 
#metadata(8): meanChiSq firstPass ... params version
#rownames: NULL
#colnames: NULL
#factorData names(0):

##We need to make 3, maybe four plots. The first is just a dotplot where 
##size=cell type correlation and color=ECS correlation
##One hot encode annotation and condition
#First, load sce subset
load('newSceSubset_geneFiltered_under20.rda')
sce.subset$cellType<-as.character(sce.subset$cellType)
sce.subset$cellType<-recode(sce.subset$cellType, 'ProS.1'='PS.1','ProS.2'='PS.2',
                            'RHP.L2/3'='L2/3','RHP.L5/Po'='L5/Po',
                            'RHP.L6/6b'='L6/6b','GABA.CGE.1'='GABA.1',
                            'GABA.CGE.2'='GABA.2','GABA.MGE.1'='GABA.3',
                            'GABA.MGE.2'='GABA.4','GABA.MGE.3'='GABA.5')

sce.subset$cellType<-factor(sce.subset$cellType)#, levels=c("DG.1","DG.2","CA4",
                                                  #        "CA1",'PS.1','PS.2',
                                                   #       'L2/3','L5/Po','L6/6b',
                                                  #        'GABA.1','GABA.2','GABA.3',
                                                  #        'GABA.4','GABA.5'))

cellType<-data.frame('ID'=colnames(sce.subset),'cellType'=sce.subset$cellType)


colData(sce.subset)<-cbind(colData(sce.subset),results@sampleFactors)

breaks<-seq(-0.1,1,by=.1)
yl<-as.vector(y)
pat_order<-c(26,46,55,23,15,19,20,51,32,56,47,
             27,10,11,57,52,48,2,33,28,31,29,
             40,49,35,30,50,12,42,53,36,43,37,
             34,41,38,39,9,13,44,54,45,1,16,17,
             6,18,14,3,8,22,21,5,24,7,25,4)
type_order<-c(6,7,1,2,3,4,5,16,17,18,13,14,15,8,9,10,11,12)
y<-cor(x,results@sampleFactors)
z<-cor(condition,results@sampleFactors)
y<-y[rev(type_order),pat_order]
z<-z[,pat_order]
melted<-melt(y)
ECS<-z[2,]
guide<-as.character(unique(melted$Var2))
ECS<-ECS[names(ECS) %in% guide]
ECS<-rep(ECS,each=18)
melted$ECS<-ECS

pdf('test_balloon.pdf',height=4,width=9)
ggplot(melted, aes(x =Pattern, y = Annotation)) + 
  geom_point(aes(size=cor_annotation,color=cor_ECS)) +
  theme(panel.background=element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA),
        axis.text.x = element_text(angle = 90),
        text = element_text(size = 10)) +
  scale_radius(breaks=c(0,.25,.5,.75,1),limits=c(-0.25,1),range=c(0,5)) +
  viridis::scale_color_viridis(option='plasma',breaks=c(-0.4,0,0.4), limits=c(-0.41,0.41)) +
  labs(size=expression(paste("Pearson's ", 
                             italic("r"), 
                             ", pattern and annotation")),
       color=expression(paste("Pearson's ", 
                              italic("r"), 
                              ", pattern and ECS")))
dev.off()

names(melted)<-c("Annotation","Pattern","cor_annotation","cor_ECS")

hex <- make_hexbin(sce.subset, nbins = 100, 
                           dimension_reduction = "UMAP", use_dims=c(1,2))

label_df <- make_hexbin_label(hex, col="cellType")


pdf('hex_23.pdf',h=4,w=4)
plot_hexbin_meta(hex,col='Pattern_23',action='median') + 
  labs(title='Pattern 23',
       x = 'UMAP1', y='UMAP2',
       color='Pattern weights') + 
  theme(text = element_text(size = 9))+ 
  theme(text = element_text(size = 9)) + ggrepel::geom_label_repel(data = label_df, aes(x=x, y=y, label = label), 
                                                                   colour="black",  label.size = NA, fill = NA)+
  scale_fill_viridis(option='turbo')
dev.off()

pdf('hex_26.pdf',h=4,w=4)
plot_hexbin_meta(hex,col='Pattern_26',action='median') + 
  labs(title='Pattern 26',
       x = 'UMAP1', y='UMAP2',
       color='Pattern weights') + 
  theme(text = element_text(size = 9))+ 
  theme(text = element_text(size = 9)) + ggrepel::geom_label_repel(data = label_df, aes(x=x, y=y, label = label), 
                                                                   colour="black",  label.size = NA, fill = NA)+
  scale_fill_viridis(option='turbo')
dev.off()

pdf('hex_46.pdf',h=4,w=4)
plot_hexbin_meta(hex,col='Pattern_46',action='median') + 
  labs(title='Pattern 46',
       x = 'UMAP1', y='UMAP2',
       color='Pattern weights') + 
  theme(text = element_text(size = 9))+ 
  theme(text = element_text(size = 9)) + ggrepel::geom_label_repel(data = label_df, aes(x=x, y=y, label = label), 
                                                                   colour="black",  label.size = NA, fill = NA)+
  scale_fill_viridis(option='turbo')
dev.off()

pdf('hex_55.pdf',h=4,w=4)
plot_hexbin_meta(hex,col='Pattern_55',action='median') + 
  labs(title='Pattern 55',
       x = 'UMAP1', y='UMAP2',
       color='Pattern weights') + 
  theme(text = element_text(size = 9))+ 
  theme(text = element_text(size = 9)) + ggrepel::geom_label_repel(data = label_df, aes(x=x, y=y, label = label), 
                                                                   colour="black",  label.size = NA, fill = NA)+
  scale_fill_viridis(option='turbo')
dev.off()





