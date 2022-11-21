
####Code for figure 5: scCoGAPS
library(CoGAPS)
library(scater)
library(scran)
library(SingleCellExperiment)
library(schex)
library(here)
library(ggpubr)
library(reshape2)
library(pheatmap)

##scCoGAPS has already finished running at this point. Let's load the results
load("processed_data/2022_scCoGAPS.rda")
load("processed_data/sce_subset.rda")

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


##make one hot encoded dataframes for celltype and condition
#celltype
data<-as.data.frame(sce.subset$cellType)
colnames(data)<-'cellType'
onehot_celltype <-  dcast(data = data, rownames(data) ~ cellType, length)
rownames(onehot_celltype)<-onehot_celltype[,1]
onehot_celltype[,1]<-as.numeric(onehot_celltype[,1])
onehot_celltype<-onehot_celltype[order(onehot_celltype[,1],decreasing=F),]
onehot_celltype[,1]<-NULL

#condition
data<-as.data.frame(sce.subset$condition)
colnames(data)<-'condition'
onehot_condition <- dcast(data = data, rownames(data) ~ condition, length)
rownames(onehot_condition)<-onehot_condition[,1]
onehot_condition[,1]<-as.numeric(onehot_condition[,1])
onehot_condition<-onehot_condition[order(onehot_condition[,1],decreasing=F),]
onehot_condition[,1]<-NULL

##set up balloon plot

pat_order<-c(26,46,55,23,15,19,20,51,32,56,47,
             27,10,11,57,52,48,2,33,28,31,29,
             40,49,35,30,50,12,42,53,36,43,37,
             34,41,38,39,9,13,44,54,45,1,16,17,
             6,18,14,3,8,22,21,5,24,7,25,4)
type_order<-c(1,2,7,6,5,4,3,8,9,10,11,12,13,14,15,16,17,18)
y<-cor(onehot_celltype,results@sampleFactors)
z<-cor(onehot_condition,results@sampleFactors)
y<-y[rev(type_order),pat_order]
z<-z[,pat_order]
melted<-melt(y)
ECS<-z[2,]
ECS<-rep(ECS,each=18)
melted$ECS<-ECS
names(melted)<-c("Annotation","Pattern","cor_annotation","cor_ECS")


pdf('plots/fig5/nmf_balloon_plot.pdf',height=4,width=9)
ggplot(melted, aes(x =Pattern, y = Annotation)) + 
  geom_point(aes(size=cor_annotation,color=cor_ECS)) +
  theme(panel.background=element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA),
        axis.text.x = element_text(angle = 90),
        text = element_text(size = 10)) +
  scale_radius(breaks=c(0,.25,.5,.75,1),limits=c(-0.25,1),range=c(0,5)) +
  scale_x_discrete(labels=pat_order)+
  viridis::scale_color_viridis(option='plasma',breaks=c(-0.4,0,0.4), limits=c(-0.41,0.41)) +
  labs(size=expression(paste("Pearson's ", 
                             italic("r"), 
                             ", pattern and annotation")),
       color=expression(paste("Pearson's ", 
                              italic("r"), 
                              ", pattern and ECS")))
dev.off()


hex <- make_hexbin(sce.subset, nbins = 100, 
                           dimension_reduction = "UMAP", use_dims=c(1,2))

label_df <- make_hexbin_label(hex, col="cellType")


pdf('hex_23.pdf',h=3,w=4)
plot_hexbin_meta(hex,col='Pattern_23',action='median') + 
  labs(title='Pattern 23',
       x = 'UMAP1', y='UMAP2',
       color='Pattern weights') + 
  theme(text = element_text(size = 9))+ 
  theme(text = element_text(size = 9)) + ggrepel::geom_label_repel(data = label_df, aes(x=x, y=y, label = label), 
                                                                   colour="black",  label.size = NA, fill = NA)+
  scale_fill_viridis(option='turbo')
dev.off()

pdf('hex_26.pdf',h=3,w=4)
plot_hexbin_meta(hex,col='Pattern_26',action='median') + 
  labs(title='Pattern 26',
       x = 'UMAP1', y='UMAP2',
       color='Pattern weights') + 
  theme(text = element_text(size = 9))+ 
  theme(text = element_text(size = 9)) + ggrepel::geom_label_repel(data = label_df, aes(x=x, y=y, label = label), 
                                                                   colour="black",  label.size = NA, fill = NA)+
  scale_fill_viridis(option='turbo')
dev.off()

pdf('hex_46.pdf',h=3,w=4)
plot_hexbin_meta(hex,col='Pattern_46',action='median') + 
  labs(title='Pattern 46',
       x = 'UMAP1', y='UMAP2',
       color='Pattern weights') + 
  theme(text = element_text(size = 9))+ 
  theme(text = element_text(size = 9)) + ggrepel::geom_label_repel(data = label_df, aes(x=x, y=y, label = label), 
                                                                   colour="black",  label.size = NA, fill = NA)+
  scale_fill_viridis(option='turbo')
dev.off()

pdf('hex_55.pdf',h=3,w=4)
plot_hexbin_meta(hex,col='Pattern_55',action='median') + 
  labs(title='Pattern 55',
       x = 'UMAP1', y='UMAP2',
       color='Pattern weights') + 
  theme(text = element_text(size = 9))+ 
  theme(text = element_text(size = 9)) + ggrepel::geom_label_repel(data = label_df, aes(x=x, y=y, label = label), 
                                                                   colour="black",  label.size = NA, fill = NA)+
  scale_fill_viridis(option='turbo')
dev.off()


heat<-results@featureLoadings
heat<-heat[,c(26,23,46,55)]
prg<-c('Atf3','Ptgs2','Baz1a','Brinp1','Bdnf','Maml3','Rheb','Nrn1','Gadd45g',
            'Arc', 'Egr3', 'Grasp', 'Sv2c', 'Sult2b1','Sstr2','Mbnl2','Cystm1',
            'Sertad1','Slc2a3','Gpr22','Errfi1','Fndc9','Osgin2','Fam126b')

heat_prg<-heat[rownames(heat) %in% prg,]
heat_prg<-heat_prg[rev(match(prg,rownames(heat_prg))),]

pdf('nmf_heatmap_prgs.pdf',h=3,w=2)
pheatmap(heat_prg,cluster_cols=F,cluster_rows=F,scale='row')
dev.off()

srg<-c('Sema3e','Nptx2','Kcna4','Acan','Kcnj4','Gpr3','Hdac9','Kcna1','Gpr63',
       'Trip10','Stac','Slitrk4','Mfap3l','Bmp3','Cd109',
       'Fezf2','Palmd','Hunk','Stc2','Kcnh8')
heat_srg<-heat[rownames(heat) %in% srg,]
heat_srg<-heat_srg[rev(match(srg,rownames(heat_srg))),]

pdf('nmf_heatmap_srgs.pdf',h=3,w=2)
pheatmap(heat_srg,cluster_cols=F,cluster_rows=F,scale='row')
dev.off()

synapse<-c('Tanc2','Acvr1','Sntb2','Nf1','Sorcs3','Epha10','Susd4','Ston2','Lgi1',
           'Efhd2','Actn4','Nectin1','Nrxn2','Lrfn2','Syn2','Adgra2','Plec',
           'Svop','Amph')
heat_synapse<-heat[rownames(heat) %in% synapse,]
heat_synapse<-heat_synapse[rev(match(synapse,rownames(heat_synapse))),]

pdf('nmf_heatmap_synapse.pdf',h=3,w=2)
pheatmap(heat_synapse,cluster_cols=F,cluster_rows=F,scale='row')
dev.off()


