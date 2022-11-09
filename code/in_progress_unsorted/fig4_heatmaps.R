library(pheatmap)

load('de_dge.rda')

##pull out the FDRs for heatmaps
heat_all<-data.frame(
  'GC'<-de.dge[[6]]$logFC,
  'CA3'<-de.dge[[3]]$logFC,
  'PS'<-de.dge[[12]]$logFC,
  'CA1'<-de.dge[[1]]$logFC,
  row.names=rownames(de.dge[[1]])
)
colnames(heat_all)<-c('GC','CA3','PS','CA1')


##gene lists for heatmaps
synapse<-c('Grin2a','Sorcs3','Dlg4','Ntng2','Baiap2',
           'Grasp','Phactr1','Svop','Sv2b','Slc8a1')


signal<-c('Ntrk2','Bdnf','Brinp1','Ptgs2','Tgfb2','Ednrb','Gpr176',
          'Sh3gl2','Gpr68','Pde4b')


transfer<-c('Plk2','Akap13','Pim1','Cop1','Nrp1','Mapk4',
            'Nedd4l','Itpk1')

##make heatmaps
heat_syn<-heat_all[rownames(heat_all) %in% synapse,]
heat_syn<-heat_syn[match(synapse,rownames(heat_syn)),]

heat_sig<-heat_all[rownames(heat_all) %in% signal,]
heat_sig<-heat_sig[match(signal,rownames(heat_sig)),]

heat_trans<-heat_all[rownames(heat_all) %in% transfer,]
heat_trans<-heat_trans[match(transfer,rownames(heat_trans)),]


pdf('synapse.pdf',h=3,w=3)
pheatmap(heat_syn,cluster_cols=F,cluster_rows=F,
         scale='row',angle_col='45',fontsize=11)
dev.off()

pdf('signal.pdf',h=3,w=3)
pheatmap(heat_sig,cluster_cols=F,cluster_rows=F,
         scale='row',angle_col='45',fontsize=11)
dev.off()

pdf('transfer.pdf',h=3,w=3)
pheatmap(heat_trans,cluster_cols=F,cluster_rows=F,
         scale='row',angle_col='45',fontsize=11)
dev.off()




