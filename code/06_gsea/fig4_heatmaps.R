library(pheatmap)
library(org.Mm.eg.db)

load('de_dge.rda')

##pull out the FDRs for heatmaps
heat_all<-data.frame(
  'GC'<-de.dge[[11]]$logFC,
  'CA3'<-de.dge[[3]]$logFC,
  'PS'<-de.dge[[15]]$logFC,
  'CA1'<-de.dge[[1]]$logFC,
  row.names=rownames(de.dge[[1]])
)
colnames(heat_all)<-c('GC','CA3','PS','CA1')
#heat_all[is.na(heat_all)]<-0

is.de <- decideTestsPerLabel(de.dge, threshold=0.05)
is.de<-as.data.frame(is.de)
DEG.GC<-rownames(is.de)[is.de$DG %in% c(1)]
DEG.CA1<-rownames(is.de)[is.de$CA1 %in% c(1)]
DEG.CA3.1<-rownames(is.de)[is.de$CA3.1 %in% c(1)]
DEG.PS.1<-rownames(is.de)[is.de$ProS.1 %in% c(1)]


syn<-rownames(cc_parent_df[[1]])
growth<-rownames(bp_parent_df[[11]])
ubiq<-rownames(mf_parent_df[[8]])
kinases<-rownames(mf_parent_df[[2]])

GC<-de.dge[[11]][rownames(de.dge[[11]]) %in% DEG.GC,]
CA1<-de.dge[[1]][rownames(de.dge[[1]]) %in% DEG.CA1,]
CA3<-de.dge[[3]][rownames(de.dge[[3]]) %in% DEG.CA3,]
PS<-de.dge[[15]][rownames(de.dge[[15]]) %in% DEG.PS,]

##function to find genes with heterogeneous DE patterns across cell types
find_genes<-function(x,y,z){
  table<-relocate(heat_all,x)
  table<-table[rownames(table) %in% y,]
  table$meanRatio<-table[,1]-pmax(table[,2],table[,3],table[,4])
  table<-table[order(table$meanRatio,decreasing=T),]
  retrieved <- AnnotationDbi::select(org.Mm.eg.db, keytype="GO", 
                                     keys=z, columns="SYMBOL")$SYMBOL
  table<-table[rownames(table) %in% retrieved,]
  return(table)
}


##gene lists for heatmaps
synapse<-c('Arc','Nptx2','Grin2a','Sorcs3','Dlg4','Baiap2','Iqsec3',
           'Grasp','Ppm1h','Phactr1','Svop','Sv2b','Sh3gl2','Slc8a1')


signal<-c('Inhba','Ptgs2','Bdnf','Ntrk2','Tgfb2','Rhoq',
          'Smad3','Gclc','Gcnt2','Twf2','Pde4b','Nos1')


transfer<-c('Cop1','Fbxw7','Ube2g1','Cul3','Ube2q2','Mib1',
            'Znrf3','Rnf38','Bcor','Rmnd5b','Rnf150','Nedd4l')

kinase<-c('Pim1','Plk2','Akap13','Mapk4','Stk40','Rps6ka3',
          'Pi4ka','Tlk1','Itpka','Igf1r','Mapk11','Itpk1')
##make heatmaps
heat_syn<-heat_all[rownames(heat_all) %in% synapse,]
heat_syn<-heat_syn[match(synapse,rownames(heat_syn)),]

heat_sig<-heat_all[rownames(heat_all) %in% signal,]
heat_sig<-heat_sig[match(signal,rownames(heat_sig)),]

heat_trans<-heat_all[rownames(heat_all) %in% transfer,]
heat_trans<-heat_trans[match(transfer,rownames(heat_trans)),]

heat_kinase<-heat_all[rownames(heat_all) %in% kinase,]
heat_kinase<-heat_kinase[match(kinase,rownames(heat_kinase)),]

heat<-rbind(heat_syn,heat_sig,heat_trans,heat_kinase)

pdf('heat_all.pdf',h=8.5,w=3)
pheatmap(heat,cluster_cols=F,cluster_rows=F,
         scale='row',angle_col='45',fontsize=11,na_col='darkblue')
dev.off()
