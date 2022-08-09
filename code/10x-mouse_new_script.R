print(
  plotExpression(sce.total,features=c("Slc17a6","S100b","Dlk1","Tpbg","Fn1","Gpc3","Cbln4","Col5a2","Ly6g6e"),
                 x="label2", colour_by="level_1", point_alpha=0.5, point_size=.7,add_legend=F)+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                 width = 0.3)) +
    ggtitle(label=paste0(names(GABA_list)[i], " markers"))
)

mossy_list<-list(
  "Mossy1" = c("Ajap1", "Ap2s1", "Calb2", "Cda", "Chga", "Cntn6"),
  "Mossy2" = c("Col12a1", "Col19a1", "Crispld1", "Csf2rb2", "Cspg5", "Cthrc1"),
  "Mossy3" = c("Dhrs3", "Dkkl1", "Drd2", "Epb41l2", "Ephb6", "Exph5"),
  "Mossy4" = c("Fgf1", "Gal", "Glp1r", "Gnb4", "Grm8", "Hpcal1"),
  "Mossy5" = c("Nmb", "Shh", "Prmt2", "Rbpms2", "Rgs12", "Rspo1"),
  "Mossy6" = c("Serpinf1", "Slc41a3", "Sod3", "Thbs2", "Tm4sf1", "Ttn")
)

pdf("n7_mouse10x_mossy_markers_residuals.pdf", height=6,width=12)
for(i in 1:length(mossy_list)){
  print(
    plotExpression(sce.total,features=c(mossy_list[[i]]),
                   x="label2", colour_by="level_1", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(mossy_list)[i], " markers"))
  )}
dev.off()

DG_list<-list(
  "DG1" = c("Ablim3",   "Akap7",    "Arhgap20", "Btg2",     "C1ql2",    "Cald1"),
  "DG2" = c("Cygb",     "Dgkh",     "Dock10",   "Fam163b",  "Glis3",    "Igfbp5"),
  "DG3" = c("Il1rap",   "Jun",      "Carmil1",  "Lrrtm4",   "Maml2",    "Marcks"),
  "DG4" = c("Mef2c",    "Npnt",     "Pcdh8",    "Pde7b",    "Pitpnc1",  "Pitpnm2"),
  "DG5" = c("Plekha2",  "Prdm5",    "Prox1",    "Rfx3",     "Scn3a",    "Shisa9"),
  "DG6" = c("Sipa1l2",  "Slc29a4",  "Slc4a4", "Stxbp6",   "Sv2c"),
  "DG7" = c("Tiam1",    "Trpc6",    "Zfp536")
)

pdf("n7_mouse10x_DG_markers_residuals.pdf", height=6,width=12)
for(i in 1:length(mossy_list)){
  print(
    plotExpression(sce.total,features=c(DG_list[[i]]),
                   x="label2", colour_by="level_1", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(DG_list)[i], " markers"))
  )}
dev.off()

CA1_list<-list(
  "CA1_hipposeq1" = c("BC030500","Cds1","Fibcd1","Adgrl2","Man1a","Mpped1"),
  "CA1_hipposeq2"= c("Nov","Pex5l","Satb2","Zfp386"),
  "CA1_divseq1" = c("Doc2b","Gabrd","Grm2","C1ql3","St3gal1","Trpc6"),
  "CA1_divseq2" = c("Fxyd6","Marcksl1","Atp2b4","Ldb2","Calb1","Cnih3")
)

pdf("n7_mouse10x_CA1_markers_residuals.pdf", height=6,width=12)
for(i in 1:length(CA1_list)){
  print(
    plotExpression(sce.total,features=c(CA1_list[[i]]),
                   x="label2", colour_by="level_1", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(CA1_list)[i], " markers"))
  )}
dev.off()

CA2_list<-list(
  "CA2_1" = c("Adcy5", "Arg2", "Cacna2d3", "Cacng5", "Camk4", "Ccdc3"),         
  "CA2_2" = c("Clmn", "Ctsc", "Dusp5", "Eepd1", "Strip2", "Fam46a"),
  "CA2_3" = c("Fgf2", "Fgf5", "Tgfb1i1", "Gpr12", "Gsto1", "Lman2"),
  "CA2_4" = c("Lrrc8c", "Mmp14", "Ndrg3", "Necab2", "Ntsr2", "Plcxd3"),
  "CA2_5"= c("Prss23", "Pygo1", "Rgs14", "Rgs5", "S100b","Slc4a8"),
  "CA2_6" = c("Sostdc1","Srgap2","St8sia6","Stard5","Stxbp5l","Syce2"),
  "CA2_7" = c("Vcan","Vit","Zfp804a")
)

pdf("n7_mouse10x_CA2_markers_residuals.pdf", height=6,width=12)
for(i in 1:length(CA2_list)){
  print(
    plotExpression(sce.total,features=c(CA2_list[[i]]),
                   x="label2", colour_by="level_1", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(CA2_list)[i], " markers"))
  )}
dev.off()

CA2_list<-list(
  "CA2_1" = c("Adcy5", "Arg2", "Cacna2d3", "Cacng5", "Camk4", "Ccdc3"),         
  "CA2_2" = c("Clmn", "Ctsc", "Dusp5", "Eepd1", "Strip2", "Fam46a"),
  "CA2_3" = c("Fgf2", "Fgf5", "Tgfb1i1", "Gpr12", "Gsto1", "Lman2"),
  "CA2_4" = c("Lrrc8c", "Mmp14", "Ndrg3", "Necab2", "Ntsr2", "Plcxd3"),
  "CA2_5"= c("Prss23", "Pygo1", "Rgs14", "Rgs5", "S100b","Slc4a8"),
  "CA2_6" = c("Sostdc1","Srgap2","St8sia6","Stard5","Stxbp5l","Syce2"),
  "CA2_7" = c("Vcan","Vit","Zfp804a")
)

pdf("n7_mouse10x_CA2_markers_residuals.pdf", height=6,width=12)
for(i in 1:length(CA2_list)){
  print(
    plotExpression(sce.total,features=c(CA2_list[[i]]),
                   x="label2", colour_by="level_1", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(CA2_list)[i], " markers"))
  )}
dev.off()

CA3_list<-list(
  "CA3_1" = c("Mndal", "Adamts9", "Clmp", "Smoc2", "Slco2a1", "Gm20751"),
  "CA3_2" = c("Slc22a4", "P2rx7", "Crlf1", "Cdh24", "Evc", "Bok"),
  "CA3_3" = c("Dpf3", "Galnt3", "Ptgs2", "Synpr", "Rhebl1", "Rerg"),
  "CA3_4" = c("Nkd2", "Car4", "Adamts14", "Lnx1", "Ifi203", "Tbata"),
  "CA3_5" = c("Homer3", "Tspan17", "Zfpm2")
)

pdf("n7_mouse10x_CA3_markers_residuals.pdf", height=6,width=12)
for(i in 1:length(CA3_list)){
  print(
    plotExpression(sce.total,features=c(CA3_list[[i]]),
                   x="label2", colour_by="level_1", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(CA3_list)[i], " markers"))
  )}
dev.off()

GABA_list<-list(
  "GABA1" = c("Slc32a1", "Dlx1", "Dlx6os1", "Kcnmb2", "Arx", "Ptprm"),
  "GABA3" = c("Btbd11", "Maf", "Nxph1", "Erbb4", "Pnoc", "Hapln1"),
  "GABA4" = c("Gm13629", "Lhx6", "Alk", "Dlx1as", "Grik1", "Grip2"),
  "GABA5" = c("Rbms3", "Kcnip1", "Gad2", "Coro6", "Reln", "Cacna2d2"),
  "GABA5" = c("Gad1", "Arl4c", "Osbpl3", "Ank1", "Rims3", "Sema4g")
)

pdf("n7_mouse10x_GABA_markers_residuals.pdf", height=6, width=12)
for(i in 1:length(GABA_list)){
  print(
    plotExpression(sce.total,features=c(GABA_list[[i]]),
                   x="label2", colour_by="level_1", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3) +
      ggtitle(label=paste0(names(GABA_list)[i], " markers"))
  )}
dev.off()

##Just some fun markers...
pdf("n7_mouse10x_exploratory_markers_residuals.pdf", height=6, width=12)
  print(
    plotExpression(sce.total,features=c('Nrgn','Slc17a6','Slc17a7','Slc17a8','Nts','Ntng2','Fn1','Sntg2'),
                   x="label2", colour_by="level_1", point_alpha=0.5, point_size=.7,add_legend=F)+
      stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                   width = 0.3)
  )
dev.off()

colData(sce.total)$level_1<-
  factor(
    ifelse(colData(sce.total)$label2 %in% c(2,3,4,5,9,11,12,13,
                                                   14,15,17,19,20,21,
                                                   22,23,24,25,26,27,28,29,
                                                   31),
           'CA_cells',
           ifelse(colData(sce.total)$label2 %in% c(6,8),
                  'DG_cells',
                         'GABA_cells')))
##From plots it looks like 20=CA2, 5,19=CA3, 26=mossy
colData(sce.total)$level_2<-as.character(colData(sce.total)$level_1)
colData(sce.total)$level_2 <- factor(
    ifelse(colData(sce.total)$label %in% c(5,19), 'CA3_cells',
           ifelse(colData(sce.total)$label==20 , 'CA2_cells',
                  ifelse(colData(sce.total)$label==26, 'mossy_cells',
                         colData(sce.total)$level_2))))

##Plot expression of subiculum projection neuron markers
subiculum_list<-list(
  "Subiculum" = c("Slc17a7", "Nts", "Fn1", "Slc17a6", "Slc17a6", "Cntn6")
)
pdf("n7_mouse10x_subiculum_markers2_residuals.pdf", height=6, width=12)
print(
  plotExpression(sce.total,features=c("Slc17a6", "Nts", "Fn1", "S100b","Tpbg", "Dlk1","Gpc3","Col5a2","Ly6g6e","Cbln4"),
                 x="label", colour_by="level_2", point_alpha=0.5, point_size=.7,add_legend=F)+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                 width = 0.3) +
    ggtitle(label="subiculum markers"))

dev.off()

colData(sce.total)$level_3<-as.character(colData(sce.total)$level_2)
colData(sce.total)$level_3 <- factor(
         ifelse(colData(sce.total)$label==11 , 'subiculum_cells',
                ifelse(colData(sce.total)$label==2, 'presubiculum_cells',
                       colData(sce.total)$level_3)))

### Doublet detection===========================================================
library(BiocSingular)
set.seed(100)

# Setting up the parameters for consistency with denoisePCA();
# this can be changed depending on your feature selection scheme.
dbl.dens <- computeDoubletDensity(sce.total, subset.row=hdg, 
                                  d=50)
summary(dbl.dens)
sce.total$doubletScore <- dbl.dens
plotUMAP(sce.total, colour_by="doubletScore", text_by="label")
##seems like low doublets/they're dispersed across clusters

g25 <- buildSNNGraph(sce.subset, k=25, use.dimred = 'PCA')
clust25 <- igraph::cluster_walktrap(g25)$membership
colData(sce.subset)$k_25_label <- factor(clust25)
table(colData(sce.subset)$k_25_label)

sce.subset$cell.type<-factor(
  ifelse(sce.subset$label==6,'CA1',
      ifelse(sce.subset$label==13,'CA2',
          ifelse(sce.subset$label==7,'CA3.1',
              ifelse(sce.subset$label==12,'CA3.2',
                  ifelse(sce.subset$label==16,'Mossy',
                      ifelse(sce.subset$label==2,'GABA.1',
                             ifelse(sce.subset$label==9,'GABA.2',
                                    ifelse(sce.subset$label==10,'GABA.3',
                                           ifelse(sce.subset$label==5,'DG.1',
                                                  ifelse(sce.subset$label==8,'DG.2',
                                                         ifelse(sce.subset$label==11,'Subiculum',
                                                                ifelse(sce.subset$label==4,'L5/6/NP.RHP',
                                                                       ifelse(sce.subset$label==3,'L5/6/CT.RHP',
                                                                              ifelse(sce.subset$label==1,'PPP.RHP',
                                                                                     ifelse(sce.subset$label==14,'L3.RHP',
                                                                                            ifelse(sce.subset$label==15,'L2/3.RHP.1',
                                                                                                   ifelse(sce.subset$label==17,'L2/3.RHP.2',
                                                                                                          'L2/3.RHP.3'
  )))))))))))))))))
)

neurogenesis.features<-c('Id4', 'Thrsp', 'Apoe', 'Rest', 'Hopx', 'Sox2', 'Gli1', 'Neurog2', 'Ascl1', 'Nes', 'Mki67', 'Elavl2', 'Cdk4', 'Mcm2', 'Pcna', 'Eomes', 'Sox11', 'Neurod1', 'Prox1', 'Dpysl3', 'Dcx', 'Ntrk2', 'Slc12a2', 'Calb2', 'Ncam1', 'Egr1', 'Calb1', 'Rbfox3', 'Creb1', 'Slc12a1')
plotHeatmap(sce.subset,features=neurogenesis.features,order_columns_by="cell.type",
            colour_columns_by="cell.type",main="Adult neurogenesis-associated gene expression")