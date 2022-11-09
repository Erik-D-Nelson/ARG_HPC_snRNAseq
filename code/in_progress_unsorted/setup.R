
bp_total<-list()
for(i in 1:length(bp)){
  bp_total[[i]]<-rownames(bp[[i]])
}

bp_total<-unique(do.call(c, bp_total))


for(i in 1:length(bp)){
  bp[[i]]<-bp[[i]][match(bp_total,rownames(bp[[i]])),]
}

fdr<-list()
for(i in 1:length(bp)){
  fdr[[i]]<-bp[[i]]$p.adjust
}

fdr_df<-do.call(cbind.data.frame, fdr)
fdr_df[is.na(fdr_df)]<-1
fdr_df= -log10(fdr_df)
fdr_df$mean<-rowMeans(fdr_df)
bp_scores<-fdr_df$mean
names(bp_scores)<-bp_total





cc_total<-list()
for(i in 1:length(cc)){
  cc_total[[i]]<-rownames(cc[[i]])
}

cc_total<-unique(do.call(c, cc_total))


for(i in 1:length(cc)){
  cc[[i]]<-cc[[i]][match(cc_total,rownames(cc[[i]])),]
}

fdr<-list()
for(i in 1:length(cc)){
  fdr[[i]]<-cc[[i]]$p.adjust
}

fdr_df<-do.call(cbind.data.frame, fdr)
fdr_df[is.na(fdr_df)]<-1
fdr_df= -log10(fdr_df)
fdr_df$mean<-rowMeans(fdr_df)
cc_scores<-fdr_df$mean
names(cc_scores)<-cc_total


mf_total<-list()
for(i in 1:length(mf)){
  mf_total[[i]]<-rownames(mf[[i]])
}

mf_total<-unique(do.call(c, mf_total))


for(i in 1:length(mf)){
  mf[[i]]<-mf[[i]][match(mf_total,rownames(mf[[i]])),]
}

fdr<-list()
for(i in 1:length(mf)){
  fdr[[i]]<-mf[[i]]$p.adjust
}

fdr_df<-do.call(cbind.data.frame, fdr)
fdr_df[is.na(fdr_df)]<-1
fdr_df= -log10(fdr_df)
fdr_df$mean<-rowMeans(fdr_df)
mf_scores<-fdr_df$mean
names(mf_scores)<-mf_total


reduced_bp <- reduceSimMatrix(bp_simmatrix,
                                bp_scores,
                                threshold=0.9,
                                orgdb="org.Mm.eg.db")

reduced_cc <- reduceSimMatrix(simMatrix_CC,
                              cc_scores,
                              threshold=0.9,
                              orgdb="org.Mm.eg.db")

reduced_mf <- reduceSimMatrix(mf_simmatrix,
                              mf_scores,
                              threshold=0.9,
                              orgdb="org.Mm.eg.db")
