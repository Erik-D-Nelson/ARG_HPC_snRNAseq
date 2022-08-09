pdf("plotExpression_bdnf_shamonly.pdf", height=6,width=12)
print(
  plotExpression(x,features='Bdnf',
                 x="level_6", colour_by="level_6", point_alpha=0.5, point_size=.7,add_legend=F)+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
                 width = 0.3))
dev.off()