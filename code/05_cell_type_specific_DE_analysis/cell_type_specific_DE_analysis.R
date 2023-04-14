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
library(viridis)

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

###store sum degs for later
sum_deg<-x$sum
names(sum_deg)<-rownames(x)

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

##pvalue histograms
dist<-list()
for(i in 1:length(de.dge)){  
  dist[[i]]<-ggplot(as.data.frame(de.dge[[i]]),aes(x=PValue))+
    geom_histogram(bins=50)+
    ggtitle(paste0('Gene PValue Histogram: ',names(de.dge)[i]))}

pdf('plots/figS4/figS4.pdf',h=15,w=15)
grid.arrange(grobs=dist)
dev.off()

####figure s5 bootstrapping and plots
set.seed(19210)
# Function to bootstrap a specific cellType_DE to a desired number
bootstrap_cellType_DE_with_sample_proportions <- function(sce, cellType_DE, target_count, n_iterations) {
  cells_of_type <- sce[, sce$cellType_DE == cellType_DE]
  sample_counts <- table(cells_of_type$Sample)
  
  # Calculate the proportion of cells in each sample for the given cellType_DE
  sample_proportions <- sample_counts / sum(sample_counts)
  
  # Calculate the target number of cells for each sample level
  target_sample_counts <- round(sample_proportions * target_count)
  
  # Bootstrap cells for each sample level for n_iterations
  bootstrapped_sce_list <- lapply(1:n_iterations, function(iter) {
    message(paste("Sampling iteration", iter))
    bootstrapped_sce_iteration <- lapply(names(target_sample_counts), function(sample_name) {
      cells_to_sample <- which(cells_of_type$Sample == sample_name)
      target_sample_count <- target_sample_counts[[sample_name]]
      
      # Sample with replacement
      cells_to_keep_after_bootstrapping <- sample(cells_to_sample, target_sample_count, replace = TRUE)
      
      sce_subset <- cells_of_type[, cells_to_keep_after_bootstrapping]
      
      return(sce_subset)
    })
    
    # Combine the bootstrapped SingleCellExperiment objects for the given cellType_DE in the current iteration
    bootstrapped_sce_combined <- do.call(cbind, bootstrapped_sce_iteration)
    
    return(bootstrapped_sce_combined)
  })
  
  return(bootstrapped_sce_list)
}

# Find the lowest number of cells in any cellType_DE
min_cell_count <- min(table(sce.subset$cellType_DE))

# Get unique cell types_DE
unique_cellTypes_DE <- unique(sce.subset$cellType_DE)

# Number of bootstrap iterations
n_iterations <- 250

# Bootstrap each cellType_DE to the lowest number of cells, preserving sample proportions, for n_iterations
bootstrapped_sce_list <- lapply(unique_cellTypes_DE, function(ct) {
  bootstrap_cellType_DE_with_sample_proportions(sce.subset, ct, min_cell_count, n_iterations)
})

library(edgeR)

# Function to run DE analysis and return counts of upregulated and downregulated genes
run_de_analysis <- function(sce) {
  # Aggregate data by Sample and cellType_DE
  message("Aggregating data by Sample and cellType_DE...")
  message(Sys.time())
  summed <- aggregateAcrossCells(sce, ids=colData(sce)[c('Sample', 'cellType_DE')])
  
  # Perform DE analysis
  message("Performing DE analysis...")
  message(Sys.time())
  de.dge <- pseudoBulkDGE(summed,
                          label=summed$cellType_DE,
                          design=~condition,
                          coef=2,
                          condition=summed$condition)
  is.de <- decideTestsPerLabel(de.dge, threshold=0.05)
  summary_per_label <- summarizeTestsPerLabel(is.de)
  
  message("DE analysis completed.")
  message(Sys.time())
  return(summary_per_label)
}


# Run DE analysis for each bootstrap iteration and store the results in a list
de_results_list <- lapply(bootstrapped_sce_list, function(cellType_bootstraps) {
  lapply(cellType_bootstraps, run_de_analysis)
})
save(de_results_list,file='bootstrapping_shit2.rda')




de_list<-de_results_list
for (i in 1:length(de_list)) {
  for (j in 1:length(de_list[[i]])) {
    de_list[[i]][[j]] <- de_results_list[[i]][[j]][, match(colnames(de_results_list[[2]][[1]]), colnames(de_results_list[[i]][[j]]))]
  }
}

# Replace all NA values with 0 in de_list
de_list <- lapply(de_list, function(list_of_tables) {
  lapply(list_of_tables, function(table) {
    table[is.na(table)] <- 0
    table
  })
})

# Sum columns 1 and 3 of each table in de_list
summed_list <- lapply(de_list, function(list_of_tables) {
  lapply(list_of_tables, function(table) {
    sum(table[c(1, 3)], na.rm = TRUE)
  })
})

# Create a new empty list object
x <- list()

# Assign values to the elements of the list using a for loop
for (i in seq_along(summed_list)) {
  x[[i]] <- unlist(summed_list[[i]])
}

names(x)<-unique_cellTypes_DE
x<-x[c(17,4,11,14,15,5,8,7,3,12,6,10,16,1,9,13,2)]

data_frame <- bind_rows(lapply(seq_along(x), function(i) {
  data.frame(Group = factor(i), Value = unlist(x[[i]]))
}))
data_frame2<-data_frame
data_frame2$Value<-data_frame2$Value+0.5


library(ggplot2)


# Set the same limits for the x and y axes
#xlims <- log2(range(unlist(x)))
#ylims <- c(0, max(sapply(x, density)$y))

# Loop over the elements of x and make a density plot for each
#for (i in seq_along(x)) {
# Calculate the density estimate for the current element
dens <- density(x[[i]])
custom_colors<-viridis_pal(option='H')(17)
# Create a ggplot2 object for the density plot
boxplot <- ggplot(data_frame, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() +
  scale_y_continuous(trans = "pseudo_log",breaks=c(1,10,100,1000)) +
  scale_x_discrete(labels = names(x)) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Distribution of DE genes by cell type after downsampling (250 bootstrapped iterations)", x = "Annotated Cell Type", y = "DEGs") +
  guides(fill = guide_legend(override.aes = list(size=1.5))) +
  theme_minimal()+
  theme( 
    axis.line = element_line(colour = "black"),
    text=element_text(colour='black'),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    legend.position='none',
    axis.text.x = element_text(angle = 45, hjust = 1,colour='black'))
pdf('DEG_boxplot_downsampling.pdf',h=4,w=6)
boxplot
dev.off()


# Prepare the data
cell_type_counts <- data.frame(CellType = names(sum_deg), DE_Genes = sum_deg)
cell_type_sample_counts <- data.frame(CellType = names(table(sce.subset$cellType_DE)), Sample_Count = as.numeric(table(sce.subset$cellType_DE)))
df <- left_join(cell_type_sample_counts, cell_type_counts, by = "CellType")
df<-df[order(df$DE_Genes),]
names(custom_colors)<-df$CellType

# Fit a linear model to the log-transformed data
fit <- lm(DE_Genes ~ Sample_Count, data = df)

# Extract the coefficients and R2 value
intercept <- coef(fit)[1]
slope <- coef(fit)[2]
r2 <- summary(fit)$r.squared

# Create the scatter plot with a regression line and overall R2 and equation
scatter<-ggplot(df, aes(x = Sample_Count, y = DE_Genes, color = CellType)) +
  geom_point(aes(color = CellType), size = 3) +
  scale_color_manual(values = custom_colors[c(match(df$CellType, names(custom_colors)))]) +
  geom_text_repel(aes(label = CellType),max.overlaps = 15) +
  labs(title = "DE Genes per cell type vs nuclei per cell type",
       x = "Number of nuclei",
       y = "Number of DEGs") +
  theme_minimal() +
  scale_y_continuous(trans = "pseudo_log", breaks = c(1, 10, 100, 1000, 10000)) +
  scale_x_continuous(trans = "log10", breaks = c(100, 750, 5000)) +
  geom_smooth(method = "glm", se = F, color = "black",linetype='dashed',size=0.5) +
  theme(
    axis.line = element_line(colour = "black"),
    text = element_text(colour = 'black'),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    legend.position = 'none') +
  annotate("text", x = 150, y = 1000, 
           label = sprintf("Equation: y = %.2f + %.2fx)", intercept, slope),
           size = 4) +
  annotate("text", x = 150, y = 500, 
           label = sprintf("R2 = %.2f", r2), 
           size = 4)
pdf('DEG_scatterplot.pdf',h=4,w=6)
scatter
dev.off()
