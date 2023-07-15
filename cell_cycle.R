setwd("./rds")
library(Seurat)
library(reshape2)
library(ggplot2)
library(ggsci)
library(RCircos)
library(Hmisc)
library(ggpubr)
library(plyr)
library(scales)
show_col(pal_simpsons(alpha = 0.8)(9))
cols <- c('G1' = '#709AE1CC', 'S' = '#FD7446CC', 'G2/M' = '#D2AF81CC')

ova <- readRDS('./OVA.rds')
s.genes <- capitalize(tolower(cc.genes$s.genes))
g2m.genes <- capitalize(tolower(cc.genes$g2m.genes))
ova <- CellCycleScoring(ova, s.features = s.genes, g2m.features = g2m.genes)
ova@meta.data$celltype <- ova@active.ident[rownames(ova@meta.data)]

cellcycle_bar <- function(ct){ #pass cell type name to this function
  sub_df <- subset(ova@meta.data, celltype == ct)
  sub_df_freq <- as.data.frame(table(sub_df$group, sub_df$Phase))
  names(sub_df_freq) <- c('group', 'cellcycle', 'number')
  sub_df_freq2 <- ddply(sub_df_freq, .(group), .fun = function(x)x[,'number']/sum(x[,'number']))
  names(sub_df_freq2) <- c('group', levels(sub_df_freq$cellcycle))
  sub_df_freq2 <- melt(sub_df_freq2)
  
  sub_df_freq2$group <- factor(sub_df_freq2$group, levels = c('OVA5D', 'OVA5W', 'OVA3M', 'OVA6M', 'OVA9M', 'OVA12M', 'OVA15M'))
  sub_df_freq2$variable <- as.character(sub_df_freq2$variable)
  sub_df_freq2$variable[sub_df_freq2$variable == 'G2M'] <- "G2/M"
  sub_df_freq2$variable <- factor(sub_df_freq2$variable, levels = c('G1', 'S', 'G2/M'))
  (p <- ggbarplot(sub_df_freq2, x = 'group', y = 'value', fill = 'variable', size = 0.2) +scale_fill_manual(values = cols)
    + theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12), 
            axis.text.y = element_text(color = 'black'),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            legend.text = element_text(size = 12))
    + labs(fill = 'Cell cycle', y = 'Cell proportion'))
  return(p)
}

