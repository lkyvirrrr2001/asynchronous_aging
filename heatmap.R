library(Seurat)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
cal_ratio2 <- function(obj,genes,g){
  ratios <- c()
  ova_tmp <- subset(obj, group == g)
  for (i in genes){
    genes_expr_cluster <- ova_tmp[i,]
    a <- as.matrix(genes_expr_cluster[["RNA"]]@data)
    gene1_exp <- as.character(a[1,])
    exp_num <- length(gene1_exp[gene1_exp!=0])
    ratio <- exp_num/length(gene1_exp)
    ratios <- c(ratios, ratio)
  }
  
  return(ratios)

}
mean_expression2 <- function(obj,g){
  genes_expr_5d <- subset(obj, group == g)
  a <- AverageExpression(genes_expr_5d, group.by = 'active.ident')#, layer = 'scale.data'
  mean_xp <- a$RNA
  return(mean_xp)
}

ova <- readRDS('../../../../ThesisofAging/公司回传的/三个大类rds/OVA.rds')
ova <- UpdateSeuratObject(object = ova)
ova$active.ident <- ova@active.ident
tfs <- read.csv('../tf_analysis2/ova_tf_freq.txt', sep = '\t')
ova <- ova[tfs$Var1,]

a <- DotPlot(ova,features = row.names(ova), assay = 'RNA', split.by = 'group', cols = brewer.pal(n = 7, 'YlGnBu'))
a_df <- a$data
row.names(a_df) <- paste(a_df$features.plot, a_df$id, sep = '_')
exp_df_hm <- data.frame()
for (g in c('OVA5D', 'OVA5W', 'OVA3M', 'OVA6M', 'OVA9M', 'OVA12M', 'OVA15M')){
  menxp <- mean_expression2(ova, g)
  menxp <- as.data.frame(menxp)
  menxp$gene <- row.names(menxp)
  tmp_df <- melt(menxp, value.name = 'gene')
  names(tmp_df)[2:3] <- c('celltype', 'mean_expression') 
  tmp_df$group <- g
  exp_df_hm <- rbind(exp_df_hm, tmp_df)
}
row.names(exp_df_hm) <- paste(paste(exp_df_hm$gene, exp_df_hm$celltype, sep = '_'), exp_df_hm$group, sep = '_')
exp_df_hm$group <- factor(exp_df_hm$group, levels = c('OVA5D', 'OVA5W', 'OVA3M', 'OVA6M', 'OVA9M', 'OVA12M', 'OVA15M'))

exp_df_hm$'Percent Expressed' <- a_df[row.names(exp_df_hm), 'pct.exp']

exp_df_hm$id <- paste(exp_df_hm$celltype, exp_df_hm$group, sep = '_')
exp_df_hm$celltype <- factor(exp_df_hm$celltype, levels = c('EC', 'Epi', 'GC', 'OO', 'IC', 'Mur', 'T&S'))

cols <- RColorBrewer::brewer.pal(n = 9, name = 'RdPu')

g <- ggplot(exp_df_hm, aes(x = group, y = gene, 
                           fill = mean_expression))  + geom_tile(color = 'white') + #+ geom_point(pch = 16, size = 3) + 
  facet_wrap(.~celltype, nrow = 1, scales = 'free_x') +
  theme_linedraw() + theme(panel.grid.major = element_blank()) + 
  theme(axis.text.x = element_text(angle = 45, hjust = NULL, 
                                   vjust = 0.5), axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
g <- g + scale_fill_gradientn(colors = alpha(colorRampPalette(cols)(99), 1)) +  guides(fill = guide_colourbar(barwidth = 0.8, 
                                                                                                              title = "Avg. exp."))
g <- g + theme(text = element_text(size = 13)) + 
  theme(legend.title = element_text(size = 11), legend.text = element_text(size = 11),
        strip.text = element_text(color = "black"), strip.background = element_rect(fill = "white"))

