library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(ggVennDiagram)
library(UpSetR)
library(ggrepel)
setwd("./venn_data/")

merge_array <- function(ar1, ar2){
  genes1 <- row.names(ar1)
  genes2 <- row.names(ar2)
  #  union_genes <- union(genes1, genes2)
  common_genes <- intersect(genes1, genes2)
  diff_genes1 <- setdiff(genes1, genes2)
  diff_genes2 <- setdiff(genes2, genes1)
  
  ar <- cbind(ar1[common_genes,], ar2[common_genes,])
  ar_tmp1 <- cbind(ar1[diff_genes1,],matrix(0, nrow = length(diff_genes1), ncol = ncol(ar2), dimnames = list(diff_genes1, colnames(ar2))))
  ar_tmp2 <- cbind(matrix(0, nrow = length(diff_genes2), ncol = ncol(ar1), dimnames = list(diff_genes2, colnames(ar1))),ar2[diff_genes2,])
  ar <- rbind(ar,ar_tmp1)
  ar <- rbind(ar,ar_tmp2)
  return(ar)
}

make_df <- function(table, geneset){
  df_melt0 <- data.frame()
  i = 1
  for (ctype in unique(table$celltype)) {
    
    deg_array <- acast(table, gene_id~compare, subset = .(celltype==ctype))
    deg_array[!is.na(deg_array)] <- 1
    deg_array[is.na(deg_array)] <- 0
    
    if (i == 1) {
      deg_array0 <- deg_array
    } else{
      deg_array0 <- merge_array(deg_array0,deg_array) 
    }
    
    tmp_gene <- setdiff(geneset, row.names(deg_array))
    df_tmp <- matrix(0, nrow=length(tmp_gene), ncol = ncol(deg_array), dimnames = list(tmp_gene, colnames(deg_array)))
    deg_array <- rbind(deg_array, df_tmp)
    df_melt <- melt(deg_array)
    names(df_melt) <- c('Gene', 'Compare', 'Exist')
    df_melt$celltype <- ctype
    df_melt0 <- rbind(df_melt0, df_melt)
    
    i = i + 1
  }
  df_melt0$Compare <- factor(df_melt0$Compare, levels = c("5Wvs5D", "3Mvs5W", "6Mvs3M", "9Mvs6M", "12Mvs9M", "15Mvs12M"))
  df_melt0$celltype <- factor(df_melt0$celltype, levels = unique(table$celltype))
  return(list(deg_array0 = deg_array0, df_melt0 = df_melt0))
}

order_genes <- function(df, arr){
  gene_set0 <- c()
  for (ctype in unique(df$celltype)) {
    genes <- unique(subset(df, celltype==ctype)$gene_id)
    gene_set0 <- c(gene_set0, genes)
  }
  df1 <-  as.data.frame(table(gene_set0))
  df1 <- df1[order(df1$Freq,decreasing = T),]
  #write.table(df1, table_name, sep = '\t', row.names = F, quote = F)
  #table(df1$Freq)
  ordered_genes <- c()
  gene_num <- c()
  s <- 0
  for(n in unique(df1$Freq)){
  
      # gene_set1 <- as.character(subset(df1, Freq %in% 2:6)$gene_set0)
      # sub_arr <- arr[gene_set1,]
      # d <- dist(sub_arr)
      # clust <-hclust(d)
      # ordered_genes <- c(ordered_genes, row.names(sub_arr)[clust$order])

      gene_set1 <- as.character(subset(df1, Freq == n)$gene_set0)
      sub_arr <- arr[gene_set1,]
      d <- dist(sub_arr)
      clust <-hclust(d)
      ordered_genes <- c(ordered_genes, row.names(sub_arr)[clust$order])
      s <- s + length(gene_set1)
      gene_num <- c(gene_num, s)

  }
  return(list(ordered_genes = ordered_genes, gene_num = gene_num))
}

####read the up-regulated genes (results of 'FindMarkers' function in Seurat package)
ova_up <- read.csv('P21062303_OVA_totalSampleCompare_diffgenes.up.xls', sep = '\t', stringsAsFactors = F)
ovi_up <- read.csv('P21062303_OVI_totalSampleCompare_diffgenes.rm.undefined.up.xls', sep = '\t', stringsAsFactors = F)
ute_up <- read.csv('P21062303_UTE_totalSampleCompare_diffgenes.up.xls', sep = '\t', stringsAsFactors = F)
cols <- c("0" = "#E8E9EB", "1" = "#FF410DFF")
########ova
gene_set <- unique(ova_up$gene_id)
tmp <- make_df(ova_up, gene_set)
df_melt0 <- tmp$df_melt0
deg_array0 <- tmp$deg_array0

# d <- dist(deg_array0)
# clust <-hclust(d)
gene_order <- order_genes(ova_up, deg_array0)$ordered_genes
df_melt0$Gene <- factor(df_melt0$Gene,level = gene_order)
num <- order_genes(ova_up, deg_array0)$gene_num
xaxis_lab <- gene_order[c(1,num)]
pdf('ova_up_genes_heatmap.pdf', height = 9.5)
(P <- ggplot(df_melt0,aes(Gene, Compare, fill = Exist)) + geom_tile() + scale_fill_manual(values = cols)
  + facet_wrap(~celltype, ncol = 1, strip.position = "right") + theme_minimal() + scale_x_discrete(breaks = xaxis_lab, labels = c(1,num))
  + theme( axis.title.y = element_blank(), legend.position="none", axis.text.x = element_text(size = 12), 
           axis.title = element_text(size = 13),  axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 13), strip.text = element_text(size = 13, angle = 90))
  + labs(x='Up-regulated genes') )
dev.off()

ova_up_inter <- subset(df_melt0, Gene %in% gene_order[1:num[1]])
########ovi
gene_set <- unique(ovi_up$gene_id)
tmp <- make_df(ovi_up, gene_set)
df_melt0 <- tmp$df_melt0
deg_array0 <- tmp$deg_array0

gene_order <- order_genes(ovi_up, deg_array0)$ordered_genes
df_melt0$Gene <- factor(df_melt0$Gene,level = gene_order)
num <- order_genes(ovi_up, deg_array0)$gene_num
xaxis_lab <- gene_order[c(1,num)]

pdf('ovi_up_genes_heatmap.pdf', height = 9.5)
(P <- ggplot(df_melt0,aes(Gene, Compare, fill = Exist)) + geom_tile() + scale_fill_manual(values = cols)
  + facet_wrap(~celltype, ncol = 1, strip.position = "right", scales = "free_y") + theme_minimal() + scale_x_discrete(breaks = xaxis_lab, labels = c(1,num))
  + theme( axis.title.y = element_blank(), legend.position="none", axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.title = element_text(size = 13), 
           axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 13), strip.text = element_text(size = 13)) + labs(x='Genes') )
dev.off()    

########ute
gene_set <- unique(ute_up$gene_id)
tmp <- make_df(ute_up, gene_set)
df_melt0 <- tmp$df_melt0
deg_array0 <- tmp$deg_array0

gene_order <- order_genes(ute_up, deg_array0)$ordered_genes
df_melt0$Gene <- factor(df_melt0$Gene,level = gene_order)
num <- order_genes(ute_up, deg_array0)$gene_num
xaxis_lab <- gene_order[c(1,num)]

pdf('ute_up_genes_heatmap.pdf', height = 9.5)
(P <- ggplot(df_melt0,aes(Gene, Compare, fill = Exist)) + geom_tile() + scale_fill_manual(values = cols)
  + facet_wrap(~celltype, ncol = 1, strip.position = "right", scales = "free_y") + theme_minimal() + scale_x_discrete(breaks = xaxis_lab, labels = c(1,num))
  + theme( axis.title.y = element_blank(), legend.position="none", axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.title = element_text(size = 13), 
           axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 13), strip.text = element_text(size = 13)) + labs(x='Genes') )

dev.off()    

####read the down-regulated genes (results of 'FindMarkers' function in Seurat package)
ova_dn <- read.csv('P21062303_OVA_totalSampleCompare_diffgenes.down.xls', sep = '\t')
ovi_dn <- read.csv('P21062303_OVI_totalSampleCompare_diffgenes.rm.undefined.down.xls', sep = '\t')
ute_dn <- read.csv('P21062303_UTE_totalSampleCompare_diffgenes.down.xls', sep = '\t')
cols <- c("0" = "#E8E9EB", "1" = "#5050FFFF") #show_col(pal_tron()(7))
########ova
gene_set <- unique(ova_dn$gene_id)
tmp <- make_df(ova_dn, gene_set)
df_melt0 <- tmp$df_melt0
deg_array0 <- tmp$deg_array0

gene_order <- order_genes(ova_dn, deg_array0)$ordered_genes
df_melt0$Gene <- factor(df_melt0$Gene,level = gene_order)
num <- order_genes(ova_dn, deg_array0)$gene_num
xaxis_lab <- gene_order[c(1,num)]

pdf('ova_dn_genes_heatmap.pdf', height = 9.5)
(P <- ggplot(df_melt0,aes(Gene, Compare, fill = Exist)) + geom_tile() + scale_fill_manual(values = cols)
  + facet_wrap(~celltype, ncol = 1, strip.position = "right") + theme_minimal() + scale_x_discrete(breaks = xaxis_lab, labels = c(1,num))
  + theme( axis.title.y = element_blank(), legend.position="none", axis.text.x = element_text(size = 12), 
           axis.title = element_text(size = 13),  axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 13), strip.text = element_text(size = 13, angle = 90))
  + labs(x='Down-regulated genes'))
dev.off()  

########ovi
gene_set <- unique(ovi_dn$gene_id)
tmp <- make_df(ovi_dn, gene_set)
df_melt0 <- tmp$df_melt0
deg_array0 <- tmp$deg_array0

gene_order <- order_genes(ovi_dn, deg_array0)$ordered_genes
df_melt0$Gene <- factor(df_melt0$Gene,level = gene_order)
num <- order_genes(ovi_dn, deg_array0)$gene_num
xaxis_lab <- gene_order[c(1,num)]


pdf('ovi_dn_genes_heatmap.pdf', height = 9.5)
(P <- ggplot(df_melt0,aes(Gene, Compare, fill = Exist)) + geom_tile() + scale_fill_manual(values = cols)
  + facet_wrap(~celltype, ncol = 1, strip.position = "right", scales = "free_y") + theme_minimal() + scale_x_discrete(breaks = xaxis_lab, labels = c(1,num))
  + theme( axis.title.y = element_blank(), legend.position="none", axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.title = element_text(size = 13), 
           axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 13), strip.text = element_text(size = 13)) + labs(x='Genes') )
dev.off()    

########ute
gene_set <- unique(ute_dn$gene_id)
tmp <- make_df(ute_dn, gene_set)
df_melt0 <- tmp$df_melt0
deg_array0 <- tmp$deg_array0

gene_order <- order_genes(ute_dn, deg_array0)$ordered_genes
df_melt0$Gene <- factor(df_melt0$Gene,level = gene_order)
num <- order_genes(ute_dn, deg_array0)$gene_num
xaxis_lab <- gene_order[c(1,num)]

pdf('ute_dn_genes_heatmap.pdf', height = 9.5)
(P <- ggplot(df_melt0,aes(Gene, Compare, fill = Exist)) + geom_tile() + scale_fill_manual(values = cols)
  + facet_wrap(~celltype, ncol = 1, strip.position = "right", scales = "free_y") + theme_minimal() + scale_x_discrete(breaks = xaxis_lab, labels = c(1,num))
  + theme( axis.title.y = element_blank(), legend.position="none", axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.title = element_text(size = 13), 
           axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 13), strip.text = element_text(size = 13)) + labs(x='Genes') )
dev.off()    
