##this script was used to generate node and edge tables to be input for Cytoscape.
setwd('D:/limo/我的坚果云/3organs/tf_analysis/')
library(plyr)
###read TFs and regulons output by pySCENIC:
tf_df1 <- read.csv('tf_gene/OVA_add/5W5D_OVA.tf_add.csv')[-1,-1]
tf_df2 <- read.csv('tf_gene/OVA_add/3M5W_OVA.tf_add.csv')[-1,-1]
tf_df3 <- read.csv('tf_gene/OVA_add/6M3M_OVA.tf_add.csv')[-1,-1]
tf_df4 <- read.csv('tf_gene/OVA_add/9M6M_OVA.tf_add.csv')[-1,-1]
tf_df5 <- read.csv('tf_gene/OVA_add/12M9M_OVA.tf_add.csv')[-1,-1]
tf_df6 <- read.csv('tf_gene/OVA_add/15M12M_OVA.tf_add.csv')[-1,-1]
tf_df7 <- read.csv('tf_gene/OVA_add/15M3M_OVA.tf_add.csv')[-1,-1]

sep_lines <- function(tmp_df){
  genes <- unlist(strsplit(tmp_df$gene, split = ',', fixed = T))
  tmp_df2 <- data.frame( gene = genes)#comparison = tmp_df$comparison,
  return(tmp_df2)
}

add_color <- function(tmp_df){
  if(tmp_df$type=='TFs'&tmp_df$regulation=='up'){
    tmp_df$color_lab <- 'TFs-up'
  }else if(tmp_df$type=='TFs'&tmp_df$regulation=='down'){
    tmp_df$color_lab <- 'TFs-down'
  }else{
    tmp_df$color_lab <- tmp_df$celltype
  }
  return(tmp_df)
}

#select TFs that significantly differentially expressed between two groups.
gene_diff <- function(sub_tf_df_sep, gene_df_dn, gene_df_up, group = '5Wvs5D'){
  sub_gene_df_dn <- subset(gene_df_dn, compare==group)
  sub_gene_df_dn$direction <- 'down'
  sub_gene_df_up <- subset(gene_df_up, compare==group)
  sub_gene_df_up$direction <- 'up'
  sub_gene_df <- rbind(sub_gene_df_dn, sub_gene_df_up)
  row.names(sub_gene_df) <- paste(sub_gene_df$celltype, sub_gene_df$gene_id, sep = '&')
  query <- paste(sub_tf_df_sep$celltype, sub_tf_df_sep$gene, sep = '&')
  direction <- sub_gene_df[query,'direction']
  sub_tf_df_sep2 <- sub_tf_df_sep[!is.na(direction),]
  sub_tf_df_sep2$direction <- direction[!is.na(direction)]
  return(sub_tf_df_sep2)
}

make_file <- function(tf_df1, gene_df_dn, gene_df_up, compare, prefix){
  row.names(tf_df1) <- tf_df1$regulon
  tf_df1_tmp <- subset(tf_df1,gene!='NA,'&gene!='0')
  tf_df1_tmp$regulation <- 'up'
  tf_df1_tmp$regulation[tf_df1_tmp$avg_logFC<0] <- 'down'
  
  tf_df_sep <- adply(tf_df1_tmp,1,sep_lines)

  sub_tf_df_sep <- gene_diff(tf_df_sep, gene_df_dn, gene_df_up, group = compare)
  
  a <- as.data.frame(table(sub_tf_df_sep$regulon), stringsAsFactors = F)
  a <- a[order(a$Freq, decreasing = T),]
  row.names(a) <- a$Var1
  
#select TFs with more than 30 target genes:
  if(nrow(a)>30){
    selected_TFs <- a$Var1[1:30]
  }else{
    selected_TFs <- a$Var1
  }
  
  sub_tf_df_sep2 <- subset(sub_tf_df_sep, regulon %in% selected_TFs)
  sub_tf_df_sep2$gene_lab <- paste('gene',1:nrow(sub_tf_df_sep2), sep = '')
  
  edge_file <- paste(prefix, 'edge.txt', sep = '_')
  node_file <- paste(prefix, 'node.txt', sep = '_')
  write.table(sub_tf_df_sep2, edge_file, sep = '\t', quote = F, row.names = F)
  node_df <- data.frame(node = selected_TFs, type = 'TFs', regulation = tf_df1_tmp[selected_TFs,]$regulation, celltype = tf_df1_tmp[selected_TFs,]$celltype, gene_number = a[selected_TFs,'Freq'])
  
  tmp_df <- data.frame(node = sub_tf_df_sep2$gene_lab, type = 'genes', regulation = sub_tf_df_sep2$direction, celltype = sub_tf_df_sep2$celltype, gene_number = 1)
  node_df <- rbind(node_df, tmp_df)
  node_df2 <- adply(node_df,1,add_color)
  write.table(node_df2,node_file, sep = '\t', quote = F, row.names = F)
  
}

gene_df_dn <- read.csv('P21062303_OVA_totalSampleCompare_diffgenes.down.xls', sep = '\t')
gene_df_up <- read.csv('P21062303_OVA_totalSampleCompare_diffgenes.up.xls', sep = '\t')
make_file(tf_df1, gene_df_dn, gene_df_up, '5Wvs5D', prefix = '5W5D')
make_file(tf_df2, gene_df_dn, gene_df_up, '3Mvs5W', prefix = '3M5W')
make_file(tf_df3, gene_df_dn, gene_df_up, '6Mvs3M', prefix = '6M3M')
make_file(tf_df4, gene_df_dn, gene_df_up, '9Mvs6M', prefix = '9M6M')
make_file(tf_df5, gene_df_dn, gene_df_up, '12Mvs9M', prefix = '12M9M')
make_file(tf_df6, gene_df_dn, gene_df_up, '15Mvs12M', prefix = '15M12M')
diff_genes_add_15Mvs3M <- read.csv('./P21062303_OVA_totalSampleCompare_diffgenes.xls', sep = '\t')
gene_df_up <- subset(diff_genes_add_15Mvs3M,avg_logFC>0)
gene_df_dn <- subset(diff_genes_add_15Mvs3M,avg_logFC<0)
make_file(tf_df7, gene_df_dn, gene_df_up, '15Mvs3M', prefix = '15M3M')

###OVI
tf_df1 <- read.csv('tf_gene/OVI_add/5W5D_OVI.tf_add.csv')[-1,-1]
tf_df2 <- read.csv('tf_gene/OVI_add/3M5W_OVI.tf_add.csv')[-1,-1]
tf_df3 <- read.csv('tf_gene/OVI_add/6M3M_OVI.tf_add.csv')[-1,-1]
tf_df4 <- read.csv('tf_gene/OVI_add/9M6M_OVI.tf_add.csv')[-1,-1]
tf_df5 <- read.csv('tf_gene/OVI_add/12M9M_OVI.tf_add.csv')[-1,-1]
tf_df6 <- read.csv('tf_gene/OVI_add/15M12M_OVI.tf_add.csv')[-1,-1]
tf_df7 <- read.csv('tf_gene/OVI_add/15M3M_OVI.tf_add.csv')[-1,-1]
gene_df_dn <- read.csv('P21062303_OVI_totalSampleCompare_diffgenes.rm.undefined.down.xls', sep = '\t')
gene_df_up <- read.csv('P21062303_OVI_totalSampleCompare_diffgenes.rm.undefined.up.xls', sep = '\t')
make_file(tf_df1, gene_df_dn, gene_df_up, '5Wvs5D', prefix = 'OVI_5W5D')
make_file(tf_df2, gene_df_dn, gene_df_up, '3Mvs5W', prefix = 'OVI_3M5W')
make_file(tf_df3, gene_df_dn, gene_df_up, '6Mvs3M', prefix = 'OVI_6M3M')
make_file(tf_df4, gene_df_dn, gene_df_up, '9Mvs6M', prefix = 'OVI_9M6M')
make_file(tf_df5, gene_df_dn, gene_df_up, '12Mvs9M', prefix = 'OVI_12M9M')
make_file(tf_df6, gene_df_dn, gene_df_up, '15Mvs12M', prefix = 'OVI_15M12M')
diff_genes_add_15Mvs3M <- read.csv('P21062303_OVI_totalSampleCompare_diffgenes.xls', sep = '\t')
gene_df_up <- subset(diff_genes_add_15Mvs3M,avg_logFC>0)
gene_df_dn <- subset(diff_genes_add_15Mvs3M,avg_logFC<0)
make_file(tf_df7, gene_df_dn, gene_df_up, '15Mvs3M', prefix = 'OVI_15M3M')

###UTE
tf_df1 <- read.csv('tf_gene/UTE_add/5W5D_UTE.tf_add.csv')[-1,-1]
tf_df2 <- read.csv('tf_gene/UTE_add/3M5W_UTE.tf_add.csv')[-1,-1]
tf_df3 <- read.csv('tf_gene/UTE_add/6M3M_UTE.tf_add.csv')[-1,-1]
tf_df4 <- read.csv('tf_gene/UTE_add/9M6M_UTE.tf_add.csv')[-1,-1]
tf_df5 <- read.csv('tf_gene/UTE_add/12M9M_UTE.tf_add.csv')[-1,-1]
tf_df6 <- read.csv('tf_gene/UTE_add/15M12M_UTE.tf_add.csv')[-1,-1]
tf_df7 <- read.csv('tf_gene/UTE_add/15M3M_UTE.tf_add.csv')[-1,-1]
gene_df_dn <- read.csv('P21062303_UTE_totalSampleCompare_diffgenes.down.xls', sep = '\t')
gene_df_up <- read.csv('P21062303_UTE_totalSampleCompare_diffgenes.up.xls', sep = '\t')
make_file(tf_df1, gene_df_dn, gene_df_up, '5Wvs5D', prefix = 'UTE_5W5D')
make_file(tf_df2, gene_df_dn, gene_df_up, '3Mvs5W', prefix = 'UTE_3M5W')
make_file(tf_df3, gene_df_dn, gene_df_up, '6Mvs3M', prefix = 'UTE_6M3M')
make_file(tf_df4, gene_df_dn, gene_df_up, '9Mvs6M', prefix = 'UTE_9M6M')
make_file(tf_df5, gene_df_dn, gene_df_up, '12Mvs9M', prefix = 'UTE_12M9M')
make_file(tf_df6, gene_df_dn, gene_df_up, '15Mvs12M', prefix = 'UTE_15M12M')
diff_genes_add_15Mvs3M <- read.csv('P21062303_UTE_totalSampleCompare_diffgenes.xls', sep = '\t')
gene_df_up <- subset(diff_genes_add_15Mvs3M,avg_logFC>0)
gene_df_dn <- subset(diff_genes_add_15Mvs3M,avg_logFC<0)
make_file(tf_df7, gene_df_dn, gene_df_up, '15Mvs3M', prefix = 'UTE_15M3M')

########select TFs OVA for visualization in heatmaps
node1 <- read.csv('5W5D_node.txt', sep = '\t')
tf1 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('3M5W_node.txt', sep = '\t')
tf2 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('6M3M_node.txt', sep = '\t')
tf3 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('9M6M_node.txt', sep = '\t')
tf4 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('12M9M_node.txt', sep = '\t')
tf5 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('15M12M_node.txt', sep = '\t')
tf6 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('15M3M_node.txt', sep = '\t')
tf7 <- node1[node1$type == 'TFs','node']
unique_tf <- unique(c(tf1,tf2,tf3,tf4,tf5,tf6,tf7))
tb_ova <- as.data.frame(table(c(tf1,tf2,tf3,tf4,tf5,tf6,tf7)))
tb_ova <- tb_ova[order(tb_ova$Freq, decreasing = T),]
a <- tb_ova[tb_ova$Freq>=4,]
write.table(a, 'ova_tf_freq.txt', sep = '\t', row.names = F, quote = F)

########select TFs OVI for visualization in heatmaps
node1 <- read.csv('OVI_5W5D_node.txt', sep = '\t')
tf1 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('OVI_3M5W_node.txt', sep = '\t')
tf2 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('OVI_6M3M_node.txt', sep = '\t')
tf3 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('OVI_9M6M_node.txt', sep = '\t')
tf4 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('OVI_12M9M_node.txt', sep = '\t')
tf5 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('OVI_15M12M_node.txt', sep = '\t')
tf6 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('OVI_15M3M_node.txt', sep = '\t')
tf7 <- node1[node1$type == 'TFs','node']
#unique_tf <- unique(c(tf1,tf2,tf3,tf4,tf5,tf6,tf7))
tb_ovi <- as.data.frame(table(c(tf1,tf2,tf3,tf4,tf5,tf6,tf7)))
tb_ovi <- tb_ovi[order(tb_ovi$Freq, decreasing = T),]
a <- tb_ovi[tb_ovi$Freq>=4,]
write.table(a, 'ovi_tf_freq.txt', sep = '\t', row.names = F, quote = F)

########select TFs UTE for visualization in heatmaps
node1 <- read.csv('UTE_5W5D_node.txt', sep = '\t')
tf1 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('UTE_3M5W_node.txt', sep = '\t')
tf2 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('UTE_6M3M_node.txt', sep = '\t')
tf3 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('UTE_9M6M_node.txt', sep = '\t')
tf4 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('UTE_12M9M_node.txt', sep = '\t')
tf5 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('UTE_15M12M_node.txt', sep = '\t')
tf6 <- node1[node1$type == 'TFs','node']
node1 <- read.csv('UTE_15M3M_node.txt', sep = '\t')
tf7 <- node1[node1$type == 'TFs','node']
#unique_tf <- unique(c(tf1,tf2,tf3,tf4,tf5,tf6,tf7))
tb_ute <- as.data.frame(table(c(tf1,tf2,tf3,tf4,tf5,tf6,tf7)))
tb_ute <- tb_ute[order(tb_ute$Freq, decreasing = T),]
a <- tb_ute[tb_ute$Freq>=4,]
write.table(a, 'ute_tf_freq.txt', sep = '\t', row.names = F, quote = F)
