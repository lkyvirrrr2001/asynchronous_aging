library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(ggVennDiagram)
library(UpSetR)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(Hmisc)
library(RColorBrewer)
library(ggpubr)
library(networkD3)
setwd("/venn_data/")
celltypes <- c('EC', 'Epi', 'GC', 'IC', 'Mur', 'T&S')
deg_df <- read.csv('OVA15MvsOVA3M.xls')  ##read the results "FindMarker()" between 15M and 3M in OVA
deg_df <- subset(deg_df, (avg_logFC>0.25|avg_logFC < -0.25) & p_val_adj<0.05)
deg_df$direction <- 'Up'
deg_df$direction[deg_df$avg_logFC < 0] <- 'Down'
up_df_ova <- subset(deg_df, direction == 'Up')
dn_df_ova <- subset(deg_df, direction == 'Down')
up_df_ova$celltype[up_df_ova$celltype == 'T&S'] <- 'SC/T&S'
dn_df_ova$celltype[dn_df_ova$celltype == 'T&S'] <- 'SC/T&S'

deg_df <- read.csv('OVI15MvsOVA3M.xls')  ##read the results "FindMarker()" between 15M and 3M in OVI
deg_df <- subset(deg_df, (avg_logFC>0.25|avg_logFC < -0.25) & p_val_adj<0.05)
deg_df$direction <- 'Up'
deg_df$direction[deg_df$avg_logFC < 0] <- 'Down'

up_df_ovi <- subset(deg_df, direction == 'Up')
dn_df_ovi <- subset(deg_df, direction == 'Down')
up_df_ovi$celltype[up_df_ovi$celltype == 'SC'] <- 'SC/T&S'
dn_df_ovi$celltype[dn_df_ovi$celltype == 'SC'] <- 'SC/T&S'

deg_df <- read.csv('UTE15MvsOVA3M.xls')  ##read the results "FindMarker()" between 15M and 3M
deg_df <- subset(deg_df, (avg_logFC>0.25|avg_logFC < -0.25) & p_val_adj<0.05)
deg_df$direction <- 'Up'
deg_df$direction[deg_df$avg_logFC < 0] <- 'Down'

up_df_ute <- subset(deg_df, direction == 'Up')
dn_df_ute <- subset(deg_df, direction == 'Down')
up_df_ute$celltype[up_df_ute$celltype == 'SC'] <- 'SC/T&S'
dn_df_ute$celltype[dn_df_ute$celltype == 'SC'] <- 'SC/T&S'

up_df_ova$organ <- "OVA"
up_df_ovi$organ <- "OVI"
up_df_ute$organ <- "UTE"
dn_df_ova$organ <- "OVA"
dn_df_ovi$organ <- "OVI"
dn_df_ute$organ <- "UTE"
up_df <- rbind(rbind(up_df_ova, up_df_ovi), up_df_ute)
up_df$celltype[up_df$celltype == "EPi"] <- "Epi"
dn_df <- rbind(rbind(dn_df_ova, dn_df_ovi), dn_df_ute)
dn_df$celltype[dn_df$celltype == "EPi"] <- "Epi"

sankey_plot <- function(up_df){
  tb <- as.data.frame(table(up_df$celltype, up_df$organ), stringsAsFactors = F)
  tb <- tb[tb$Freq != 0,]
  node <- data.frame(name = unique(c("OVA", "OVI", "UTE", tb$Var1)))
  #node$group <- as.character(1:14)
  node$group <- node$name
  link <- tb
  colnames(link) <- c("target", "source", "value")
  link$IDsource = match(link$source, node$name)-1
  link$IDtarget = match(link$target, node$name)-1
  link$link_group <- link$source
  
  link_up <- data.frame()
  #link_dn <- data.frame()
  common_genes <- list()
  OVA_specific <- list()
  OVI_specific <- list()
  UTE_specific <- list()
  for(ct in c("EC", "Epi", "IC", "SC/T&S")){
    up_df2 <- subset(up_df, celltype == ct)
    up_df3 <- as.data.frame(table(up_df2$X), stringsAsFactors = F)
    common_num <- nrow(up_df3[up_df3$Freq == 3,])
    common_df <- data.frame('Var1' = "", Freq = common_num*3, celltype = ct, type = "common")
    tmp_gene <- up_df3[up_df3$Freq == 3, "Var1"]
    common_genes[[ct]] <- tmp_gene
    
    tmp <- up_df3[up_df3$Freq == 1, 'Var1']
    specific_df <- subset(up_df2, X %in% tmp)
    specific_df2 <- as.data.frame(table(specific_df$organ), stringsAsFactor = F)
    specific_df2$celltype <- ct
    specific_df2$type <- "specific"
    tmp <- rbind(common_df, specific_df2)
    link_up <- rbind(link_up, tmp)
    
    tmp_gene <- specific_df[specific_df$organ == "OVA", "X"]
    OVA_specific[[ct]] <- tmp_gene
    tmp_gene <- specific_df[specific_df$organ == "OVI", "X"]
    OVI_specific[[ct]] <- tmp_gene
    tmp_gene <- specific_df[specific_df$organ == "UTE", "X"]
    UTE_specific[[ct]] <- tmp_gene
    
  }
  
  for(ct in c("MC", "PC", "GC", "Mur", "SMC", "Myo")){
    tmp <-  subset(up_df, celltype == ct)
    for (organ in unique(tmp$organ)){
      tmp2 <- data.frame('Var1' = organ, Freq = nrow(tmp[tmp$organ == organ,]), celltype = ct, type = "specific")
      link_up <- rbind(link_up, tmp2)
      
      tmp_gene <- tmp[tmp$organ == organ, "X"]
      if(organ == "OVA"){
        OVA_specific[[ct]] <- tmp_gene
      }
      if(organ == "OVI"){
        OVI_specific[[ct]] <- tmp_gene
      }
      if(organ == "UTE"){
        UTE_specific[[ct]] <- tmp_gene
      }
    }
    print(ct)
  }
  
  link_up$deg <- paste(link_up$Var1, link_up$type, sep = "_")

  
  link_tmp <- link_up[, c("deg", "Freq", "celltype")]
  link_tmp$link_group <- link_tmp$celltype
  colnames(link_tmp) <- c("target", "value", "source", "link_group")
  
  node2 <- rbind(node, data.frame(name = unique(link_tmp$target), group = 'DEG'))
  
  link_tmp$IDsource = match(link_tmp$source, node2$name)-1
  link_tmp$IDtarget = match(link_tmp$target, node2$name)-1
 
  link2 <- rbind(link, link_tmp)
  return(list(node2, link2, common_genes, OVA_specific, OVI_specific, UTE_specific))
  
}
up <- sankey_plot(up_df)
node2 <- up[[1]]   # input of sankeyNetwork
link2 <- up[[2]]   # input of sankeyNetwork

####enrichment analysis for up-regulated DEGs
common_genes <- up[[3]]
OVA_specific <- up[[4]]
OVI_specific <- up[[5]]
UTE_specific <- up[[6]]
genelist <- c()
for (ct in names(common_genes)){
  genelist <- union(genelist, common_genes[[ct]])
  genelist <- unique(genelist)
}
hg<-bitr(genelist,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Mm.eg.db")
go <- enrichGO(hg$ENTREZID,OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
go_df <- go@result
go_df <- go_df[order(go_df$p.adjust, decreasing = F),]
go_df <- subset(go_df, p.adjust < 0.0505)
go_df$classification <- "common"
go_df2 <- go_df

## file containing GO terms that you want to remove in the results.
protein_related <- read.csv("ribosomal_related_pathways.txt", sep = '\t', stringsAsFactors = F)[,1]

for (lst in c("OVA_specific", "OVI_specific", "UTE_specific")){
  l <- get(lst)
  for (name in names(l)){
    genelist <- l[[name]]
    hg<-bitr(genelist,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Mm.eg.db")
    go <- enrichGO(hg$ENTREZID,OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    go_df <- go@result
    go_df <- go_df[order(go_df$p.adjust, decreasing = F),]
    go_df <- subset(go_df, p.adjust < 0.0505)
    go_df$classification <- paste(lst, name, sep = '_')
    go_df2 <- rbind(go_df2, go_df)
  }
}
go_df2 <- go_df2[!(capitalize(as.character(go_df2$Description)) %in% protein_related),]
write.table(go_df2, "sankey_go_up.txt", sep = '\t', quote = F, row.names = F)

####enrichment analysis for down-regulated DEGs
up <- sankey_plot(dn_df)
sankey <- up[[1]]
common_genes <- up[[2]]
OVA_specific <- up[[3]]
OVI_specific <- up[[4]]
UTE_specific <- up[[5]]

####enrichment

genelist <- c()
for (ct in names(common_genes)){
  genelist <- union(genelist, common_genes[[ct]])
  genelist <- unique(genelist)
}
hg<-bitr(genelist,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Mm.eg.db")
go <- enrichGO(hg$ENTREZID,OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
go_df <- go@result
go_df <- go_df[order(go_df$p.adjust, decreasing = F),]
go_df <- subset(go_df, p.adjust < 0.0505)
go_df$classification <- "common"
go_df2 <- go_df

protein_related <- read.csv("ribosomal_related_pathways.txt", sep = '\t', stringsAsFactors = F)[,1]

for (lst in c("OVA_specific", "OVI_specific", "UTE_specific")){
  l <- get(lst)
  for (name in names(l)){
    genelist <- l[[name]]
    hg<-bitr(genelist,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Mm.eg.db")
    go <- enrichGO(hg$ENTREZID,OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    go_df <- go@result
    go_df <- go_df[order(go_df$p.adjust, decreasing = F),]
    go_df <- subset(go_df, p.adjust < 0.0505)
    go_df$classification <- paste(lst, name, sep = '_')
    go_df2 <- rbind(go_df2, go_df)
  }
}
go_df2 <- go_df2[!(capitalize(as.character(go_df2$Description)) %in% protein_related),]
write.table(go_df2, "sankey_go_dn.txt", sep = '\t', quote = F, row.names = F)

##Venn diagram
venn_list1 <- list("OVA" = unique(up_df[up_df$organ == "OVA", "X"]), "OVI" = unique(up_df[up_df$organ == "OVI", "X"]), 
                   "UTE" = unique(up_df[up_df$organ == "UTE", "X"]))
venn_list2 <- list("OVA" = unique(dn_df[dn_df$organ == "OVA", "X"]), "OVI" = unique(dn_df[dn_df$organ == "OVI", "X"]), 
                   "UTE" = unique(dn_df[dn_df$organ == "UTE", "X"]))
palette <- c("#FC8D62", "#66C2A5", "#8DA0CB") #set2 color

upset(fromList(venn_list1),
      sets = c("OVA", "OVI", "UTE"),
      keep.order = T,
      point.size = 2.5,
      line.size = 1,
      mainbar.y.label = 'Numbers of intersections of DEGs',
      sets.x.label = 'Number of DEGs',
      sets.bar.color = palette,
      main.bar.color = "#FF410DFF", #竖直柱图颜色
      matrix.color = "#FF410DFF", # 点和线的颜色
      text.scale = c(1.5,1.5,1.3,1.5,1.5,1.5)

)

upset(fromList(venn_list2),
      sets = c("OVA", "OVI", "UTE"),
      keep.order = T,
      point.size = 2.5,
      line.size = 1,
      mainbar.y.label = 'Numbers of intersections of DEGs',
      sets.x.label = 'Number of DEGs',
      sets.bar.color = palette,
      main.bar.color = "#5050FFFF", #竖直柱图颜色
      matrix.color = "#5050FFFF", # 点和线的颜色
      text.scale = c(1.5,1.5,1.3,1.5,1.5,1.5)
)

###volcano plot
up_df2 <- ddply(up_df, .(organ),reframe, DEG = unique(X))
gene.freq <- as.data.frame(table(up_df2$DEG), stringsAsFactors = F)
gene.freq <- gene.freq[order(gene.freq$Freq, decreasing = T),]
common_genes <- gene.freq[gene.freq$Freq == 3, 'Var1']
up_df2 <- subset(up_df, X %in% common_genes)
up_df2$X <- factor(up_df2$X, levels = unique(up_df2$X))
up_df3 <- ddply(up_df2, .(X), summarize, avg_logFC = mean(avg_logFC), p_val_adj = mean(p_val_adj))
up_df3$significant <- 'Up'

dn_df2 <- ddply(dn_df, .(organ),reframe, DEG = unique(X))
gene.freq <- as.data.frame(table(dn_df2$DEG), stringsAsFactors = F)
gene.freq <- gene.freq[order(gene.freq$Freq, decreasing = T),]
common_genes <- gene.freq[gene.freq$Freq == 3, 'Var1']
dn_df2 <- subset(dn_df, X %in% common_genes)
dn_df2$X <- factor(dn_df2$X, levels = unique(dn_df2$X))
dn_df3 <- ddply(dn_df2, .(X), summarize, avg_logFC = mean(avg_logFC), p_val_adj = mean(p_val_adj))
dn_df3$significant <- 'Down'

df_volcano <- rbind(up_df3, dn_df3)
(p <- ggplot(
  df_volcano, aes(x = avg_logFC, y = -log10(p_val_adj))) +
    geom_point(aes(color = significant), size=2) +
    scale_color_manual(values = c("Down" = "#5050FFFF","Up" = "#FF410DFF", "Stable" = "grey")))

