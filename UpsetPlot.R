library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(UpSetR)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(Hmisc)
library(ggpubr)
setwd("D:/limo/我的坚果云/李墨修改的/venn_data/")
filepath1 <- 'D:/ThesisofAging/公司回传的/alldiffgene2/alldiffgene/OVA'
celltypes <- c('EC', 'Epi', 'GC', 'IC', 'Mur', 'T&S')

deg_df <- data.frame()
for (ct in celltypes){
  tmpname <- paste('/OVA15MvsOVA3M/marker_df_OVA15MvsOVA3M',paste(ct,'.xls', sep = ''), sep = '_')
  filename <- paste(filepath1, paste(ct, tmpname, sep = ''), sep = '/')
  dt <- read.csv(filename, sep = '\t', stringsAsFactors = F)
  dt2 <- subset(dt, (avg_logFC>0.25|avg_logFC < -0.25) & p_val_adj<0.05)
  dt2$celltype <- ct
  dt2$direction <- 'Up'
  dt2$direction[dt2$avg_logFC < 0] <- 'Down'
  deg_df <- rbind(deg_df, dt2)
}

up_df_ova <- subset(deg_df, direction == 'Up')
dn_df_ova <- subset(deg_df, direction == 'Down')
up_df_ova$celltype[up_df_ova$celltype == 'T&S'] <- 'SC/T&S'
dn_df_ova$celltype[dn_df_ova$celltype == 'T&S'] <- 'SC/T&S'

filepath1 <- 'D:/ThesisofAging/公司回传的/alldiffgene2/alldiffgene/OVI'
celltypes <- c('EC', 'EPi', 'MC', 'IC', 'PC', 'SMC', 'SC')

deg_df <- data.frame()
for (ct in celltypes){
  tmpname <- paste('/OVI15MvsOVI3M/marker_df_OVI15MvsOVI3M',paste(ct,'.xls', sep = ''), sep = '_')
  filename <- paste(filepath1, paste(ct, tmpname, sep = ''), sep = '/')
  dt <- read.csv(filename, sep = '\t', stringsAsFactors = F)
  dt2 <- subset(dt, (avg_logFC>0.25|avg_logFC < -0.25) & p_val_adj<0.05)
  dt2$celltype <- ct
  dt2$direction <- 'Up'
  dt2$direction[dt2$avg_logFC < 0] <- 'Down'
  deg_df <- rbind(deg_df, dt2)
}

up_df_ovi <- subset(deg_df, direction == 'Up')
dn_df_ovi <- subset(deg_df, direction == 'Down')
up_df_ovi$celltype[up_df_ovi$celltype == 'SC'] <- 'SC/T&S'
dn_df_ovi$celltype[dn_df_ovi$celltype == 'SC'] <- 'SC/T&S'

filepath1 <- 'D:/ThesisofAging/公司回传的/alldiffgene2/alldiffgene/UTE'
celltypes <- c('EC', 'EPi', 'MC', 'IC', 'PC', 'Myo', 'SC')

deg_df <- data.frame()
for (ct in celltypes){
  tmpname <- paste('/UTE15MvsUTE3M/marker_df_UTE15MvsUTE3M',paste(ct,'.xls', sep = ''), sep = '_')
  filename <- paste(filepath1, paste(ct, tmpname, sep = ''), sep = '/')
  dt <- read.csv(filename, sep = '\t', stringsAsFactors = F)
  dt2 <- subset(dt, (avg_logFC>0.25|avg_logFC < -0.25) & p_val_adj<0.05)
  dt2$celltype <- ct
  dt2$direction <- 'Up'
  dt2$direction[dt2$avg_logFC < 0] <- 'Down'
  deg_df <- rbind(deg_df, dt2)
}

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

# length(unique(up_df[up_df$organ == "OVA", "X"]))
# length(unique(dn_df[dn_df$organ == "OVA", "X"]))
# up_df$label <- paste(up_df$X, up_df$celltype, sep = '_')
# dn_df$label <- paste(dn_df$X, dn_df$celltype, sep = '_')

venn_list1 <- list("OVA" = unique(up_df[up_df$organ == "OVA", "X"]), "OVI" = unique(up_df[up_df$organ == "OVI", "X"]), 
                   "UTE" = unique(up_df[up_df$organ == "UTE", "X"]))
venn_list2 <- list("OVA" = unique(dn_df[dn_df$organ == "OVA", "X"]), "OVI" = unique(dn_df[dn_df$organ == "OVI", "X"]), 
                   "UTE" = unique(dn_df[dn_df$organ == "UTE", "X"]))
palette <- c("#FC8D62", "#66C2A5", "#8DA0CB") #set2 color

pdf('./3organs/upset_3organs_up.pdf', width = 5.5, height = 4.5)
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
      # queries = list(list(query = intersects, params = list("OVA", "OVI", "UTE"), color = "#D7191C", active = T, query.name = "Intersection of three organs")),
      #                 query.legend = "bottom"
)
dev.off()

pdf('./3organs/upset_3organs_dn.pdf', width = 5.5, height = 4.5)
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
      # queries = list(list(query = intersects, params = list("OVA", "OVI", "UTE"), color = "#D7191C", active = T, query.name = "Intersection of three organs")),
      #                 query.legend = "bottom"
)
dev.off()

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
  # 数据、映射、颜色
  df_volcano, aes(x = avg_logFC, y = -log10(p_val_adj))) +
    geom_point(aes(color = significant), size=2) +
    #   scale_y_continuous(limits = (c(-5,300)))+
    scale_color_manual(values = c("Down" = "#5050FFFF","Up" = "#FF410DFF", "Stable" = "grey")))
(p2 <- p+geom_text_repel(
  data = df_volcano,
  aes(label = X),
  linewidth = 5,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  label.size = 0.08,
  max.overlaps = 20) +
    # 辅助线
    geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = 1,lty=4,col="black",lwd=0.8) +
    # 坐标轴
    labs(x="Average log2(fold change)",
         y="Average -log10 (P.adj)") +
    # 图例
    theme_light() +
    theme( legend.title = element_blank(), #legend.position = "bottom",
           legend.text = element_text(size = 13),
           axis.text = element_text(color = 'black', size = 14), axis.title = element_text(color = 'black', size = 14))
)

pdf('./3organs/commonDEG_volcano.pdf', width = 7, height = 5)
p2
dev.off()

enrichment_analysis <- function(genelist){
  hg<-bitr(genelist,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Mm.eg.db")
  go <- enrichGO(hg$ENTREZID,OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go_df <- go@result
  go_df <- go_df[order(go_df$p.adjust, decreasing = F),]
  go_df <- subset(go_df, p.adjust < 0.0505)
  go_df <- go_df[!(capitalize(as.character(go_df$Description)) %in% protein_related),]
  if(nrow(go_df) > 10){
    go_top <- go_df[1:10,]
  }else{
    go_top <- go_df
  }
  tmp <- strsplit(go_top$GeneRatio, split = '/')
  ratio <- unlist(lapply(tmp, function(x)as.numeric(x[1])/as.numeric(x[2])))
  go_top$GeneRatio <- round(ratio, 3)
  go_top <- go_top[order(go_top$p.adjust, decreasing = T),]
  go_top <- go_top[order(go_top$GeneRatio, decreasing = F),]
  #  go_top <- go_top[order(go_top$p.adjust, decreasing = F),]
  go_top$Description <- factor(capitalize(go_top$Description), levels = capitalize(go_top$Description))
  
  return(go_top)
}

enrichment_plot <- function(genelist, direction){
  go_top <- enrichment_analysis(genelist)
  if(direction == 'Up'){
    cols <- cols1
  }else{
    cols <- cols2
  }
  (p <- ggdotchart(go_top, y = 'Description', x = 'GeneRatio', color = "p.adjust", dot.size = "Count" , size = 1, ylab = '', sorting = 'none', legend = 'top') + theme_light()
    + theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 10,angle = 45, hjust = 1), legend.text = element_text(size = 12), 
            legend.title = element_text(size = 12))
    + scale_color_gradientn(colors = rev(cols))) 
  return(p)
}
cols1 <- get_palette('YlOrRd',10)
cols2 <- brewer.pal(9, 'YlGnBu')
protein_related <- read.csv("ribosomal_related_pathways.txt", sep = '\t', stringsAsFactors = F)[,1]
common_up <- unique(up_df2$X)
common_dn <- unique(dn_df2$X)
p1 <- enrichment_plot(common_up, 'Up')
p2 <- enrichment_plot(common_dn, 'Down')
pdf('./3organs/common_up_enrich.pdf', width = 12.5, height = 4.5)
p1
dev.off()
pdf('./3organs/common_dn_enrich.pdf', width = 6.45, height = 4.2)
p2
dev.off()

