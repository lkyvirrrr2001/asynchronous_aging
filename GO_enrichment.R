library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(Hmisc)
library(ggpubr)

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