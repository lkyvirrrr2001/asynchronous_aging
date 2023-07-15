library(ggplot2)
library(ggpubr)
library(cowplot)
setwd("./6enrich_diffgenes/")
picked_go <- read.table('picked_GO.txt', sep = '\t', stringsAsFactors = F)
picked_kegg <- read.table('picked_KEGG.txt', sep = '\t', stringsAsFactors = F)
filter_updn <- function(ec_up_go, ec_dn_go){
  row.names(ec_dn_go) <- ec_dn_go$ID
  row.names(ec_up_go) <- ec_up_go$ID
  tmp_df1 <- subset(ec_up_go, !(ID %in% ec_dn_go$ID))
  overlap_go <- intersect(ec_up_go$ID, ec_dn_go$ID)
  overlap_df <- ec_up_go[overlap_go,]
  if(length(overlap_go) != 0){
    tmp_df2 <- overlap_df[ec_up_go[overlap_go,'p.adjust'] > ec_dn_go[overlap_go,'p.adjust'],]
    if(nrow(tmp_df2)!=0){
      df <- rbind(tmp_df1, tmp_df2)
    }else{
      df <- tmp_df1
    }
  }else{
    df <- tmp_df1
  }

  return(df)
}
top_enriched <- function(df){
  if(nrow(df)>5){
    df2 <- df[order(df$p.adjust),][1:5,]
    return(df2)
  }else{
    return(df)
  }
}
