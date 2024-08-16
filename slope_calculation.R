library(ggplot2)
library(ggpubr)
library(ggsci)
library(rstatix)
library(cowplot)
library(plyr)
library(scales)
library(pheatmap)
library(reshape2)
library(wesanderson)
library(lsmeans)
setwd("./geneset_score/ucell_slope/")
#######calculate Ucell slope using linear regression model in each organ.
dir1_1 <- c("1.HALLMARK_DNA_REPAIR", "12.HALLMARK_APOPTOSIS", "4.REACTOME_CELLULAR_SENESCENCE", "5.REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
            "17.GOBP_DNA_METHYLATION", "7.HALLMARK_INFLAMMATORY_RESPONSE", "11.BIOCARTA_IGF1_PATHWAY", "26.GOBP_NEGATIVE_REGULATION_OF_CELL_CYCLE")
groups <- c('3M', '6M', '9M', '12M', '15M')

calculate_slope <- function(dir1_1){
  ucell_df1 <- data.frame()
  slope_df <- data.frame()
  signif_df <- data.frame() 
  for (dir0 in dir1_1){
    dir2 <- paste(paste("../pathgroup_0727/pathgroup/", dir0, sep = ''), 'merge', sep = '/')
    pathway <- strsplit(dir0, '\\.')[[1]][2]
    ucell_df1 <- data.frame()
    for(g in groups){
      dt <- paste(dir2, paste(g, '_data.xls', sep = ''), sep = '/')
      score_df <- read.csv(dt, sep = '\t')
      score_df$pathway <- pathway
      score_df <- score_df[,c('orig.ident', 'group', 'UCell', 'pathway')]
      ucell_df1 <- rbind(ucell_df1, score_df)
      
      
    }
    ucell_df1$age <- 3
    ucell_df1$age[ucell_df1$group == '3M'] <- 3
    ucell_df1$age[ucell_df1$group == '6M'] <- 6
    ucell_df1$age[ucell_df1$group == '9M'] <- 9
    ucell_df1$age[ucell_df1$group == '12M'] <- 12
    ucell_df1$age[ucell_df1$group == '15M'] <- 15
    ucell_df1$orig.ident <- factor(ucell_df1$orig.ident, levels = c('OVA', 'OVI', 'UTE'))
    lm_ucell <- lm(UCell ~ age+orig.ident+age*orig.ident , data = ucell_df1)

    lst_ucell <- lstrends(lm_ucell, 'orig.ident', var = 'age')
    pair_df <- as.data.frame(pairs(lst_ucell))
    pair_df$pathway <- pathway
    signif_df <- rbind(signif_df, pair_df)
    lst_ucell <- as.data.frame(lst_ucell)
    lst_ucell$pathway <- pathway
    slope_df <- rbind(slope_df, as.data.frame(lst_ucell))
  }
  
  

  slope_df$pathway[slope_df$pathway=="REACTOME_CELLULAR_SENESCENCE"] <- "Cellular senescence"
  slope_df$pathway[slope_df$pathway=="REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP"] <- 'SASP'
  slope_df$pathway[slope_df$pathway=="HALLMARK_APOPTOSIS"] <- "Apoptosis"
  slope_df$pathway[slope_df$pathway=="HALLMARK_DNA_REPAIR"] <- "DNA repair"
  slope_df$pathway[slope_df$pathway=="HALLMARK_INFLAMMATORY_RESPONSE"] <- "Inflammatory response"
  slope_df$pathway[slope_df$pathway=="BIOCARTA_NFKB_PATHWAY"] <- "NF-κB pathway"
  slope_df$pathway[slope_df$pathway=="BIOCARTA_INSULIN_PATHWAY"] <- "Insulin pathway"
  slope_df$pathway[slope_df$pathway=="BIOCARTA_IGF1_PATHWAY"] <- "Insulin/IGF1 pathway"
  slope_df$pathway[slope_df$pathway=="WP_AGERAGE_PATHWAY"] <- "AGE/RAGE pathway"
  slope_df$pathway[slope_df$pathway=="KEGG_MTOR_SIGNALING_PATHWAY"] <- "mTOR signaling pathway"
  slope_df$pathway[slope_df$pathway=="HALLMARK_OXIDATIVE_PHOSPHORYLATION"] <- "oxidative phosphorylation"
  slope_df$pathway[slope_df$pathway=="HALLMARK_FATTY_ACID_METABOLISM"] <- "fatty acid metabolism"
  slope_df$pathway[slope_df$pathway=="GOBP_DNA_METHYLATION"] <- "DNA methylation"
  slope_df$pathway[slope_df$pathway=="GOBP_NEGATIVE_REGULATION_OF_CELL_CYCLE"] <- "Negative regulation of cell cycle"
  slope_df$abs.trends <- abs(slope_df$age.trend)
  return(list(slope_df, signif_df))
}
a <- calculate_slope(dir1_1)

sl_df <- a[[1]]
cols <- c('OVA' = '#FD7446CC','OVI' = '#709AE1CC', 'UTE'='#FED439CC')
p <- ggplot(sl_df, aes(x = orig.ident, y = abs.trends*1000, fill = orig.ident)) + geom_bar(stat = 'identity', position = position_dodge(),width = 0.7) +
  facet_wrap(~pathway, nrow = 1, strip.position = 'right') + scale_fill_manual(values = cols) #+ coord_flip()
(p2 <- p + theme_half_open() + 
    geom_errorbar(aes(ymin = (abs.trends-SE)*1000, ymax = (abs.trends+SE)*1000), width = 0.16, linewidth = 0.3, color = 'gray') +
    theme(axis.text.x = element_text(size = 10,  color = 'black'), axis.text.y = element_text(size = 11,color = 'black'), legend.title = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_text(size = 11, color = 'black'), axis.line = element_line(size = 0.05), legend.text = element_text(size = 11), 
          legend.position = 'top', strip.background = element_rect(fill = 'white'), strip.text = element_text(size = 10))
  + labs(y = 'UCell slope'))

pdf('plots_1012/Ucell_slope.pdf',width = 12, height = 3)
p2
dev.off()

write.table(sl_df, 'Ucell_slope.txt', sep = '\t', row.names = F)
signif_df <- a[[2]]
write.table(signif_df, 'Ucell_slope_significance.txt', sep = '\t', row.names = F)

#regression between ages and average Ucell scores of three organs and compare their difference:
data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)

  colnames(data_sum)[colnames(data_sum) == 'mean'] <- varname
  return(data_sum)
}
ucell_df1 <- data.frame()
for (dir0 in dir1_1){
  dir2 <- paste(paste("../pathgroup_0727/pathgroup/", dir0, sep = ''), 'merge', sep = '/')
  pathway <- strsplit(dir0, '\\.')[[1]][2]
  for(g in groups){
    dt <- paste(dir2, paste(g, '_data.xls', sep = ''), sep = '/')
    score_df <- read.csv(dt, sep = '\t')
    score_df$pathway <- pathway
    score_df <- score_df[,c('orig.ident', 'group', 'UCell', 'pathway')]
    ucell_df1 <- rbind(ucell_df1, score_df) 
  }
  
}

ucell_df1$pathway[ucell_df1$pathway=="REACTOME_CELLULAR_SENESCENCE"] <- "Cellular senescence"
ucell_df1$pathway[ucell_df1$pathway=="REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP"] <- 'SASP'
ucell_df1$pathway[ucell_df1$pathway=="HALLMARK_APOPTOSIS"] <- "Apoptosis"
ucell_df1$pathway[ucell_df1$pathway=="HALLMARK_DNA_REPAIR"] <- "DNA repair"
ucell_df1$pathway[ucell_df1$pathway=="HALLMARK_INFLAMMATORY_RESPONSE"] <- "Inflammatory response"
ucell_df1$pathway[ucell_df1$pathway=="BIOCARTA_NFKB_PATHWAY"] <- "NF-κB pathway"
ucell_df1$pathway[ucell_df1$pathway=="BIOCARTA_INSULIN_PATHWAY"] <- "Insulin pathway"
ucell_df1$pathway[ucell_df1$pathway=="BIOCARTA_IGF1_PATHWAY"] <- "Insulin/IGF1 pathway"
ucell_df1$pathway[ucell_df1$pathway=="WP_AGERAGE_PATHWAY"] <- "AGE/RAGE pathway"
ucell_df1$pathway[ucell_df1$pathway=="KEGG_MTOR_SIGNALING_PATHWAY"] <- "mTOR signaling pathway"
ucell_df1$pathway[ucell_df1$pathway=="HALLMARK_OXIDATIVE_PHOSPHORYLATION"] <- "oxidative phosphorylation"
ucell_df1$pathway[ucell_df1$pathway=="HALLMARK_FATTY_ACID_METABOLISM"] <- "fatty acid metabolism"
ucell_df1$pathway[ucell_df1$pathway=="GOBP_DNA_METHYLATION"] <- "DNA methylation"
ucell_df1$pathway[ucell_df1$pathway=="GOBP_NEGATIVE_REGULATION_OF_CELL_CYCLE"] <- "Negative regulation of cell cycle"

stat.test <- ucell_df1 %>% 
  group_by(pathway, group) %>%
  t_test(UCell ~ orig.ident, p.adjust.method = 'fdr')

stat.test <- stat.test %>%
  add_y_position(step.increase=0.06)


ucell_df1$pathway <- factor(ucell_df1$pathway, levels = unique(ucell_df1$pathway))
ucell_df1$group  <- factor(ucell_df1$group,levels = c("3M", "6M", "9M", "12M", "15M"))
levels(ucell_df1$pathway) <- sub('_', '_\n', levels(ucell_df1$pathway))

ucell_df3 <- data_summary(ucell_df1, varname = 'UCell', groupnames = c("orig.ident","group", "pathway"))

ucell_df3$age <- 3
ucell_df3$age[ucell_df3$group == '3M'] <- 3
ucell_df3$age[ucell_df3$group == '6M'] <- 6
ucell_df3$age[ucell_df3$group == '9M'] <- 9
ucell_df3$age[ucell_df3$group == '12M'] <- 12
ucell_df3$age[ucell_df3$group == '15M'] <- 15


lm_ucell.p <- dlply(ucell_df3, .variables = .(orig.ident, pathway), .fun = function(x) {
  summary(lm(UCell ~ age, data = x))$coefficients['age', 'Pr(>|t|)'];
}) ###calculate p-values for each organ and pathway.
lm_ucell.p <- as.data.frame(unlist(lm_ucell.p))
names(lm_ucell.p) <- 'lm.p'
write.table(lm_ucell.p, 'ucell_mean_lm_p.txt', sep = '\t', quote = F)
