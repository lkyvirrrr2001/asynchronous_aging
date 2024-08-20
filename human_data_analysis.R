####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ovary~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#####1.DEG#####
deglist<-list()
Idents(human_seurat_harmony)<-'celltype_3'
#ident2对照组
for (j in celltype3) {
  deglist[[j]] <- FindMarkers(human_seurat_harmony, ident.1 = "old", ident.2 = "young",group.by = "group", subset.ident = j,logfc.threshold=0)
  print(j)
}

deglist[[1]]
deglist_filter<-list()
for (i in names(deglist)) {
  deglist_filter[[i]]<-filter(deglist[[i]],p_val_adj<0.05 & abs(avg_log2FC)>0.25)
}

table(HOM_MouseHumanSequence.rpt$Common.Organism.Name)
mouse_gene_df<-HOM_MouseHumanSequence.rpt[which(HOM_MouseHumanSequence.rpt$Common.Organism.Name %in% "mouse, laboratory"),]
human_gene_df<-HOM_MouseHumanSequence.rpt[which(HOM_MouseHumanSequence.rpt$Common.Organism.Name %in% "human"),]
colnames(human_2_mouse)
human_2_mouse<-merge(mouse_gene_df,human_gene_df,by="DB.Class.Key")
human_2_mouse<-human_2_mouse[,c("DB.Class.Key","Symbol.x","Symbol.y")]

deglist_filter_trans<-deglist_filter
for (i in names(deglist)) {
  deglist_filter_trans[[i]]$mouse<-ifelse(rownames(deglist_filter_trans[[i]]) %in% human_2_mouse$Symbol.y, 
                                          human_2_mouse$Symbol.x[match(rownames(deglist_filter_trans[[i]]),human_2_mouse$Symbol.y)],NA)
}
####COMMON DEGS WITH THE SAME DIRECTION####


for (i in names(deglist_filter_trans)) {
  df_human<-deglist_filter_trans[[i]]
  df_mouse<-mouse_degs_list_filter[[i]]
  intersect<-intersect(df_mouse$X,df_human$mouse)
  
  df_mouse_common<-df_mouse[which(df_mouse$X %in% intersect),]
  df_human_common<-df_human[which(df_human$mouse %in% intersect),]
  df_human_common<-rownames_to_column(df_human_common,"Symbol")
  colnames(df_mouse_common)<-c("symbol","mouse_p_val", "mosue_avg_logFC", "mouse_pct.1","mouse_pct.2","mouse_p_val_adj")
  colnames(df_human_common)<-c("Symbol","human_p_val", "human_avg_logFC", "human_pct.1","human_pct.2","human_p_val_adj","symbol")
  common_gene<-merge(df_mouse_common,df_human_common,by="symbol")
  common_gene$coincident<-ifelse((common_gene$mosue_avg_logFC<0 &common_gene$human_avg_logFC<0)|
                                   (common_gene$mosue_avg_logFC>0 &common_gene$human_avg_logFC>0),T,F)
  
  common_gene$cell_type<-c(rep(i,nrow(common_gene)))
  common_gene_8samples[[i]]<-common_gene
  
  common_gene_human_samedire<-df_human_common[which(df_human_common$Symbol %in% common_gene$Symbol[which(common_gene$coincident==T)]),]
  
  common_gene_human_samedire$direction<-ifelse(common_gene_human_samedire$human_avg_logFC>0,"Up","Down")
  common_gene_human_samedire$active.ident<-c(rep(i,nrow(common_gene_human_samedire)))
  common_gene_human_samedire<-common_gene_human_samedire[,c(2:6,9,1,8,7)]
  colnames(common_gene_human_samedire)<-c("p_val",	"avg_log2FC","pct.1","pct.2","p_val_adj","active.ident","gene","direction","gene_ms")
  common_gene_8samples_samedire[[i]]<-common_gene_human_samedire
  
}

Output<-data.frame()
for (i in names(common_gene_8samples_samedire)[2:6]) {
  Output<-rbind(Output,common_gene_8samples_samedire[[i]])
}
Output<-rbind(Output,common_gene_8samples_samedire[["T&S"]])

DEGs_common<-Output$gene_ms
DEGs_common_HUMAN<-Output$gene
DefaultAssay(Seurat_mouse)<-'RNA'
mouse_exp<-FetchData(Seurat_mouse,DEGs_common)
human_exp<-FetchData(human_seurat_harmony,DEGs_common_HUMAN)
human_exp_averahe<-human_exp
#存在两个人基因对应一个小鼠基因的情况，取均值#
which(colnames(human_exp_averahe) %in% "MYL12A"|colnames(human_exp_averahe) %in% "MYL12B")
which(colnames(human_exp_averahe) %in% "HSPA1A"|colnames(human_exp_averahe) %in% "HSPA1B")

mean_hspa1<-apply(human_exp_averahe[,c(12,117)],1,mean)
mean_myl12b<-apply(human_exp_averahe[,c(18,22)],1,mean)
human_exp_averahe<-human_exp_averahe[,-c(12,18,22,117)]
colnames(human_exp_averahe)<-Output$gene_ms[match(colnames(human_exp_averahe),Output$gene)]
human_exp_averahe$Hspa1a<-mean_hspa1
human_exp_averahe$Myl12b<-mean_myl12b
####整理annotate信息
human_exp_t<-as.data.frame(t(human_exp_averahe))
mouse_exp_t<-as.data.frame(t(mouse_exp))
metadata_mouse<-data.frame("Sample"=Seurat_mouse$sample[match(colnames(mouse_exp_t),rownames(Seurat_mouse@meta.data))],
                           "Cell_type"=Seurat_mouse$cluster_id[match(colnames(mouse_exp_t),rownames(Seurat_mouse@meta.data))],
                           "Species"=rep("mouse",length(colnames(mouse_exp_t))))

metadata_human<-data.frame("Sample"=human_seurat_harmony$group[match(colnames(human_exp_t),rownames(human_seurat_harmony@meta.data))],
                           "Cell_type"=human_seurat_harmony$celltype_3[match(colnames(human_exp_t),rownames(human_seurat_harmony@meta.data))],
                           "Species"=rep("human",length(colnames(human_exp_t))))

##去掉小鼠的卵母细胞
oo<-rownames(metadata_mouse[which(metadata_mouse$Cell_type %in% "OO"),])
index_oo<-which(colnames(mouse_exp_t) %in% oo)
mouse_exp_ooremove<-mouse_exp_t[,-index_oo]
metadata_mouse_ooremove<-metadata_mouse[!metadata_mouse$Cell_type %in% "OO",]
metadata_mouse_ooremove$Sample<-ifelse(metadata_mouse_ooremove$Sample %in% "OVA3M","young","old")

##merge
metadata_mushuman<-rbind(metadata_mouse_ooremove,metadata_human)
expdata_mushuman<-cbind(mouse_exp_ooremove,human_exp_t)


metadata_mushuman$Cell_type<-factor(metadata_mushuman$Cell_type,levels=c("GC","T&S","EC","Epi","IC","Mur"))
metadata_mushuman$Sample<-factor(metadata_mushuman$Sample,levels=c("young","old"))
metadata_mushuman_ord<-metadata_mushuman[order(metadata_mushuman$Cell_type,metadata_mushuman$Sample),]
expdata_mushuman2 = as.data.frame(cbind(ScaleData(mouse_exp_ooremove,features = rownames(expdata_mushuman)),
                                        ScaleData(human_exp_t,features = rownames(expdata_mushuman))))

annot_mouse_new<-data.frame(row.names = colnames(mouse_aver_new_removeoo),
                            "celltype"=rep(c("EC","Epi","GC","IC","Mur","T&S"),2),
                            "age"=c(rep("old",6),
                                    rep("young",6)),
                            "species"=c(rep("mouse",12))
)

annot_human_new<-data.frame(row.names = colnames(human_exp_averagesamegene_new),
                            "celltype"=rep(c("EC","Epi","GC","IC","Mur","T&S"),2),
                            "age"=c(rep("old",6),
                                    rep("young",6)),
                            "species"=c(rep("human",12))
)

human_exp_averagesamegene_new<-human_exp_averagesamegene_new[rownames(mouse_aver_new_removeoo),]
annot_new<-rbind(annot_mouse_new,annot_human_new)
averexp_humanandmouse_new<-as.data.frame(cbind(mouse_aver_new_removeoo,human_exp_averagesamegene_new))
annot_new$celltype<-factor(annot_new$celltype,levels=c("GC","T&S","EC","Epi","IC","Mur"))
annot_new$age<-factor(annot_new$age,levels=c("young","old"))
annot_new$species<-factor(annot_new$species,levels = c("mouse","human"))
##调整细胞顺序####
annot_new<-annot_new[order(annot_new$celltype,annot_new$age),]
annot_new<-arrange(annot_new,desc(species))
annot_new<-annot_new[,c(2,1,3)]
averexp_humanandmouse_new<-averexp_humanandmouse_new[,c(rownames(annot_new))]
##调整基因的顺序####
Output_removesample$active.ident<-factor(Output_removesample$active.ident,levels=c("GC","T&S","EC","Epi","IC","Mur"))
Output_removesample<-Output_removesample[order(Output_removesample$active.ident,Output_removesample$direction),]
gene_order_new<-Output_removesample$gene_ms
gene_order_new<-unique(gene_order_new)
averexp_humanandmouse_new<-averexp_humanandmouse_new[gene_order_new,]

####complexheatmap画图####
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(scales)
show_col(pal_nejm(palette = "default")(5))
annotation_col = HeatmapAnnotation(df = annot[,c(3,2,1)],
                                   col=list(celltype =c("GC"="#E64B35FF",
                                                        "T&S"="#4DBBD5FF", 
                                                        "EC"="#00A087FF" ,
                                                        "Epi"="#3C5488FF" ,
                                                        "IC"="#F39B7FFF" ,
                                                        "Mur"="#8491B4FF") ,
                                            species = c("mouse"="#20854EFF",
                                                        "human"="#DF8F44FF"),
                                            age = c("young"="#BC3C29FF",
                                                    "old"="#00468BFF")))

##注释bar图


df1_new<-data.frame()
for (i in rownames(averexp_humanandmouse_new)) {
  df1_new[i,1]<-sum(Output_removesample$gene_ms == i)
}
mark_new<-which(df1_new$V1==6)
markGENe_new<-rownames(df1_new)[mark_new]
df_new<-data.frame()
table(Output_removesample$active.ident)
for (i in rownames(averexp_humanandmouse_new)) {
  subset<-subset(Output_removesample,gene_ms == i)
  df_new[i,1]<-ifelse('GC' %in% subset$active.ident,1,0)
  df_new[i,2]<-ifelse('T&S' %in% subset$active.ident,1,0)
  df_new[i,3]<-ifelse('EC' %in% subset$active.ident,1,0)
  df_new[i,4]<-ifelse('Epi' %in% subset$active.ident,1,0)
  df_new[i,5]<-ifelse('IC' %in% subset$active.ident,1,0)
  df_new[i,6]<-ifelse('Mur' %in% subset$active.ident,1,0)
  
}
colnames(df_new)<-c("GC","T&S","EC","Epi","IC","Mur")

match(markGENe_new,rownames(averexp_humanandmouse_new))
####
Heatmap(averexp_humanandmouse_new,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = colorRamp2(c(seq(-2, 2, length.out = 9), 4), rev(brewer.pal(9, "RdBu"))[c(1:9,9)]),
        top_annotation = annotation_col,
        #rect_gp=gpar(col = "white", lty = 0.6, lwd = 0.6),
        column_split =rep(c("young","old"),each=12),
        show_column_names = F,
        show_row_names = F,
        row_names_gp = gpar(fontsize=6),
        #row_names_rot = -45,
        #row_names_side = "right",
        right_annotation = rowAnnotation(frequency = anno_barplot(cbind(df_new),#注释bar图
                                                                  bar_width = 0.5,#bar的宽度相对于单元格宽度的比值
                                                                  width = unit(2,"cm"),#整体bar图的宽
                                                                  axis_param=list(at=c(0,1,2,3,4,5,6),
                                                                                  labels=c(0,1,2,3,4,5,6)),#BAR的坐标label
                                                                  gp=gpar(fill=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF"),
                                                                          col=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF"))
        ), 
        foo=anno_mark(at=c(6, 12,21,26,31),#展示重点标签
                      labels = c("Serf2","Dynll1","Sec61b","Uqcrq","Uqcr10"))),
        row_title = NULL,
        name = "Expression",#颜色代表什么
        column_title = NULL)

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~oviduct~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
setwd("F:\\ThesisofAging/Nature aging/ovi/")
load("new_ovi.RData")

rm(list = ls())
library(SeuratObject)
library(Seurat)
library(tibble)
library(ggsci)
library(scales)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

Seurat_newhuman$age<-ifelse(Seurat_newhuman$reportname %in% "1714-YS-1"|Seurat_newhuman$reportname %in% "1714-YS-2"|Seurat_newhuman$reportname %in% "1714-YS-3",46,
                            ifelse(Seurat_newhuman$reportname %in% "1427-NU-1"|Seurat_newhuman$reportname %in% "1427-NU-2"|Seurat_newhuman$reportname %in% "1427-NU-3",30,
                                   ifelse(Seurat_newhuman$reportname %in% "747-YS-1"|Seurat_newhuman$reportname %in% "747-YS-2"|Seurat_newhuman$reportname %in% "747-YS-3",52,31)))


table(Seurat_newhuman$sample,Seurat_newhuman$age)
table(Seurat_newhuman$sample,Seurat_newhuman$reportname)


Seurat_newhuman<-Seurat_newhuman[,!Seurat_newhuman$age %in% "46"]
Seurat_newhuman$group<-ifelse(Seurat_newhuman$age  %in% "30"|Seurat_newhuman$age  %in% "31","young","old")
table(Seurat_newhuman$group)
human_seurat_harmony<-Seurat_newhuman
head(human_seurat_harmony@meta.data)
table(human_seurat_harmony$sample,human_seurat_harmony$age)
table(human_seurat_harmony$sample)
Idents(human_seurat_harmony)<-'cluster_standard'
####DEGS####

deglist<-list()
#ident2对照组
for (j in unique(human_seurat_harmony$cluster_standard)) {
  deglist[[j]] <- FindMarkers(human_seurat_harmony, 
                              ident.1 = "old", ident.2 = "young",
                              group.by = "group", subset.ident = j,logfc.threshold=0)
  print(j)
}

deglist_filter<-list()
for (i in names(deglist)) {
  deglist_filter[[i]]<-filter(deglist[[i]],p_val_adj<0.05 & abs(avg_log2FC)>0.25)
}


table(HOM_MouseHumanSequence.rpt$Common.Organism.Name)
mouse_gene_df<-HOM_MouseHumanSequence.rpt[which(HOM_MouseHumanSequence.rpt$Common.Organism.Name %in% "mouse, laboratory"),]
human_gene_df<-HOM_MouseHumanSequence.rpt[which(HOM_MouseHumanSequence.rpt$Common.Organism.Name %in% "human"),]
#colnames(human_2_mouse)
human_2_mouse<-merge(mouse_gene_df,human_gene_df,by="DB.Class.Key")
human_2_mouse<-human_2_mouse[,c("DB.Class.Key","Symbol.x","Symbol.y")]






deglist_filter_trans<-deglist_filter
for (i in names(deglist_filter)) {
  deglist_filter_trans[[i]]$mouse<-ifelse(rownames(deglist_filter_trans[[i]]) %in% human_2_mouse$Symbol.y, 
                                          human_2_mouse$Symbol.x[match(rownames(deglist_filter_trans[[i]]),human_2_mouse$Symbol.y)],NA)
}
#save.image(file="new_ovi.RData")
mouse_degs_list<-list.files(path = "F:\\ThesisofAging\\Nature aging\\ovi_deg",pattern = "marker_df_OVI15MvsOVI3M_")
mouse_degs_list1<-lapply(mouse_degs_list, function(x){a=read_tsv(file.path('F:\\ThesisofAging\\Nature aging\\ovi_deg',x))})
table(Seurat_newhuman$cluster_standard)
names(mouse_degs_list1)<-mouse_degs_list
mouse_degs_list1<-mouse_degs_list1[c(1,2,3,5,6,7)]
names(mouse_degs_list1)<-c("EC","Epi","IC","PC","SC","SMC")

mouse_degs_list_filter<-list()
for (i in names(mouse_degs_list1)) {
  mouse_degs_list_filter[[i]]<-filter(mouse_degs_list1[[i]],p_val_adj<0.05 & abs(avg_logFC)>0.25)
  colnames(mouse_degs_list_filter[[i]])[1]<-"X"
}



####COMMON DEGS####
commongenes<-list()
commongenes_samedire<-list()
for (i in names(deglist_filter_trans)) {
  df_human<-deglist_filter_trans[[i]]
  df_mouse<-mouse_degs_list_filter[[i]]
  intersect<-intersect(df_mouse$X,df_human$mouse)
  
  df_mouse_common<-df_mouse[which(df_mouse$X %in% intersect),]
  df_human_common<-df_human[which(df_human$mouse %in% intersect),]
  df_human_common<-rownames_to_column(df_human_common,"Symbol")
  colnames(df_mouse_common)<-c("symbol","mouse_p_val", "mosue_avg_logFC", "mouse_pct.1","mouse_pct.2","mouse_p_val_adj")
  colnames(df_human_common)<-c("Symbol","human_p_val", "human_avg_logFC", "human_pct.1","human_pct.2","human_p_val_adj","symbol")
  common_gene<-merge(df_mouse_common,df_human_common,by="symbol")
  common_gene$coincident<-ifelse((common_gene$mosue_avg_logFC<0 &common_gene$human_avg_logFC<0)|
                                   (common_gene$mosue_avg_logFC>0 &common_gene$human_avg_logFC>0),T,F)
  
  common_gene$cell_type<-c(rep(i,nrow(common_gene)))
  commongenes[[i]]<-common_gene
  
  common_gene_human_samedire<-df_human_common[which(df_human_common$Symbol %in% common_gene$Symbol[which(common_gene$coincident==T)]),]
  
  common_gene_human_samedire$direction<-ifelse(common_gene_human_samedire$human_avg_logFC>0,"Up","Down")
  common_gene_human_samedire$active.ident<-c(rep(i,nrow(common_gene_human_samedire)))
  common_gene_human_samedire<-common_gene_human_samedire[,c(2:6,9,1,8,7)]
  colnames(common_gene_human_samedire)<-c("p_val",	"avg_log2FC","pct.1","pct.2","p_val_adj","active.ident","gene","direction","gene_ms")
  commongenes_samedire[[i]]<-common_gene_human_samedire
  
}

Output<-data.frame()
for (i in names(commongenes_samedire)) {
  Output<-rbind(Output,commongenes_samedire[[i]])
}
write.csv(Output,"human_mouse_gene_ovi_0621.csv")

unique(Output$gene_ms)#374
DEGs_common<-Output$gene_ms
DEGs_common_human<-Output$gene



head(Seurat_mouse@meta.data)
head(human_seurat_harmony@meta.data)
table(Seurat_mouse$active.ident)
table(human_seurat_harmony$celltype)
table(Output$active.ident)



Seurat_mouse$age_celltype<-paste0(Seurat_mouse$sample,"_",Seurat_mouse$cluster)
table(Seurat_mouse$age_celltype)


mouse_aver<- as.data.frame(AverageExpression(Seurat_mouse,
                                             features = DEGs_common,
                                             group.by = "age_celltype",
                                             assays = "RNA",
                                             slot = 'data'))

mouse_aver<-as.data.frame(ScaleData(mouse_aver,features = rownames(mouse_aver) ))



human_seurat_harmony$age_celltype<-paste0(human_seurat_harmony$group,"_",human_seurat_harmony$cluster_standard) 
table(human_seurat_harmony$age_celltype)
human_aver<- as.data.frame(AverageExpression(human_seurat_harmony,
                                             features = DEGs_common_human,
                                             group.by = "age_celltype",
                                             assays = "RNA",
                                             slot = 'data')) 




human_trans<-human_aver
rownames(human_trans)<-Output$gene_ms[match(rownames(human_trans),Output$gene)]
Output$gene[which(Output$gene_ms %in% c("Ccl3","Fcgr2b","Hspa1a","Mt1","Pirb","Trbc1"))]
Output$gene[which(Output$gene_ms %in% c(#"Ccl3","Fcgr2b""Hspa1a"#,"Mt1"#"Pirb",
  "Trbc1"
))]


#Ccl3
which(rownames(human_aver) %in% "CCL3"|rownames(human_aver) %in% "CCL3L1")
mean_CCL3<-t(as.data.frame(apply(human_aver[c(147,209),],2,mean)))
rownames(mean_CCL3)<-"Ccl3"

#Fcgr2b
which(rownames(human_aver) %in% "FCGR2A"|rownames(human_aver) %in% "FCGR2B")
mean_Fcgr2b<-t(as.data.frame(apply(human_aver[c(180,308),],2,mean)))
rownames(mean_Fcgr2b)<-"Fcgr2b"
#Hspa1a
which(rownames(human_aver) %in% "HSPA1A"|rownames(human_aver) %in% "HSPA1B")
mean_Hspa1a<-t(as.data.frame(apply(human_aver[c(124,125),],2,mean)))
rownames(mean_Hspa1a)<-"Hspa1a"
#mt1
which(rownames(human_aver) %in% "MT1E"|rownames(human_aver) %in% "MT1M"|rownames(human_aver) %in% "MT1F"|
        rownames(human_aver) %in% "MT1X"|rownames(human_aver) %in% "MT1A"|rownames(human_aver) %in% "MT1G")

mean_Mt1<-t(as.data.frame(apply(human_aver[c(103,105,110,353),],2,mean)))
rownames(mean_Mt1)<-"Mt1"
#Pirb

which(rownames(human_aver) %in% "LILRB2"|rownames(human_aver) %in% "LILRB1")
mean_Pirb<-t(as.data.frame(apply(human_aver[c(190,222),],2,mean)))
rownames(mean_Pirb)<-"Pirb"

#Trbc1


which(rownames(human_aver) %in% "TRBC2"|rownames(human_aver) %in% "TRBC1")
mean_Trbc1<-t(as.data.frame(apply(human_aver[c(153,310),],2,mean)))
rownames(mean_Trbc1)<-"Trbc1"


human_exp_averagesamegene<-human_aver[-c(147,209,180,308,124,125,103,105,110,353,190,222,153,310),]
rownames(human_exp_averagesamegene)<-Output$gene_ms[match(rownames(human_exp_averagesamegene),Output$gene)]

human_exp_averagesamegene<-rbind(human_exp_averagesamegene,mean_CCL3)
human_exp_averagesamegene<-rbind(human_exp_averagesamegene,mean_Fcgr2b)
human_exp_averagesamegene<-rbind(human_exp_averagesamegene,mean_Hspa1a)
human_exp_averagesamegene<-rbind(human_exp_averagesamegene,mean_Pirb)
human_exp_averagesamegene<-rbind(human_exp_averagesamegene,mean_Mt1)
human_exp_averagesamegene<-rbind(human_exp_averagesamegene,mean_Trbc1)

human_exp_averagesamegene<-as.data.frame(ScaleData(human_exp_averagesamegene,features = rownames(human_exp_averagesamegene)) )
###注释信息
annot_mouse<-data.frame(row.names = colnames(mouse_aver),
                        "celltype"=rep(c("EC","Epi","IC","PC","SC","SMC"),2),
                        "age"=c(rep("old",6),
                                rep("young",6)),
                        "species"=c(rep("mouse",12))
)

annot_human<-data.frame(row.names = colnames(human_exp_averagesamegene),
                        "celltype"=rep(c("EC","Epi","IC","PC","SC","SMC"),2),
                        "age"=c(rep("old",6),
                                rep("young",6)),
                        "species"=c(rep("human",12))
)

annot<-rbind(annot_mouse,annot_human)
human_exp_averagesamegene<-human_exp_averagesamegene[rownames(mouse_aver),]
averexp_humanandmouse<-as.data.frame(cbind(mouse_aver,human_exp_averagesamegene))
annot$celltype<-factor(annot$celltype,levels=c("EC","Epi","IC","PC","SC","SMC"))
annot$age<-factor(annot$age,levels=c("young","old"))
annot$species<-factor(annot$species,levels = c("mouse","human"))


annot<-annot[order(annot$celltype,annot$age),]
annot<-arrange(annot,desc(species))
annot<-annot[,c(2,1,3)]
averexp_humanandmouse<-averexp_humanandmouse[,c(rownames(annot))]
##调整基因的顺序####
Output$active.ident<-factor(Output$active.ident,levels=c("EC","Epi","IC","PC","SC","SMC"))
Output<-Output[order(Output$active.ident,Output$direction),]
gene_order<-Output$gene_ms
gene_order<-unique(gene_order)
averexp_humanandmouse<-averexp_humanandmouse[gene_order,]
table(human_seurat_harmony$age)
###画图####
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(scales)
show_col(pal_igv(palette = "default")(10))

annotation_col = HeatmapAnnotation(df = annot[,c(3,2,1)],
                                   col=list(celltype =c("EC"="#802268ff",
                                                        "Epi"="#ce3d32ff", 
                                                        "IC"="#749b58ff" ,
                                                        "PC"="#f0e685ff" ,
                                                        "SC"="#d595a7ff" ,
                                                        "SMC"="#5db1ddff") ,
                                            species = c("mouse"="#20854EFF",
                                                        "human"="#DF8F44FF"),
                                            age = c("young"="#BC3C29FF",
                                                    "old"="#00468BFF")))

##注释bar图

df1<-data.frame()
for (i in rownames(averexp_humanandmouse)) {
  df1[i,1]<-sum(Output$gene_ms == i)
}
mark<-which(df1$V1==6)
markGENe<-rownames(df1)[mark]
df<-data.frame()
for (i in rownames(averexp_humanandmouse)) {
  subset<-subset(Output,gene_ms == i)
  df[i,1]<-ifelse('EC' %in% subset$active.ident,1,0)
  df[i,2]<-ifelse('Epi' %in% subset$active.ident,1,0)
  df[i,3]<-ifelse('IC' %in% subset$active.ident,1,0)
  df[i,4]<-ifelse('PC' %in% subset$active.ident,1,0)
  df[i,5]<-ifelse('SC' %in% subset$active.ident,1,0)
  df[i,6]<-ifelse('SMC' %in% subset$active.ident,1,0)
  
}
colnames(df)<-c("EC","Epi","IC","PC","SC","SMC")

match(markGENe,rownames(averexp_humanandmouse))
####
Heatmap(averexp_humanandmouse,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = colorRamp2(c(seq(-2, 2, length.out = 9), 4), rev(brewer.pal(9, "RdBu"))[c(1:9,9)]),
        top_annotation = annotation_col,
        #rect_gp=gpar(col = "white", lty = 0.6, lwd = 0.6),
        column_split =rep(c("young","old"),each=12),
        show_column_names = F,
        show_row_names = F,
        row_names_gp = gpar(fontsize=6),
        #row_names_rot = -45,
        #row_names_side = "right",
        right_annotation = rowAnnotation(frequency = anno_barplot(cbind(df),#注释bar图
                                                                  ylim = c(0,6),#坐标轴范围
                                                                  bar_width = 0.5,#bar的宽度相对于单元格宽度的比值
                                                                  width = unit(2,"cm"),#整体bar图的宽
                                                                  axis_param=list(at=c(0,1,2,3,4,5,6),
                                                                                  labels=c(0,1,2,3,4,5,6)),#BAR的坐标label
                                                                  gp=gpar(fill=c("#802268ff","#ce3d32ff","#749b58ff","#f0e685ff","#d595a7ff","#5db1ddff"),
                                                                          col=c("#802268ff","#ce3d32ff","#749b58ff","#f0e685ff","#d595a7ff","#5db1ddff"))
        ), 
        foo=anno_mark(at=c(6,23),#展示重点标签
                      labels = c("Serpinh1","Rps27"))),
        row_title = NULL,
        name = "Expression",#颜色代表什么
        column_title = NULL)


#####~~~~~~~~~~~~~~~~~~~~~~~~~~~uterus~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

library(SeuratObject)
library(Seurat)
library(tibble)
library(ggsci)
library(scales)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
head(ute_hm@meta.data)
table(ute_hm$celltype)

Idents(ute_hm)<-'celltype'


table(ute_hm$age,ute_hm$celltype)
####common degs####
deglist_filter<-list()
for (i in names(deglist)) {
  deglist_filter[[i]]<-filter(deglist[[i]],p_val_adj<0.05 & abs(avg_log2FC)>0.25)
}


table(HOM_MouseHumanSequence.rpt$Common.Organism.Name)
mouse_gene_df<-HOM_MouseHumanSequence.rpt[which(HOM_MouseHumanSequence.rpt$Common.Organism.Name %in% "mouse, laboratory"),]
human_gene_df<-HOM_MouseHumanSequence.rpt[which(HOM_MouseHumanSequence.rpt$Common.Organism.Name %in% "human"),]
colnames(human_2_mouse)
human_2_mouse<-merge(mouse_gene_df,human_gene_df,by="DB.Class.Key")
human_2_mouse<-human_2_mouse[,c("DB.Class.Key","Symbol.x","Symbol.y")]

deglist_filter_trans<-deglist_filter
for (i in names(deglist_filter)) {
  deglist_filter_trans[[i]]$mouse<-rownames( deglist_filter_trans[[i]])
}


mouse_degs_list<-list.files(path = "F:\\ThesisofAging\\Nature aging\\ute_deg/",pattern = "marker_df_UTE15MvsUTE3M_")
mouse_degs_list1<-lapply(mouse_degs_list, function(x){a=read_tsv(file.path('F:\\ThesisofAging\\Nature aging\\ute_deg/',x))})
table(human_seurat_harmony$celltype)
names(mouse_degs_list1)<-mouse_degs_list
#mouse_degs_list1<-mouse_degs_list1[c(1,2,3,5,6,7)]
names(mouse_degs_list1)<-c("EC","IC","Myo","PC","SC")

mouse_degs_list_filter<-list()
for (i in names(mouse_degs_list1)) {
  mouse_degs_list_filter[[i]]<-filter(mouse_degs_list1[[i]],p_val_adj<0.05 & abs(avg_logFC)>0.25)
  colnames(mouse_degs_list_filter[[i]])[1]<-"X"
}




commongenes<-list()
commongenes_samedire<-list()
for (i in names(deglist_filter_trans)) {
  df_human<-deglist_filter_trans[[i]]
  df_mouse<-mouse_degs_list_filter[[i]]
  intersect<-intersect(df_mouse$X,df_human$mouse)
  
  df_mouse_common<-df_mouse[which(df_mouse$X %in% intersect),]
  df_human_common<-df_human[which(df_human$mouse %in% intersect),]
  df_human_common<-rownames_to_column(df_human_common,"Symbol")
  colnames(df_mouse_common)<-c("symbol","mouse_p_val", "mosue_avg_logFC", "mouse_pct.1","mouse_pct.2","mouse_p_val_adj")
  colnames(df_human_common)<-c("Symbol","human_p_val", "human_avg_logFC", "human_pct.1","human_pct.2","human_p_val_adj","symbol")
  common_gene<-merge(df_mouse_common,df_human_common,by="symbol")
  common_gene$coincident<-ifelse((common_gene$mosue_avg_logFC<0 &common_gene$human_avg_logFC<0)|
                                   (common_gene$mosue_avg_logFC>0 &common_gene$human_avg_logFC>0),T,F)
  
  common_gene$cell_type<-c(rep(i,nrow(common_gene)))
  commongenes[[i]]<-common_gene
  
  common_gene_human_samedire<-df_human_common[which(df_human_common$Symbol %in% common_gene$Symbol[which(common_gene$coincident==T)]),]
  
  common_gene_human_samedire$direction<-ifelse(common_gene_human_samedire$human_avg_logFC>0,"Up","Down")
  common_gene_human_samedire$active.ident<-c(rep(i,nrow(common_gene_human_samedire)))
  common_gene_human_samedire<-common_gene_human_samedire[,c(2:6,9,1,8,7)]
  colnames(common_gene_human_samedire)<-c("p_val",	"avg_log2FC","pct.1","pct.2","p_val_adj","active.ident","gene","direction","gene_ms")
  commongenes_samedire[[i]]<-common_gene_human_samedire
  
}

Output<-data.frame()
for (i in names(commongenes_samedire)) {
  Output<-rbind(Output,commongenes_samedire[[i]])
}

table(Output$active.ident)
write.csv(Output,file = "ute_commongene.csv")

#save.image("ute_commongene.RData")

Seurat_mouse <- readRDS("F:\\RCD data\\RDS/UTE.rds")
Seurat_mouse<-UpdateSeuratObject(Seurat_mouse)
Seurat_mouse<-Seurat_mouse[,Seurat_mouse$sample=="UTE3M"|Seurat_mouse$sample=="UTE15M"]
table(Output$active.ident)
Seurat_mouse<-Seurat_mouse[,Seurat_mouse$cluster %in% "EC"|
                             Seurat_mouse$cluster %in% "IC"|
                             Seurat_mouse$cluster %in% "Myo"|
                             Seurat_mouse$cluster %in% "PC"|
                             Seurat_mouse$cluster %in% "SC"]

human_seurat_harmony<-ute_hm[,ute_hm$celltype %in% "EC"|
                               ute_hm$celltype %in% "IC"|
                               ute_hm$celltype %in% "Myo"|
                               ute_hm$celltype %in% "PC"|
                               ute_hm$celltype %in% "SC"]

head(Seurat_mouse@meta.data)
head(human_seurat_harmony@meta.data)
table(Seurat_mouse$active.ident)
table(human_seurat_harmony$celltype)
table(Output$active.ident)



Seurat_mouse$age_celltype<-paste0(Seurat_mouse$sample,"_",Seurat_mouse$cluster)
table(Seurat_mouse$age_celltype)
DEGs_common<-Output$gene_ms
library(Seurat)
mouse_aver<- as.data.frame(AverageExpression(Seurat_mouse,
                                             features = DEGs_common,
                                             group.by = "age_celltype",
                                             assays = "RNA",
                                             slot = 'data'))

mouse_aver<-as.data.frame(ScaleData(mouse_aver,features = rownames(mouse_aver) ))



human_seurat_harmony$age_celltype<-paste0(human_seurat_harmony$age,"_",human_seurat_harmony$celltype)              
human_aver<- as.data.frame(AverageExpression(human_seurat_harmony,
                                             features = DEGs_common,
                                             group.by = "age_celltype",
                                             assays = "RNA",
                                             slot = 'data')) 





human_exp_averagesamegene<-as.data.frame(ScaleData(human_aver,features = rownames(human_aver)) )
#save.image("ute_commongene.RData")

###注释信息####
annot_mouse<-data.frame(row.names = colnames(mouse_aver),
                        "celltype"=rep(c("EC","IC","Myo","PC","SC"),2),
                        "age"=c(rep("old",5),
                                rep("young",5)),
                        "species"=c(rep("mouse",10))
)

annot_human<-data.frame(row.names = colnames(human_exp_averagesamegene),
                        "celltype"=rep(c("EC","IC","Myo","PC","SC"),2),
                        "age"=c(rep("young",5),
                                rep("old",5)),
                        "species"=c(rep("human",10))
)

annot<-rbind(annot_mouse,annot_human)
human_exp_averagesamegene<-human_exp_averagesamegene[rownames(mouse_aver),]
averexp_humanandmouse<-as.data.frame(cbind(mouse_aver,human_exp_averagesamegene))
annot$celltype<-factor(annot$celltype,levels=c("EC","IC","Myo","PC","SC"))
annot$age<-factor(annot$age,levels=c("young","old"))
annot$species<-factor(annot$species,levels = c("mouse","human"))


annot<-annot[order(annot$celltype,annot$age),]
annot<-arrange(annot,desc(species))
annot<-annot[,c(2,1,3)]
averexp_humanandmouse<-averexp_humanandmouse[,c(rownames(annot))]
##调整基因的顺序####
Output$active.ident<-factor(Output$active.ident,levels=c("EC","IC","Myo","PC","SC"))
Output<-Output[order(Output$active.ident,Output$direction),]
gene_order<-Output$gene_ms
gene_order<-unique(gene_order)
averexp_humanandmouse<-averexp_humanandmouse[gene_order,]

###画图####
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(scales)
show_col(pal_nejm(palette = "default")(10))

annotation_col = HeatmapAnnotation(df = annot[,c(3,2,1)],
                                   col=list(celltype =c("EC"="#BC3C29FF",
                                                        "IC"="#0072B5FF", 
                                                        "Myo"="#FFDC91FF" ,
                                                        "PC"="#20854EFF" ,
                                                        "SC"="#7876B1FF") ,
                                            species = c("mouse"="#20854EFF",
                                                        "human"="#DF8F44FF"),
                                            age = c("young"="#BC3C29FF",
                                                    "old"="#00468BFF")))

##注释bar图

df1<-data.frame()
for (i in rownames(averexp_humanandmouse)) {
  df1[i,1]<-sum(Output$gene_ms == i)
}

mark<-which(df1$V1==5|df1$V1==4)
markGENe<-rownames(df1)[mark]
df<-data.frame()
for (i in rownames(averexp_humanandmouse)) {
  subset<-subset(Output,gene_ms == i)
  df[i,1]<-ifelse('EC' %in% subset$active.ident,1,0)
  df[i,2]<-ifelse('IC' %in% subset$active.ident,1,0)
  df[i,3]<-ifelse('Myo' %in% subset$active.ident,1,0)
  df[i,4]<-ifelse('PC' %in% subset$active.ident,1,0)
  df[i,5]<-ifelse('SC' %in% subset$active.ident,1,0)
  
}
colnames(df)<-c("EC","IC","Myo","PC","SC")

match(markGENe,rownames(averexp_humanandmouse))
####
Heatmap(averexp_humanandmouse,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = colorRamp2(c(seq(-2, 2, length.out = 9), 4), rev(brewer.pal(9, "RdBu"))[c(1:9,9)]),
        top_annotation = annotation_col,
        #rect_gp=gpar(col = "white", lty = 0.6, lwd = 0.6),
        column_split =rep(c("young","old"),each=10),
        show_column_names = F,
        show_row_names = F,
        row_names_gp = gpar(fontsize=6),
        #row_names_rot = -45,
        #row_names_side = "right",
        right_annotation = rowAnnotation(frequency = anno_barplot(cbind(df),#注释bar图
                                                                  ylim = c(0,5),#坐标轴范围
                                                                  bar_width = 0.5,#bar的宽度相对于单元格宽度的比值
                                                                  width = unit(2,"cm"),#整体bar图的宽
                                                                  axis_param=list(at=c(0,1,2,3,4,5),
                                                                                  labels=c(0,1,2,3,4,5)),#BAR的坐标label
                                                                  gp=gpar(fill=c("#BC3C29FF","#0072B5FF","#FFDC91FF","#20854EFF","#7876B1FF"),
                                                                          col=c("#BC3C29FF","#0072B5FF","#FFDC91FF","#20854EFF","#7876B1FF"))
        ), 
        foo=anno_mark(at=c(23,27,31,40,49,61),#展示重点标签
                      labels = c("Mt1","Gem","Eif1","Cebpb","Cdkn1a","Gadd45b"))),
        row_title = NULL,
        name = "Expression",#颜色代表什么
        column_title = NULL)

save.image("ute_commongene.RData")

