library(plyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(Seurat)
library(HDF5Array)
library(rhdf5)
library(hdf5r)
library(DelayedArray)
library(DelayedMatrixStats)
require(gridExtra)
library(tidyr)
library(dendextend)
library(DESeq2)
#see https://github.com/Berlin-Hamster-Single-Cell-Consortium/Single-cell-sequencing-of-COVID-19-pathogenesis-in-golden-Hamsters-
source("./smooth_DimPlot.R")
library(akima)
library(pheatmap)
#see http://www.cookbook-r.com/Graphs/ Plotting_means_and_error_bars_(ggplot2)
source("./summarySE.R")
library(cowplot)
library(SeuratObject)
library(lme4)
library(ComplexHeatmap)
library(patchwork)
library(forcats)
library(DoubletFinder)



#For more detailed thresholding, start with combined, non-integrated object with very low threshold (nGene=250)
#This was. done in step 2
exp1_nasal <- readRDS("./seu_nasal_exp1_combined_250.rds")

exp1_nasal = NormalizeData(exp1_nasal)
exp1_nasal = FindVariableFeatures(exp1_nasal, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(exp1_nasal)
exp1_nasal <- ScaleData(exp1_nasal, features = all.genes)
exp1_nasal <- RunPCA(exp1_nasal, features = VariableFeatures(object = exp1_nasal))
pdf("exp1_nasal_combined_250_elbow.pdf")
ElbowPlot(exp1_nasal)
dev.off()
exp1_nasal <- RunUMAP(exp1_nasal, dims = 1:19)
pdf("exp1_nasal_combined_250_UMAPPlot.pdf")
UMAPPlot(object=exp1_nasal)
dev.off()
exp1_nasal <- FindNeighbors(exp1_nasal, dims = 1:19)
exp1_nasal <- FindClusters(exp1_nasal, resolution = 0.9)
pdf("exp1_nasal_combined_250_clusters.pdf")
DimPlot(exp1_nasal, reduction = "umap", label=TRUE)
dev.off()

exp1_nasal@meta.data = cbind(exp1_nasal@meta.data, exp1_nasal@reductions$umap@cell.embeddings)

exp1_nasal@meta.data$treatment <- gsub("_[1-4].*","",exp1_nasal@meta.data$orig.ident)
exp1_nasal@meta.data$hamster <- gsub(".*_([1-4].*)","Ha\\1",exp1_nasal@meta.data$orig.ident)


ggplot()+geom_violin(data=exp1_nasal@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)))
ggsave("exp1_nasal_violin_nUMI_log2_in_clusters.pdf", width=12, height=7)
ggplot()+geom_violin(data=exp1_nasal@meta.data, aes(x=seurat_clusters, y=log2(nFeature_RNA)))
ggsave("exp1_nasal_violin_nGene_log2_in_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_nasal@meta.data, aes(x=seurat_clusters, y=log2(nFeature_RNA)), stat = "summary", fun = "mean")
ggsave("exp1_nasal_barplot_log2_nGene_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_nasal@meta.data, aes(x=seurat_clusters, y=log2(nFeature_RNA)), stat = "summary", fun = "median")
ggsave("exp1_nasal_barplot_log2_nGene_median_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_nasal@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)), stat = "summary", fun = "mean")
ggsave("exp1_nasal_barplot_log2_nUMI_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_nasal@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)), stat = "summary", fun = "median")
ggsave("exp1_nasal_barplot_log2_nUMI_median_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_nasal@meta.data, aes(x=seurat_clusters, y=nFeature_RNA), stat = "summary", fun = "mean")
ggsave("exp1_nasal_barplot_nGene_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_nasal@meta.data, aes(x=seurat_clusters, y=nFeature_RNA), stat = "summary", fun = "median")
ggsave("exp1_nasal_barplot_nGene_median_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_nasal@meta.data, aes(x=seurat_clusters, y=nCount_RNA), stat = "summary", fun = "mean")
ggsave("exp1_nasal_barplot_nUMI_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_nasal@meta.data, aes(x=seurat_clusters, y=nCount_RNA), stat = "summary", fun = "median")
ggsave("exp1_nasal_barplot_nUMI_median_clusters.pdf", width=12, height=7)




###############
#cell type annotation
#LOC101833790 = Ccr5

#olfactory sensory neurons (OSN): Gng13, Ppa1, 
#Immature neurons (INP): Gng8, Sox11, Elavl4
#Horizontal basal cells (HBC): Krt5, Krt7
#Globose basal cells (GBC): Ezh2, Hes6
#Sustencular cells (SUS): Lypd2, Scgb1a1, Cyp2g1 (LOC121135503), Cyp1a2 (LOC101825565)
#Microglia: C1qa, Cd74, Aif1
#Olfactory Glia: Gde1
#Glia: Apoc1, Gpm6b, Mia, Cryab, Plp1, Apoe, S100b, Fabp7
#Microvillous type 1 (MV1): Slc12a2, Hepacam2
#Microvillous type 2 (MV1): Avil, Hepacam2
#Fibroblasts: Dcn, Lum
#TNKcells: Cd3e, Cd8a, Cd4, Nkg7
#Mast cells: Ltc4s
#Top2a, Mki67
#Bcells: Cd79a, Ms4a1
#Neutrophils: Camp, S100a8, Ltf, Retn
#Macrophages: Cd14, Vcan, Cd68, Ccr2, Cx3cr1, Adgre1
#SmoothMuscle: Myl9, Myh11

markers <- c("Gng13", "Ppa1", "Gng8", "Sox11", "Elavl4", "Krt5", "Krt7", "Ezh2", "Hes6", "Lypd2", "Scgb1a1", "LOC121135503", "LOC101825565",
             "C1qa", "Cd74", "Aif1", "Gde1", "Apoc1", "Gpm6b", "Mia", "Cryab", "Plp1", "Apoe", "S100b", "Fabp7",
             "Slc12a2", "Hepacam2", "Avil", "Dcn", "Lum", "Cd3e", "Cd8a", "Cd4", "Nkg7", "Ltc4s", "Top2a", "Mki67",
             "Cd79a", "Ms4a1", "Camp", "S100a8", "Ltf", "Retn", "Cd14", "Vcan", "Cd68", "Ccr2", "Cx3cr1", "Adgre1", "Myl9", "Myh11")
gene <- "C3"
for (gene in c("Tff3", "Wfdc2", "Muc5ac")) {
  df <- FetchData(exp1_nasal, gene)
  gene <- str_replace_all(gene, "-", "_")
  colnames(df) <- gene
  df = cbind(df, exp1_nasal@reductions$umap@cell.embeddings)
  ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
    geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
    scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    coord_fixed(ratio=1)+
    ggtitle(gene)
  ggsave(paste("exp1_nasal_combined_250", gene, "pdf", sep="."), useDingbats=FALSE)
}

avg <- AverageExpression(exp1_nasal, assays=c("RNA"), features = markers, return.seurat = T, slot="data") 
avg2 <- AverageExpression(exp1_nasal, assays=c("RNA"), features = markers, return.seurat = F, slot="data")
avg2 <- as.matrix(avg2$RNA)
hc <- hclust(dist(avg2))
gene.order <- row.names(avg2)[order.hclust(hc)]
DoHeatmap(avg, size=5, features=gene.order)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("exp1_nasal_cluster_heatmap.pdf")


exp1_nasal.cluster.markers <- FindAllMarkers(exp1_nasal, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- exp1_nasal.cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>% dplyr::pull(gene)
pdf("exp1_nasal_celltypemarkers.pdf", height=6, width=36)
DotPlot(exp1_nasal, features = rev(unique(top10)))+ theme(axis.text.x = element_text(angle = 90, hjust=1))
dev.off()

#Include doublet information by reading int the DoubSing file and attach to meta data
exp1_nasal@meta.data$cell_id <- rownames(exp1_nasal@meta.data)
DoubSing <- read.table("./DoubSing_nasal_exp1.txt", header=TRUE)
exp1_nasal@meta.data <- merge(exp1_nasal@meta.data, DoubSing, by="cell_id")
rownames(exp1_nasal@meta.data) <- exp1_nasal@meta.data$cell_id


Idents(exp1_nasal) <- exp1_nasal@meta.data$seurat_clusters

#This is the annotation for the object with more than 250 genes
#pDC are likely mixed into the myDC cluster (24)
exp1_nasal <- RenameIdents(exp1_nasal, 
                          '0'='Neutrophils',
                          '1'='Macrophages',
                          '2'='mixed1',
                          '3'='Glia',
                          '4'='mixed2',
                          '5'='Bcells',
                          '6'='Neutrophils',
                          '7'='mixed3',
                          '8'='HorBasC',
                          '9'='Tcells',
                          '10'='Macrophages',
                          '11'='Macrophages',
                          '12'='mixed4',
                          '13'='OlfSensNeu',
                          '14'='GloBasC',
                          '15'='Mast',
                          '16'='ImmatureNeutrophils',
                          '17'='Fibroblasts',
                          '18'='ImmSensNeu',
                          '19'='NK',
                          '20'='Microglia',
                          '21'='Fibroblasts',
                          '22'='Glia',
                          '23'='mixed5',
                          '24'='SustC',
                          '25'='mixed6',
                          '26'='SmoothMuscle')
exp1_nasal@meta.data$celltype <- Idents(exp1_nasal)


#save at this point
saveRDS(exp1_nasal, "./exp1_nasal_combined_250_ann.rds")



#Order cell types
the_celltypes = c("Macrophages1",
                  "Macrophages2",
                  "Neutrophils",
                  "ImmatureNeutrophils",
                  "NK",
                  "Tcells",
                  "Bcells",
                  "Mast",
                  "GloBasC",
                  "HorBasC" ,
                  "SustC",
                  "OlfSensNeu",
                  "ImmSensNeu",
                  "Glia",
                  "Microglia",
                  "Fibroblasts",
                  "Microvili",
                  "SmoothMuscle",
                  "mixed1",
                  "mixed1",
                  "mixed2",
                  "mixed3",
                  "mixed4",
                  "mixed5",
                  "mixed6"
)

celltypecolors = c(
  "Macrophages" = "#B97C9D",
  "Macrophages1" = "#B97C9D",
  "Macrophages2" = "#CA8DAE",
  "Neutrophils" = "#0081AF",
  "ImmatureNeutrophils" = "#00ABE7",
  "NK" = "#4F6D7A",
  "Treg" = "#7C6A0A",
  "Tcells" = "#368F8B",
  "Bcells" = "#62C370",
  "GloBasC" = "#F7C548",
  "HorBasC" = "#F97E44",
  "SustC" = "#FF3322",
  "Mast" = "#B2675E",
  "olfBasal" = "#FB3640",
  "Fibroblasts" = "#0D3B66",
  "Microvili" = "#FFBA73",
  "OlfSensNeu" = "#C4A381",
  "ImmSensNeu" = "#644536",
  "Glia" = "#A424FF",
  "Microglia" = "#DAA9FC",
  "SmoothMuscle" = "#0D3BFF",
  "mixed1" = "#CAD2C5",
  "mixed1" = "#CAD2C5",
  "mixed2" = "#CAD2C5",
  "mixed3" = "#CAD2C5",
  "mixed4" = "#CAD2C5",
  "mixed5" = "#CAD2C5",
  "mixed6" = "#CAD2C5"
)




legendcolors <- c("aaUntr" ="gray80", "adeno"="gray65", "attenuated"="gray50", "mRNA"="gray35")

expr <- list()
for (ct in names(celltypecolors)) {
  cbright = celltypecolors[[ct]]
  r=(col2rgb(cbright)-40)[[1]]
  r=ifelse(r<0,0,r)
  g=(col2rgb(cbright)-40)[[2]]
  g=ifelse(g<0,0,g)
  b=(col2rgb(cbright)-40)[[3]]
  b=ifelse(b<0,0,b)
  cdark = rgb(r, g, b, maxColorValue = 255)
  the_scale <- scales::seq_gradient_pal(cbright, cdark, "Lab")(seq(0,1,length.out=4))
  expr[[paste(ct, "aaUntr", sep="_")]] <- the_scale[[1]]
  expr[[paste(ct, "adeno", sep="_")]] <- the_scale[[2]]
  expr[[paste(ct, "attenuated", sep="_")]] <- the_scale[[3]]
  expr[[paste(ct, "mRNA",  sep="_")]] <- the_scale[[4]]
}
the_colors = unlist(expr)


#UMAP with cell types
means <- exp1_nasal@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(mean_U1 = mean(UMAP_1), mean_U2 = mean(UMAP_2))

ggplot()+
  geom_point(data=exp1_nasal@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=celltype), size=0.5, shape=16)+
  geom_text(data=means, aes(x=mean_U1, y=mean_U2, label=celltype), size=3)+
  coord_fixed(ratio=1)+
  ggtitle("celltype")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  scale_colour_manual(values = celltypecolors)
ggsave("exp1_nasal_combined_250_celltypes.pdf", useDingbats=FALSE, width=14, height=7)




#doublets by celltype and treatment
df1 <- cbind.data.frame(exp1_nasal@meta.data$DoubSing,
                        exp1_nasal@meta.data$celltype,
                        exp1_nasal@meta.data$treatment,
                        exp1_nasal@meta.data$hamster)
colnames(df1) <- c("DoubSing", "celltype", "treatment", "hamster")

a = df1 %>% group_by(celltype, treatment, hamster) %>% tally(name="tot") 
b = df1 %>% filter(DoubSing == "Doublet") %>% group_by(celltype, treatment, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('celltype', 'treatment', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
#write.table(c, paste("days", gene, "tsv", sep="."), sep="\t", quote=FALSE, row.names = FALSE)
tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("celltype", "treatment"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")))+
  geom_bar(aes(colour=treatment), position=position_dodge(.8), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("DoubSing per treatment", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_colors)+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))
ggsave("exp1_nasal_combined_250_DoubSing.pdf", useDingbats=FALSE, width=14, height=7)

a = df1 %>% group_by(celltype, treatment, hamster) %>% tally(name="tot") 
b = df1 %>% filter(endocytotox == "doublet") %>% group_by(celltype, treatment, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('celltype', 'treatment', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
#write.table(c, paste("days", gene, "tsv", sep="."), sep="\t", quote=FALSE, row.names = FALSE)
tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("celltype", "treatment"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")))+
  geom_bar(aes(colour=treatment), position=position_dodge(.8), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("endo/cytotox per treatment", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_colors)+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))
ggsave("exp1_nasal_combined_250_endocytotox.pdf", useDingbats=FALSE, width=14, height=7)


ggplot()+geom_violin(data=exp1_nasal@meta.data %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes)), aes(x=celltype, y=log2(nCount_RNA), fill=celltype))+scale_fill_manual(values = celltypecolors)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10))
ggsave("exp1_nasal_combined_250_violin_nUMI_log2_in_celltype.pdf", width=12, height=7)
ggplot()+geom_violin(data=exp1_nasal@meta.data %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes)), aes(x=celltype, y=log2(nFeature_RNA), fill=celltype))+scale_fill_manual(values = celltypecolors)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10))
ggsave("exp1_nasal_combined_250_violin_nGene_log2_in_celltype.pdf", width=12, height=7)





#############
#final filtering

#Most cell types get normal filtering
#mixed1 and 4 are low quality, omit

expr <- list()
for (ct in c("Neutrophils", 
              "Macrophages", "Glia", "mixed2", "Bcells", "mixed3", 
              "HorBasC", "Tcells", "OlfSensNeu", "GloBasC", "Mast", 
              "ImmatureNeutrophils", "Fibroblasts", "ImmSensNeu", "NK", "Microglia", 
              "mixed5", "SustC", "mixed6", "SmoothMuscle")) {
  df <- subset(exp1_nasal@meta.data, celltype==ct)
  expr[[ct]] <- as.data.frame(subset(df, nCount_RNA>quantile(df$nCount_RNA)[[2]] & nCount_RNA < quantile(df$nCount_RNA)[[4]]+3*IQR(df$nCount_RNA))$cell_id)
  colnames(expr[[ct]]) <- "cell_id"
}  

for (ct in c("mixed1", "mixed2", "mixed3")) {
  df <- subset(exp1_nasal@meta.data, celltype==ct)
  expr[[ct]] <- as.data.frame(subset(df, nCount_RNA>quantile(subset(exp1_nasal@meta.data, celltype=="Neutrophils")$nCount_RNA)[[2]])$cell_id)
  colnames(expr[[ct]]) <- "cell_id"
}  




cells_to_keep <- do.call(rbind,expr)
write.table(cells_to_keep, "./exp1_nasal_cells_to_keep.txt", sep="\t", quote=FALSE)

write.table(exp1_nasal@meta.data, "./exp1_nasal_combined_250_metadata.txt", sep="\t", quote=FALSE)

