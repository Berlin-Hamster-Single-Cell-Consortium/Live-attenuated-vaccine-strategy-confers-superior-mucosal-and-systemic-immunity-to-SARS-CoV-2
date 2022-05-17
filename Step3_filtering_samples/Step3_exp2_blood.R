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
#This was done in step 2
exp2_blood <- readRDS("./seu_blood_exp2_combined_250.rds")

exp2_blood = NormalizeData(exp2_blood)
exp2_blood = FindVariableFeatures(exp2_blood, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(exp2_blood)
exp2_blood <- ScaleData(exp2_blood, features = all.genes)
exp2_blood <- RunPCA(exp2_blood, features = VariableFeatures(object = exp2_blood))
pdf("exp2_blood_combined_250_elbow.pdf")
ElbowPlot(exp2_blood)
dev.off()
exp2_blood <- RunUMAP(exp2_blood, dims = 1:19)
pdf("exp2_blood_combined_250_UMAPPlot.pdf")
UMAPPlot(object=exp2_blood)
dev.off()
exp2_blood <- FindNeighbors(exp2_blood, dims = 1:19)
exp2_blood <- FindClusters(exp2_blood, resolution = 0.9)
pdf("exp2_blood_combined_250_clusters.pdf")
DimPlot(exp2_blood, reduction = "umap", label=TRUE)
dev.off()

exp2_blood@meta.data = cbind(exp2_blood@meta.data, exp2_blood@reductions$umap@cell.embeddings)

exp2_blood@meta.data$treatment <- gsub("_[1-4].*","",exp2_blood@meta.data$orig.ident)
exp2_blood@meta.data$hamster <- gsub(".*_([1-4].*)","Ha\\1",exp2_blood@meta.data$orig.ident)


ggplot()+geom_violin(data=exp2_blood@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)))
ggsave("exp2_blood_violin_nUMI_log2_in_clusters.pdf", width=12, height=7)
ggplot()+geom_violin(data=exp2_blood@meta.data, aes(x=seurat_clusters, y=log2(nFeature_RNA)))
ggsave("exp2_blood_violin_nGene_log2_in_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp2_blood@meta.data, aes(x=seurat_clusters, y=log2(nFeature_RNA)), stat = "summary", fun = "mean")
ggsave("exp2_blood_barplot_log2_nGene_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp2_blood@meta.data, aes(x=seurat_clusters, y=log2(nFeature_RNA)), stat = "summary", fun = "median")
ggsave("exp2_blood_barplot_log2_nGene_median_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp2_blood@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)), stat = "summary", fun = "mean")
ggsave("exp2_blood_barplot_log2_nUMI_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp2_blood@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)), stat = "summary", fun = "median")
ggsave("exp2_blood_barplot_log2_nUMI_median_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp2_blood@meta.data, aes(x=seurat_clusters, y=nFeature_RNA), stat = "summary", fun = "mean")
ggsave("exp2_blood_barplot_nGene_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp2_blood@meta.data, aes(x=seurat_clusters, y=nFeature_RNA), stat = "summary", fun = "median")
ggsave("exp2_blood_barplot_nGene_median_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp2_blood@meta.data, aes(x=seurat_clusters, y=nCount_RNA), stat = "summary", fun = "mean")
ggsave("exp2_blood_barplot_nUMI_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp2_blood@meta.data, aes(x=seurat_clusters, y=nCount_RNA), stat = "summary", fun = "median")
ggsave("exp2_blood_barplot_nUMI_median_clusters.pdf", width=12, height=7)




###############
#cell type annotation
#LOC101833790 = Ccr5

markers <- c("Ltf", "Adgre1", "Ccr2", "Cx3cr1", "Tcf4", "Irf8", "Flt3", "Retn", "Jchain", "Sdc1", "Cd38", "Slc3a2", "Slc7a5", "Camp", "S100a8", "Cxcr2", "Gng11", "Ppbp", "Cd79b", "Ms4a1", "Nkg7", "Gzma", "Il7r", "Cd8a", "Cd4", "Cd3e", "Cd27", "Tnfrsf17", "Prdm1", "Xbp1", "Irf4", "Sec11c", "Fkbp11", "Mki67", "Top2a")
for (gene in markers) {
  df <- FetchData(exp2_blood, gene)
  gene <- str_replace_all(gene, "-", "_")
  colnames(df) <- gene
  df = cbind(df, exp2_blood@reductions$umap@cell.embeddings)
  ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
    geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
    scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    coord_fixed(ratio=1)+
    ggtitle(gene)
  ggsave(paste("exp2_blood_combined_250", gene, "pdf", sep="."), useDingbats=FALSE)
}

avg <- AverageExpression(exp2_blood, assays=c("RNA"), features = markers, return.seurat = T, slot="data") 
avg2 <- AverageExpression(exp2_blood, assays=c("RNA"), features = markers, return.seurat = F, slot="data")
avg2 <- as.matrix(avg2$RNA)
hc <- hclust(dist(avg2))
gene.order <- row.names(avg2)[order.hclust(hc)]
DoHeatmap(avg, size=5, features=gene.order)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("exp2_blood_combined_250_cluster_heatmap.pdf", width=12, height=7)

cluster.markers <- FindAllMarkers(exp2_blood, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("exp2_blood_combined_250_cluster_markers_heatmap.pdf", width=24, height=16)
DoHeatmap(exp2_blood, features = top10$gene) + NoLegend()
dev.off()


#Include doublet information by reading int the DoubSing file and attach to meta data
exp2_blood@meta.data$cell_id <- rownames(exp2_blood@meta.data)
DoubSing <- read.table("./DoubSing_blood_exp2.txt", header=TRUE)
df <- dplyr::left_join(exp2_blood@meta.data, DoubSing, by="cell_id")
exp2_blood@meta.data$DoubSing <- df$DoubSing

Idents(exp2_blood) <- exp2_blood@meta.data$seurat_clusters

cluster.markers <- FindAllMarkers(exp2_blood, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("exp2_blood_combined_250_cluster_markers_heatmap.pdf", width=24, height=16)
DoHeatmap(exp2_blood, features = top10$gene) + NoLegend()
dev.off()


Idents(exp2_blood) <- exp2_blood@meta.data$seurat_clusters

#This is the annotation for the object with more than 250 genes
exp2_blood <- RenameIdents(exp2_blood, 
                          '0'='Neutrophil',
                          '1'='B',
                          '2'='Neutrophil',
                          '3'='Immature neutrophil',
                          '4'='Classical monocyte',
                          '5'='T',
                          '6'='Platelet',
                          '7'='B',
                          '8'='Neutrophil',
                          '9'='mixed1',
                          '10'='T',
                          '11'='mixed2',
                          '12'='B',
                          '13'='Non-classical monocyte',
                          '14'='mixed3',
                          '15'='B',
                          '16'='mixed4',
                          '17'='NK',
                          '18'='mixed5',
                          '19'='Platelet',
                          '20'='mDC',
                          '21'='Activated T',
                          '22'='mixed6',
                          '23'='Neutrophil',
                          '24'='mixed6',
                          '25'='B',
                          '26'='T',
                          '27'='pDC')
exp2_blood@meta.data$celltype <- Idents(exp2_blood)


#save at this point
saveRDS(exp2_blood, "./exp2_blood_combined_250_ann.rds")



#Order cell types
the_celltypes = c("Classical monocyte",
                  "Non-classical monocyte",
                  "Neutrophil",
                  "Immature neutrophil",
                  "mDC",
                  "pDC",
                  "NK",
                  "T",
                  "Activated T",
                  "B",
                  "Platelet")


celltypecolors<- c("T" = "#368F8B", "Activated T"= "#5CC1BC", "B" = "#62C370", "Classical monocyte" = "#B7245C",
                   "Non-classical monocyte" ="#3E2F5B", "Neutrophil" = "#0081AF",
                   "Immature neutrophil" = "#00ABE7", "NK"= "#246A73", "pDC" = "#7C6A0A",
                   "mDC" = "#4F6D7A", "Platelet" = "#832D52", "unknown" = "#CAD2C5", "mixed1" = "#CAD2C5",
                   "mixed1" = "#CAD2C5",
                   "mixed2" = "#CAD2C5",
                   "mixed3" = "#CAD2C5",
                   "mixed4" = "#CAD2C5",
                   "mixed5" = "#CAD2C5",
                   "mixed6" = "#CAD2C5",
                   "mixed7" = "#CAD2C5",
                   "mixed8" = "#CAD2C5")


legendcolors <- c("aaUntr" ="gray80", "adeno2x"="gray65", "att2x"="gray50", "mRNA2x"="gray35", "mRNAatt"="gray20")

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
  the_scale <- scales::seq_gradient_pal(cbright, cdark, "Lab")(seq(0,1,length.out=5))
  expr[[paste(ct, "aaUntr", sep="_")]] <- the_scale[[1]]
  expr[[paste(ct, "adeno2x", sep="_")]] <- the_scale[[2]]
  expr[[paste(ct, "att2x", sep="_")]] <- the_scale[[3]]
  expr[[paste(ct, "mRNA2x",  sep="_")]] <- the_scale[[4]]
  expr[[paste(ct, "mRNAatt",  sep="_")]] <- the_scale[[5]]
}
the_colors = unlist(expr)


#UMAP with cell types
means <- exp2_blood@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(mean_U1 = mean(UMAP_1), mean_U2 = mean(UMAP_2))

ggplot()+
  geom_point(data=exp2_blood@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=celltype), size=0.5, shape=16)+
  geom_text(data=means, aes(x=mean_U1, y=mean_U2, label=celltype), size=3)+
  coord_fixed(ratio=1)+
  ggtitle("celltype")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  scale_colour_manual(values = celltypecolors)
ggsave("exp2_blood_combined_250_celltypes.pdf", useDingbats=FALSE, width=14, height=7)




#doublets by celltype and treatment
df1 <- cbind.data.frame(exp2_blood@meta.data$DoubSing,
                        exp2_blood@meta.data$celltype,
                        exp2_blood@meta.data$treatment,
                        exp2_blood@meta.data$hamster)
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
ggsave("exp2_blood_combined_250_DoubSing.pdf", useDingbats=FALSE, width=14, height=7)



ggplot()+geom_violin(data=exp2_blood@meta.data %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes)), aes(x=celltype, y=log2(nCount_RNA), fill=celltype))+scale_fill_manual(values = celltypecolors)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10))
ggsave("exp2_blood_combined_250_violin_nUMI_log2_in_celltype.pdf", width=12, height=7)
ggplot()+geom_violin(data=exp2_blood@meta.data %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes)), aes(x=celltype, y=log2(nFeature_RNA), fill=celltype))+scale_fill_manual(values = celltypecolors)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10))
ggsave("exp2_blood_combined_250_violin_nGene_log2_in_celltype.pdf", width=12, height=7)




#############
#final filtering

#Most cell types get normal filtering
#mixed3 and 5 are doublet cells, i.e. no outlier removal
#mixed 1 and 2 are low quality, remove entirely

expr <- list()
for (ct in c("Classical monocyte",
             "Non-classical monocyte",
             "Neutrophil",
             "Immature neutrophil",
             "mDC",
             "pDC",
             "NK",
             "T",
             "Activated T",
             "B",
             "Platelet",
             "mixed4",
             "mixed6")) {
  df <- subset(exp2_blood@meta.data, celltype==ct)
  expr[[ct]] <- as.data.frame(subset(df, nCount_RNA>median(df$nCount_RNA) & nCount_RNA < quantile(df$nCount_RNA)[[4]]+3*IQR(df$nCount_RNA))$cell_id)
  colnames(expr[[ct]]) <- "cell_id"
}  

for (ct in c("mixed3", "mixed5")) {
  df <- subset(exp2_blood@meta.data, celltype==ct)
  expr[[ct]] <- as.data.frame(subset(df, nCount_RNA>quantile(subset(exp2_blood@meta.data, celltype=="Neutrophils")$nCount_RNA)[[2]])$cell_id)
  colnames(expr[[ct]]) <- "cell_id"
}  



cells_to_keep <- do.call(rbind,expr)
write.table(cells_to_keep, "./exp2_blood_cells_to_keep.txt", sep="\t", quote=FALSE)

write.table(exp2_blood@meta.data, "./exp2_blood_combined_250_metadata.txt", sep="\t", quote=FALSE)

