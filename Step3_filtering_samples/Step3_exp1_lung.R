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
exp1_lung <- readRDS("./seu_lung_exp1_combined_250.rds")

exp1_lung = NormalizeData(exp1_lung)
exp1_lung = FindVariableFeatures(exp1_lung, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(exp1_lung)
exp1_lung <- ScaleData(exp1_lung, features = all.genes)
exp1_lung <- RunPCA(exp1_lung, features = VariableFeatures(object = exp1_lung))
pdf("exp1_lung_combined_250_elbow.pdf")
ElbowPlot(exp1_lung)
dev.off()
exp1_lung <- RunUMAP(exp1_lung, dims = 1:19)
pdf("exp1_lung_combined_250_UMAPPlot.pdf")
UMAPPlot(object=exp1_lung)
dev.off()
exp1_lung <- FindNeighbors(exp1_lung, dims = 1:19)
exp1_lung <- FindClusters(exp1_lung, resolution = 0.9)
pdf("exp1_lung_combined_250_clusters.pdf")
DimPlot(exp1_lung, reduction = "umap", label=TRUE)
dev.off()

exp1_lung@meta.data = cbind(exp1_lung@meta.data, exp1_lung@reductions$umap@cell.embeddings)

exp1_lung@meta.data$treatment <- gsub("_[1-4].*","",exp1_lung@meta.data$orig.ident)
exp1_lung@meta.data$hamster <- gsub(".*_([1-4].*)","Ha\\1",exp1_lung@meta.data$orig.ident)


ggplot()+geom_violin(data=exp1_lung@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)))
ggsave("exp1_lung_violin_nUMI_log2_in_clusters.pdf", width=12, height=7)
ggplot()+geom_violin(data=exp1_lung@meta.data, aes(x=seurat_clusters, y=log2(nFeature_RNA)))
ggsave("exp1_lung_violin_nGene_log2_in_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_lung@meta.data, aes(x=seurat_clusters, y=log2(nFeature_RNA)), stat = "summary", fun = "mean")
ggsave("exp1_lung_barplot_log2_nGene_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_lung@meta.data, aes(x=seurat_clusters, y=log2(nFeature_RNA)), stat = "summary", fun = "median")
ggsave("exp1_lung_barplot_log2_nGene_median_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_lung@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)), stat = "summary", fun = "mean")
ggsave("exp1_lung_barplot_log2_nUMI_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_lung@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)), stat = "summary", fun = "median")
ggsave("exp1_lung_barplot_log2_nUMI_median_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_lung@meta.data, aes(x=seurat_clusters, y=nFeature_RNA), stat = "summary", fun = "mean")
ggsave("exp1_lung_barplot_nGene_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_lung@meta.data, aes(x=seurat_clusters, y=nFeature_RNA), stat = "summary", fun = "median")
ggsave("exp1_lung_barplot_nGene_median_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_lung@meta.data, aes(x=seurat_clusters, y=nCount_RNA), stat = "summary", fun = "mean")
ggsave("exp1_lung_barplot_nUMI_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=exp1_lung@meta.data, aes(x=seurat_clusters, y=nCount_RNA), stat = "summary", fun = "median")
ggsave("exp1_lung_barplot_nUMI_median_clusters.pdf", width=12, height=7)




###############
#cell type annotation
#LOC101833790 = Ccr5

markers <- c("Ace2", "Ackr1", "Acta2", "Adgre1", "Arg1", "Aspn", "C1qb", "Camp", "Ccr2", "LOC101833790", "Cd3e", "Cd4", "Cd79b", "Cd8a", "Cldn5", "Col1a2", "Cx3cr1", "Cxcr2", "Dcn", "Ednrb", "Fcgr4", "Flt3", "Foxj1", "Gzma", "Il7r", "Irf8", "Lamp3", "Marco", "Mrc1", "Ms4a1", "Nkg7", "Pdpn", "Ppbp", "Plvap", "Retn", "Rtkn2", "S100a8", "Siglecf", "Tagln", "Tcf4", "Treml4", "Pf4")
for (gene in markers) {
  df <- FetchData(exp1_lung, gene)
  gene <- str_replace_all(gene, "-", "_")
  colnames(df) <- gene
  df = cbind(df, exp1_lung@reductions$umap@cell.embeddings)
  ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
    geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
    scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    coord_fixed(ratio=1)+
    ggtitle(gene)
  ggsave(paste("exp1_lung_combined_250", gene, "pdf", sep="."), useDingbats=FALSE)
}

avg <- AverageExpression(exp1_lung, assays=c("RNA"), features = markers, return.seurat = T, slot="data") 
avg2 <- AverageExpression(exp1_lung, assays=c("RNA"), features = markers, return.seurat = F, slot="data")
avg2 <- as.matrix(avg2$RNA)
hc <- hclust(dist(avg2))
gene.order <- row.names(avg2)[order.hclust(hc)]
DoHeatmap(avg, size=5, features=gene.order)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("exp1_lung_cluster_heatmap.pdf")


#Myofibroblasts, fibroblasts, smooth muscle
#LOC101844074 is Cox4i2 (smooth muscle marker)
musclemarkers <- c("LOC101844074", "Aspn", "Dcn", "Acta2", "Tagln", "Wif1", "Fgf18", "Col1a2", "Bsg", "Cnn1", "Myh11", "Actg2", "Gpc3", "Apoe", "Serpinf1", "Gpc3")
avg <- AverageExpression(exp1_lung, assays=c("RNA"), features = musclemarkers, return.seurat = T, slot="data") 
avg2 <- AverageExpression(exp1_lung, assays=c("RNA"), features = musclemarkers, return.seurat = F, slot="data")
avg2 <- as.matrix(avg2$RNA)
hc <- hclust(dist(avg2))
gene.order <- row.names(avg2)[order.hclust(hc)]
DoHeatmap(avg, size=5, features=gene.order)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("exp1_lung_combined_250_musclemarkers_heatmap.pdf")


#Include doublet information by reading int the DoubSing file and attach to meta data
exp1_lung@meta.data$cell_id <- rownames(exp1_lung@meta.data)
DoubSing <- read.table("./DoubSing_lung_exp1.txt", header=TRUE)
exp1_lung@meta.data <- merge(exp1_lung@meta.data, DoubSing, by="cell_id")
rownames(exp1_lung@meta.data) <- exp1_lung@meta.data$cell_id

#Find the doublets as in the Nat Comms paper
genes <- c("Cldn5", "Plvap", "Gzma", "Cd8a", "Nkg7", "Cd3e", "Cd3d", "Cd3g", "Cd4", "Ifng", "Ccl5", "Gzmk")
df <- FetchData(exp1_lung, genes)
df <- cbind(df, exp1_lung@meta.data)
exp1_lung@meta.data$endocytotox <- ifelse((df$Cldn5 > 0 | df$Plvap > 0) & (df$Gzma > 0 | df$Cd8a > 0), "doublet", "singlet")


Idents(exp1_lung) <- exp1_lung@meta.data$seurat_clusters

#This is the annotation for the object with more than 250 genes
#pDC are likely mixed into the myDC cluster (24)
exp1_lung <- RenameIdents(exp1_lung, 
                       '0'='Neutrophils',
                       '1'='MonocyticMacrophages',
                       '2'='AT2',
                       '3'='MonocyticMacrophages',
                       '4'='Endothelial',
                       '5'='AlveolarMacrophages',
                       '6'='Bcells',
                       '7'='Treml4+Macrophages',
                       '8'='TNKcells',
                       '9'='AlveolarMacrophages',
                       '10'='Neutrophils',
                       '11'='AT1',
                       '12'='mixed1',
                       '13'='InterstitialMacrophages',
                       '14'='Platelets',
                       '15'='Endothelial',
                       '16'='Myofibroblast',
                       '17'='mixed2',
                       '18'='Treml4+Macrophages',
                       '19'='Endothelial',
                       '20'='Treml4+Macrophages',
                       '21'='mixed3',
                       '22'='Platelets',
                       '23'='TNKcells',
                       '24'='Neutrophils',
                       '25'='pDC',
                       '26'='Fibroblasts',
                       '27'='AlveolarMacrophages',
                       '28'='Bcells',
                       '29'='AlveolarMacrophages',
                       '30'='Bcells',
                       '31'='Ciliated',
                       '32'='InterstitialMacrophages',
                       '33'='Bcells',
                       '34'='AlveolarMacrophages',
                       '35'='Bcells')
exp1_lung@meta.data$celltype <- Idents(exp1_lung)



#save at this point
saveRDS(exp1_lung, "./exp1_lung_combined_250_ann.rds")



#Order cell types
the_celltypes = c("AlveolarMacrophages",
                  "InterstitialMacrophages",
                  "MonocyticMacrophages",
                  "Treml4+Macrophages",
                  "Neutrophils",
                  "MyeloidDendritic",
                  "pDC",
                  "TNKcells",
                  "Bcells",
                  "AT1",
                  "AT2",
                  "Fibroblasts",
                  "Ciliated",
                  "Endothelial",
                  "Myofibroblast",
                  "SmoothMuscle",
                  "Platelets",
                  "mixed1",
                  "mixed2",
                  "mixed3",
                  "mixed4",
                  "mixed5",
                  "mixed6",
                  "mixed7",
                  "mixed8",
                  "mixed9",
                  "mixed10"
)

celltypecolors = c(
  "AlveolarMacrophages" = "#DFACC4",
  "InterstitialMacrophages" = "#B97C9D",
  "MonocyticMacrophages" = "#B7245C",
  "Treml4+Macrophages" = "#3E2F5B",
  "Neutrophils" = "#0081AF",
  "MyeloidDendritic" = "#4F6D7A",
  "pDC" = "#7C6A0A",
  "TNKcells" = "#368F8B",
  "Bcells" = "#62C370",
  "AT1" = "#F7C548",
  "AT2" = "#F97E44",
  "Fibroblasts" = "#B2675E",
  "Ciliated" = "#FB3640",
  "Endothelial" = "#0D3B66",
  "Myofibroblast" = "#C4A381",
  "SmoothMuscle" = "#644536",
  "Platelets" = "#BF3E00",
  "mixed1" = "#CAD2C5",
  "mixed1" = "#CAD2C5",
  "mixed2" = "#CAD2C5",
  "mixed3" = "#CAD2C5",
  "mixed4" = "#CAD2C5",
  "mixed5" = "#CAD2C5",
  "mixed6" = "#CAD2C5",
  "mixed7" = "#CAD2C5",
  "mixed8" = "#CAD2C5",
  "mixed9" = "#CAD2C5",
  "mixed10" = "#CAD2C5"
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
means <- exp1_lung@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(mean_U1 = mean(UMAP_1), mean_U2 = mean(UMAP_2))

ggplot()+
  geom_point(data=exp1_lung@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=celltype), size=0.5, shape=16)+
  geom_text(data=means, aes(x=mean_U1, y=mean_U2, label=celltype), size=3)+
  coord_fixed(ratio=1)+
  ggtitle("celltype")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  scale_colour_manual(values = celltypecolors)
ggsave("exp1_lung_combined_250_celltypes.pdf", useDingbats=FALSE, width=14, height=7)



#doublets by celltype and treatment
df1 <- cbind.data.frame(exp1_lung@meta.data$DoubSing,
                        exp1_lung@meta.data$endocytotox,
                        exp1_lung@meta.data$celltype,
                        exp1_lung@meta.data$treatment,
                        exp1_lung@meta.data$hamster)
colnames(df1) <- c("DoubSing", "endocytotox", "celltype", "treatment", "hamster")

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
ggsave("exp1_lung_combined_250_DoubSing.pdf", useDingbats=FALSE, width=14, height=7)

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
ggsave("exp1_lung_combined_250_endocytotox.pdf", useDingbats=FALSE, width=14, height=7)


ggplot()+geom_violin(data=exp1_lung@meta.data %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes)), aes(x=celltype, y=log2(nCount_RNA), fill=celltype))+scale_fill_manual(values = celltypecolors)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10))
ggsave("exp1_lung_combined_250_violin_nUMI_log2_in_celltype.pdf", width=12, height=7)
ggplot()+geom_violin(data=exp1_lung@meta.data %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes)), aes(x=celltype, y=log2(nFeature_RNA), fill=celltype))+scale_fill_manual(values = celltypecolors)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10))
ggsave("exp1_lung_combined_250_violin_nGene_log2_in_celltype.pdf", width=12, height=7)



#Check AT1 cells, which have two clusters and very broad nUMI distribution

seu_AT1 <- subset(exp1_lung, subset = (celltype == "AT1"))
seu_AT1 <- RunPCA(seu_AT1, verbose = FALSE)
seu_AT1 <- RunUMAP(seu_AT1, dims = 1:30)
seu_AT1 <- FindNeighbors(seu_AT1, dims = 1:30)
seu_AT1 <- FindClusters(seu_AT1, resolution = 0.9)
seu_AT1@meta.data$UMAP_1 <- NULL
seu_AT1@meta.data$UMAP_2 <- NULL
seu_AT1@meta.data = cbind(seu_AT1@meta.data, seu_AT1@reductions$umap@cell.embeddings)
Idents(seu_AT1) <- seu_AT1@meta.data$seurat_clusters
pdf("seu_AT1_UMAPPlot.pdf")
UMAPPlot(seu_AT1)
dev.off()
pdf("seu_AT1_clusters.pdf")
DimPlot(seu_AT1, reduction = "umap", label=TRUE)
dev.off()

cluster.markers.neu <- FindAllMarkers(seu_AT1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers.neu %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
top10 <- cluster.markers.neu %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("seu_AT1_clustermarkers.pdf", height=12, width=7)
DoHeatmap(seu_AT1, features = top10$gene) + NoLegend()
dev.off()
ggplot()+geom_violin(data=seu_AT1@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)))
ggsave("AT1_violin_nUMI_log2_in_clusters.pdf", width=12, height=7)
dev.off()
ggplot()+geom_violin(data=seu_AT1@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA), fill=timepoint))
ggsave("AT1_violin_nUMI_log2_in_clusters_timepoint.pdf", width=12, height=7)


seu_AT1[["percent.mt"]] <- PercentageFeatureSet(object = seu_AT1, features=c("ATP6", "COX1", "COX2", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"))
pdf("AT1_violin_percent.mt.pdf", width=12, height=7)
VlnPlot(seu_AT1, features = "percent.mt", pt.size = 0, y.max = 20)
dev.off()

seu_AT2 <- subset(exp1_lung, subset = (celltype == "AT2"))
seu_AT2 <- RunPCA(seu_AT2, verbose = FALSE)
seu_AT2 <- RunUMAP(seu_AT2, dims = 1:30)
seu_AT2 <- FindNeighbors(seu_AT2, dims = 1:30)
seu_AT2 <- FindClusters(seu_AT2, resolution = 0.9)
seu_AT2@meta.data$UMAP_1 <- NULL
seu_AT2@meta.data$UMAP_2 <- NULL
seu_AT2@meta.data = cbind(seu_AT2@meta.data, seu_AT2@reductions$umap@cell.embeddings)
Idents(seu_AT2) <- seu_AT2@meta.data$seurat_clusters
pdf("seu_AT2_UMAPPlot.pdf")
UMAPPlot(seu_AT2)
dev.off()
pdf("seu_AT2_clusters.pdf")
DimPlot(seu_AT2, reduction = "umap", label=TRUE)
dev.off()

cluster.markers.neu <- FindAllMarkers(seu_AT2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers.neu %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
top10 <- cluster.markers.neu %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("seu_AT2_clustermarkers.pdf", height=12, width=7)
DoHeatmap(seu_AT2, features = top10$gene) + NoLegend()
dev.off()
ggplot()+geom_violin(data=seu_AT2@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)))
ggsave("AT2_violin_nUMI_log2_in_clusters.pdf", width=12, height=7)
dev.off()
ggplot()+geom_violin(data=seu_AT2@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA), fill=timepoint))
ggsave("AT2_violin_nUMI_log2_in_clusters_timepoint.pdf", width=12, height=7)


seu_AT2[["percent.mt"]] <- PercentageFeatureSet(object = seu_AT2, features=c("ATP6", "COX1", "COX2", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"))
pdf("AT2_violin_percent.mt.pdf", width=12, height=7)
VlnPlot(seu_AT2, features = "percent.mt", pt.size = 0, y.max = 20)
dev.off()


#Endothelial subclustering
seu_endo <- subset(exp1_lung, subset = (celltype == "Endothelial"))
seu_endo <- RunPCA(seu_endo, verbose = FALSE)
seu_endo <- RunUMAP(seu_endo, dims = 1:30)
seu_endo <- FindNeighbors(seu_endo, dims = 1:30)
seu_endo <- FindClusters(seu_endo, resolution = 0.9)
seu_endo@meta.data$UMAP_1 <- NULL
seu_endo@meta.data$UMAP_2 <- NULL
seu_endo@meta.data = cbind(seu_endo@meta.data, seu_endo@reductions$umap@cell.embeddings)
Idents(seu_endo) <- seu_endo@meta.data$seurat_clusters
pdf("seu_endo_UMAPPlot.pdf")
UMAPPlot(seu_endo)
dev.off()
pdf("seu_endo_clusters.pdf")
DimPlot(seu_endo, reduction = "umap", label=TRUE)
dev.off()

cluster.markers.neu <- FindAllMarkers(seu_endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers.neu %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
top10 <- cluster.markers.neu %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("seu_endo_clustermarkers.pdf", height=12, width=7)
DoHeatmap(seu_endo, features = top10$gene) + NoLegend()
dev.off()
ggplot()+geom_violin(data=seu_endo@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)))
ggsave("endo_violin_nUMI_log2_in_clusters.pdf", width=12, height=7)
dev.off()
ggplot()+geom_violin(data=seu_endo@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA), fill=timepoint))
ggsave("endo_violin_nUMI_log2_in_clusters_timepoint.pdf", width=12, height=7)


seu_endo[["percent.mt"]] <- PercentageFeatureSet(object = seu_endo, features=c("ATP6", "COX1", "COX2", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"))
pdf("endo_violin_percent.mt.pdf", width=12, height=7)
VlnPlot(seu_endo, features = "percent.mt", pt.size = 0, y.max = 20)
dev.off()

for (gene in c("Ednrb", "Ackr1", "Ccl21", "Gja5", "Dkk2", "Pdpn")) {
  df <- FetchData(seu_endo, gene)
  gene <- str_replace_all(gene, "-", "_")
  colnames(df) <- gene
  df = cbind(df, seu_endo@reductions$umap@cell.embeddings)
  ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
    geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
    scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    coord_fixed(ratio=1)+
    ggtitle(gene)
  ggsave(paste("seu_endo", gene, "pdf", sep="."), useDingbats=FALSE)
}

ggplot()+
  geom_point(data=seu_endo@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=log2(nCount_RNA)), shape=16, size=0.5)+
  scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  coord_fixed(ratio=1)+
  ggtitle("nCountRNA_UMAP")
ggsave(paste("seu_endo_nCountRNA_UMAP", "pdf", sep="."), useDingbats=FALSE)

#Check median of Ccl21/Pdpn expressing cells (i.e. lymphatic endothelium)
#Set a more stringent threshold for defining a cell as lymphatic, and use this to separate the endothelial cells when filtering by median
df <- FetchData(seu_endo, c("Ccl21", "Pdpn"))
df <- df %>% mutate(lymphatic = if_else((Ccl21 > 2 | Pdpn > 1), "TRUE", "FALSE"))
df = cbind(df, seu_endo@meta.data)
median(subset(df, lymphatic==TRUE)$nCount_RNA)
median(subset(df, lymphatic==FALSE)$nCount_RNA)
median(subset(df, lymphatic==TRUE)$nFeature_RNA)
median(subset(df, lymphatic==FALSE)$nFeature_RNA)
median(subset(df, lymphatic==TRUE)$percent.mt)
median(subset(df, lymphatic==FALSE)$percent.mt)
df2 <- cbind.data.frame(df$nCount_RNA, df$Ccl21)
colnames(df2) <- c("nCount_RNA", "expression")
df2$gene <- "Ccl21"
df3 <- cbind.data.frame(df$nCount_RNA, df$Pdpn)
colnames(df3) <- c("nCount_RNA", "expression")
df3$gene <- "Pdpn"
df4 <- rbind.data.frame(df2, df3)
df4 <- subset(df4, expression>0)
ggplot()+geom_violin(data=df4, aes(x=gene, y=expression))
ggsave("endo_violin_Ccl21_Pdpn_expr_lymphatic.pdf", width=12, height=7)
ggplot()+geom_violin(data=df, aes(x=seurat_clusters, y=log2(nCount_RNA), colour=lymphatic))
ggsave("endo_violin_nUMI_log2_in_clusters_lymphatic.pdf", width=12, height=7)
ggplot()+geom_bar(data=df, aes(x=lymphatic, y=nCount_RNA), position=position_dodge(), stat = "summary", fun = "median")
ggsave("barplot_log2_nUMI_endo_lymph_median.pdf", width=12, height=7)
ggplot()+geom_bar(data=df, aes(x=seurat_clusters, y=nCount_RNA, colour=lymphatic), position=position_dodge(), stat = "summary", fun = "median")
ggsave("barplot_log2_nUMI_endo_lymph_median_clusters.pdf", width=12, height=7)
ggplot()+geom_violin(data=df, aes(x=seurat_clusters, y=log2(nFeature_RNA), colour=lymphatic))
ggsave("endo_violin_nUMI_log2_in_clusters_lymphatic.pdf", width=12, height=7)


#T/NK cells

seu_TNK <- subset(exp1_lung, subset = (celltype == "TNKcells"))
seu_TNK <- RunPCA(seu_TNK, verbose = FALSE)
seu_TNK <- RunUMAP(seu_TNK, dims = 1:30)
seu_TNK <- FindNeighbors(seu_TNK, dims = 1:30)
seu_TNK <- FindClusters(seu_TNK, resolution = 0.9)
seu_TNK@meta.data$UMAP_1 <- NULL
seu_TNK@meta.data$UMAP_2 <- NULL
seu_TNK@meta.data = cbind(seu_TNK@meta.data, seu_TNK@reductions$umap@cell.embeddings)
Idents(seu_TNK) <- seu_TNK@meta.data$seurat_clusters
pdf("seu_TNK_UMAPPlot.pdf")
UMAPPlot(seu_TNK)
dev.off()
pdf("seu_TNK_clusters.pdf")
DimPlot(seu_TNK, reduction = "umap", label=TRUE)
dev.off()

cluster.markers.neu <- FindAllMarkers(seu_TNK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers.neu %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
top10 <- cluster.markers.neu %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("seu_TNK_clustermarkers.pdf", height=12, width=7)
DoHeatmap(seu_TNK, features = top10$gene) + NoLegend()
dev.off()
ggplot()+geom_violin(data=seu_TNK@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)))
ggsave("TNK_violin_nUMI_log2_in_clusters.pdf", width=12, height=7)
dev.off()
ggplot()+geom_violin(data=seu_TNK@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA), fill=timepoint))
ggsave("TNK_violin_nUMI_log2_in_clusters_timepoint.pdf", width=12, height=7)


seu_TNK[["percent.mt"]] <- PercentageFeatureSet(object = seu_TNK, features=c("ATP6", "COX1", "COX2", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"))
pdf("TNK_violin_percent.mt.pdf", width=12, height=7)
VlnPlot(seu_TNK, features = "percent.mt", pt.size = 0, y.max = 20)
dev.off()

for (gene in c("Cd4", "Cd8a", "Nkg7", "Cd3e", "Mki67")) {
  df <- FetchData(seu_TNK, gene)
  gene <- str_replace_all(gene, "-", "_")
  colnames(df) <- gene
  df = cbind(df, seu_TNK@reductions$umap@cell.embeddings)
  ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
    geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
    scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    coord_fixed(ratio=1)+
    ggtitle(gene)
  ggsave(paste("seu_TNK", gene, "pdf", sep="."), useDingbats=FALSE)
}


#Check median of Mki67/Top2a expressing cells (i.e. dividing)
#Set a more stringent threshold for defining a cell as dividing, and use this to separate the T/NK cells when filtering by median
df <- FetchData(seu_TNK, c("Mki67", "Top2a"))
df = cbind(df, seu_TNK@meta.data)
df2 <- cbind.data.frame(df$nCount_RNA, df$Mki67)
colnames(df2) <- c("nCount_RNA", "expression")
df2$gene <- "Mki67"
df3 <- cbind.data.frame(df$nCount_RNA, df$Top2a)
colnames(df3) <- c("nCount_RNA", "expression")
df3$gene <- "Top2a"
df4 <- rbind.data.frame(df2, df3)
df4 <- subset(df4, expression>0)
ggplot()+geom_violin(data=df4, aes(x=gene, y=expression))
ggsave("TNK_violin_Mki67_Top2a_expr.pdf", width=12, height=7)

df <- df %>% mutate(dividing = if_else((Mki67 > 0.5 | Top2a > 1), "TRUE", "FALSE"))
median(subset(df, dividing==TRUE)$nCount_RNA)
median(subset(df, dividing==FALSE)$nCount_RNA)
median(subset(df, dividing==TRUE)$nFeature_RNA)
median(subset(df, dividing==FALSE)$nFeature_RNA)
ggplot()+geom_violin(data=df, aes(x=seurat_clusters, y=log2(nCount_RNA), colour=dividing))
ggsave("TNK_violin_nUMI_log2_in_clusters_dividing.pdf", width=12, height=7)




#############
#final filtering

#Most cell types get standard filtering, mixed without outliers, endothelial and TNK with a subset differentiation
#mixed1-3 are doublet cells, i.e. no outlier removal
#smooth muscle and myeloid dendritic not defined (likely in other celltypes)

expr <- list()
for (ct in c("Bcells", "Neutrophils", "AT2", "AlveolarMacrophages", "Myofibroblast", "MonocyticMacrophages", "AT1", "Platelets", "InterstitialMacrophages", "Ciliated", "Treml4+Macrophages", "pDC", "Fibroblasts")) {
  df <- subset(exp1_lung@meta.data, celltype==ct)
  expr[[ct]] <- as.data.frame(subset(df, nCount_RNA>quantile(df$nCount_RNA)[[2]] & nCount_RNA < quantile(df$nCount_RNA)[[4]]+3*IQR(df$nCount_RNA))$cell_id)
  colnames(expr[[ct]]) <- "cell_id"
}  

for (ct in c("mixed1", "mixed2", "mixed3")) {
  df <- subset(exp1_lung@meta.data, celltype==ct)
  expr[[ct]] <- as.data.frame(subset(df, nCount_RNA>quantile(subset(exp1_lung@meta.data, celltype=="Neutrophils")$nCount_RNA)[[2]])$cell_id)
  colnames(expr[[ct]]) <- "cell_id"
}  


#For endothelial cells, split up in lymphatic/others
df <- FetchData(exp1_lung, c("Ccl21", "Pdpn"))
df <- df %>% mutate(lymphatic = if_else((Ccl21 > 2 | Pdpn > 1), "TRUE", "FALSE"))
df = cbind(df, exp1_lung@meta.data)
df <- subset(df, celltype=="Endothelial")
df1 <- subset(df, lymphatic == TRUE)
df2 <- subset(df1, nCount_RNA>quantile(df1$nCount_RNA)[[2]] & nCount_RNA < quantile(df1$nCount_RNA)[[4]]+3*IQR(df1$nCount_RNA))
df3 <- subset(df, lymphatic == FALSE)
df4 <- subset(df3, nCount_RNA>quantile(df3$nCount_RNA)[[2]] & nCount_RNA < quantile(df3$nCount_RNA)[[4]]+3*IQR(df3$nCount_RNA))
expr[["Endothelial"]] <- as.data.frame(rbind.data.frame(df2, df4)$cell_id)
colnames(expr[["Endothelial"]]) <- "cell_id"

#For T/NK cells, split up in dividing/others
df <- FetchData(exp1_lung, c("Mki67", "Top2a"))
df <- df %>% mutate(dividing = if_else((Mki67 > 0.5 | Top2a > 1), "TRUE", "FALSE"))
df = cbind(df, exp1_lung@meta.data)
df <- subset(df, celltype=="TNKcells")
df1 <- subset(df, dividing == TRUE)
df2 <- subset(df1, nCount_RNA>quantile(df1$nCount_RNA)[[2]] & nCount_RNA < quantile(df1$nCount_RNA)[[4]]+3*IQR(df1$nCount_RNA))
df3 <- subset(df, dividing == FALSE)
df4 <- subset(df3, nCount_RNA>quantile(df3$nCount_RNA)[[2]] & nCount_RNA < quantile(df3$nCount_RNA)[[4]]+3*IQR(df3$nCount_RNA))
expr[["TNKcells"]] <- as.data.frame(rbind.data.frame(df2, df4)$cell_id)
colnames(expr[["TNKcells"]]) <- "cell_id"


cells_to_keep <- do.call(rbind,expr)
write.table(cells_to_keep, "./exp1_lung_cells_to_keep.txt", sep="\t", quote=FALSE)

write.table(exp1_lung@meta.data, "./exp1_lung_combined_250_metadata.txt", sep="\t", quote=FALSE)

