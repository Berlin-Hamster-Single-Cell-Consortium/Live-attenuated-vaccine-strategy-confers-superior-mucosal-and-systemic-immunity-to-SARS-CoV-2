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
source("~/Documents/Largescale-data/notes and scripts/smooth_DimPlot.R")
library(akima)
library(pheatmap)
#see http://www.cookbook-r.com/Graphs/ Plotting_means_and_error_bars_(ggplot2)
#source("~/Documents/Largescale-data/notes and scripts/summarySE.R")
source("/fast/AG_Landthaler/scripts/summarySE.R")
library(cowplot)
library(SeuratObject)
library(lme4)
library(ComplexHeatmap)
library(patchwork)
library(forcats)
library(DoubletFinder)

#Fig. 3L, S9

exp2_blood  <- readRDS("./exp2_blood_combined_integrated_annotated.rds")

treatment_levels <- c("att2x", "mRNAatt", "mRNA2x", "adeno2x", "aaUntr")

#subset the B cells from the blood Seurat objects, or download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200596
seu_B_blood <- subset(exp2_blood, subset = (celltype %in% c("B")))
DefaultAssay(seu_B_Blood) <- "integrated"                  
seu_B_blood <- RunPCA(seu_B_Blood, verbose = FALSE)
seu_B_Blood <- RunUMAP(seu_B_Blood, dims = 1:30)
seu_B_Blood <- FindNeighbors(seu_B_Blood, dims = 1:30)
seu_B_Blood <- FindClusters(seu_B_Blood, resolution = 1)
seu_B_Blood@meta.data$UMAP_1 <- NULL
seu_B_Blood@meta.data$UMAP_2 <- NULL
seu_B_Blood@meta.data = cbind(seu_B_Blood@meta.data, seu_B_Blood@reductions$umap@cell.embeddings)
DefaultAssay(seu_B_Blood) <- "SCT"

UMAPPlot(seu_B_blood, label=TRUE)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("Clusters")
ggsave(paste("exp2_blood_B", "clusters", "pdf", sep="."), useDingbats=FALSE)

#Broad B cell marker list
markers <- c("Cd79b", "Cd79a", "Cd19", "Sec11c", "Jchain", "Cd24", "Irf4", "Tox", "Mki67", "Aim2", "Tnfrsf17", "Cd80", "Runx2", "Ikzf2", "Rora", "Zbtb16", "Fkbp11", "Top2a", "Foxm1", "Pcdh9", "Gdpd5", "Prdm1", "Spry1", "Trerf1", "Pbx3", "Bcl6", "E2f1", "Stat2", "Spib", "Aicda", "Gabarapl1", "Id2", "Samsn1", "Xbp1", "Batf", "Bach2", "Zhx2", "Tex9", "Znf32", "Sox5", "Sox4", "Tnfrsf13b", "Irf7", "Cd38", "Gtf2i", "Ciita", "Cd27", "Cd83", "Pigr", "Tcf4", "Cd22", "Pou2af1", "Cr2", "Irf8", "Tnfrsf13c", "Ms4a1", "Pax5", "Ebf1")

#Kassambara et al. 2021, Fig. 5A
markers <- c("Batf2", "Irf2", "Dnmt3b", "Zscan20", "Prdm1", "Preb", "Homez", "Idh1", "Nr1h3", "Eya2", "Znf691", "Snai3", "Znf2", "Maf", "Zfp64", "Xbp1", "Stat1", "Fli1", "Bhlhe41", "Crebl2", "Gfi1", "Cxxc1", "Znf583", "Znf133", "Idh2", "Mlx", "E2f3", "Arid3a", "Brd8", "Eya3", "Tet1", "Alkbh3", "Zbtb38", "Irf4", "Irf7", "Stat2", "Prdm15", "Znf710", "Lcorl", "Hcfc2", "Drap1")

#Kassambara et al. 2021, Fig. 5B
markers <- c("Aicda", "Myb", "Batf3", "Foxm1", "Arntl2", "Pcna", "Tcf19", "Batf", "Suv39h2", "Rad54l", "Id2", "E2f1", "Tfdp1", "Rfx2", "Mllt3", "Atf5", "Dnmt1", "Hells", "Prmt1", "Ttf2", "Klf5", "Mybl2", "Tp53", "Prmt3", "Ezh2", "Suv39h1", "Mbd2", "Nfe2l3", "Sox4")

#Heatmap without clustering 
avg <- AverageExpression(seu_B_blood, assays=c("SCT"), features = markers_red, return.seurat = T, slot="data") 
DoHeatmap(avg, size=5, features=markers, raster = FALSE, draw.lines = FALSE)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")
ggsave("exp2_blood_integrated_B_cluster_markers_heatmap.pdf", width=8, height=7)

#Heatmap with clustering 
avg <- AverageExpression(seu_B_blood, assays=c("SCT"), features = markers, return.seurat = T, slot="data") 
avg2 <- AverageExpression(seu_B_blood, assays=c("SCT"), features = markers, return.seurat = F, slot="data")
avg2 <- as.matrix(avg2$SCT)
hc <- hclust(dist(avg2))
gene.order <- row.names(avg2)[order.hclust(hc)]
DoHeatmap(avg, size=5, features=gene.order, raster = FALSE, draw.lines = FALSE)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")
ggsave("exp2_blood_integrated_B_cluster_markers_heatmap.pdf", width=8, height=7)

#Look at expression of specific genes per cluster
expr <- list()
markers <- c("Aicda", "Bach2", "Batf", "Cd27", "Irf4", "Pax5", "Tnfrsf13b", "Tnfrsf13c")
clusters_to_check <- unique(seu_B_blood@meta.data$seurat_clusters)

for (gene in markers) {
  df1 <- cbind.data.frame(df <- FetchData(seu_B_blood, gene),
                          seu_B_blood@meta.data$seurat_clusters,
                          seu_B_blood@meta.data$treatment,
                          seu_B_blood@meta.data$hamster)
  colnames(df1) <- c("gene", "clusters", "treatment", "hamster")
  df1 <- subset (df1, clusters %in% clusters_to_check)
  for (the_cluster in clusters_to_check) {
    for (the_treatment in unique(seu_B_blood@meta.data$treatment)) {
      a <- dim(subset(df1, clusters==the_cluster & treatment == the_treatment))[[1]]
      if (a == 0) {
        expr[[paste(gene, the_cluster, the_treatment, sep="_")]] <- "0_0"
      } else {
        b <- dim(subset(df1, clusters==the_cluster & treatment == the_treatment & gene >0))[[1]]
        if (b == 0) {
          c <- 0
        } else {
          c <- mean(subset(df1, clusters==the_cluster & treatment == the_treatment & gene >0)$gene)
        }
        expr[[paste(gene, the_cluster, the_treatment, sep="_")]] <- paste(b/a, c, sep="_")
      }
    }
  }
}
positives <- as.data.frame(do.call(rbind,expr))
positives$name <- rownames(positives)
positives <- positives %>% separate("name", into=c("gene", "cluster", "treatment"),  sep="_") %>% separate("V1", into=c("pos", "expr"),  sep="_")
positives$pos <- as.numeric(positives$pos)
positives$expr <- as.numeric(positives$expr)
positives$cluster <- as.integer(positives$cluster)

ggplot(positives,
       aes(x=gene,y=factor(treatment, rev(treatment_levels)),size=pos,color=expr)) +
  geom_point(shape=16) +
  scale_size_continuous(name='fraction positive',
                        breaks=c(0.25,0.5,0.75),
                        labels=c("0.25","0.5", "0.75"),
                        limits=c(0,1)) +
  scale_color_gradient(low='gray80',high='red3') +
  facet_wrap(~cluster, ncol=1) +
  coord_cartesian(clip = 'off') +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
  labs(x='',y='')

ggsave("exp2_blood_B_markersFig3_expr_in_clusters.pdf", width=3.5, height = 4, units="in", useDingbats=FALSE)


