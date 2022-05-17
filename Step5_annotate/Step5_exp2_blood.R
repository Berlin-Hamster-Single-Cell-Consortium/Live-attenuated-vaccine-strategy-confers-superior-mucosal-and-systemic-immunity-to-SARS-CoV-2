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


#Define cell type colors



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



#read integrated object from cluster
exp2_blood <- readRDS("./seu_exp2_blood_new_combined_integrated.rds")
exp2_blood@meta.data = cbind(exp2_blood@meta.data, exp2_blood@reductions$umap@cell.embeddings)
exp2_blood@meta.data$treatment <- gsub("_[1-4].*","",exp2_blood@meta.data$orig.ident)
exp2_blood@meta.data$hamster <- gsub(".*_([1-4].*)","Ha\\1",exp2_blood@meta.data$orig.ident)

DefaultAssay(exp2_blood) <- 'RNA'
SCoV2_rawcounts <- FetchData(exp2_blood, grep("SCoV2", exp2_blood@assays$RNA@counts@Dimnames[[1]], value="TRUE"), slot="counts")
exp2_blood@meta.data$SCoV2_sum <- rowSums(SCoV2_rawcounts)
SCoV2_rawcounts$SCoV2sum <- rowSums(SCoV2_rawcounts)
exp2_blood@meta.data$SCoV2_load <- SCoV2_rawcounts$SCoV2sum/exp2_blood@meta.data$nCount_RNA*100
DefaultAssay(exp2_blood) <- 'SCT'

ggplot()+
  geom_point(data=exp2_blood@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=treatment), size=0.5, shape=16, alpha=0.5)+
  ggtitle("treatment")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("exp2_blood_integrated_treatment.pdf", useDingbats=FALSE)
ggplot()+
  geom_point(data=exp2_blood@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=hamster), size=0.2, shape=16, alpha=0.5)+
  ggtitle("hamster")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("exp2_blood_integrated_hamsters.pdf", useDingbats=FALSE)
UMAPPlot(exp2_blood, label=TRUE)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("Clusters")
ggsave("exp2_blood_integrated_clusters.pdf", useDingbats=FALSE)

#If granularity needs to be adjusted
DefaultAssay(exp2_blood) <- 'integrated'
exp2_blood <- FindClusters(exp2_blood, resolution = 0.7)
DefaultAssay(exp2_blood) <- 'SCT'



ggplot()+geom_point(data=subset(exp2_blood@meta.data, SCoV2_load<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", size=0.5, shape = 16)+
  geom_point(data=subset(exp2_blood@meta.data, SCoV2_load>0), aes(x=UMAP_1, y=UMAP_2, colour=log10(SCoV2_load)), size=0.5, shape = 16)+
  scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("log10 SCoV2 load")
ggsave(paste("exp2_blood_integrated_log10SCoV2_load", "pdf", sep="."), useDingbats=FALSE)

unfilt_metadata <- read.table("./exp2_blood_combined_250_metadata.txt", sep="\t")
exp2_blood@meta.data$cell_id <- rownames(exp2_blood@meta.data)
df <- dplyr::left_join(exp2_blood@meta.data, unfilt_metadata, by="cell_id")
exp2_blood@meta.data$DoubSing <- df$DoubSing

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
  ggsave(paste("exp2_blood_integrated", gene, "pdf", sep="."), useDingbats=FALSE)
}

avg <- AverageExpression(exp2_blood, assays=c("SCT"), features = markers, return.seurat = T, slot="data") 
avg2 <- AverageExpression(exp2_blood, assays=c("SCT"), features = markers, return.seurat = F, slot="data")
avg2 <- as.matrix(avg2$SCT)
hc <- hclust(dist(avg2))
gene.order <- row.names(avg2)[order.hclust(hc)]
DoHeatmap(avg, size=5, features=gene.order)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("exp2_blood_integrated_cluster_heatmap.pdf", width=12, height=7)

cluster.markers <- FindAllMarkers(exp2_blood, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("exp2_blood_integrated_cluster_markers_heatmap.pdf", width=24, height=16)
DoHeatmap(exp2_blood, features = top10$gene) + NoLegend()
dev.off()

#percent clusters per treatment
df1 <- cbind.data.frame(exp2_blood@meta.data$SCoV2_load,
                        exp2_blood@meta.data$seurat_clusters,
                        exp2_blood@meta.data$treatment,
                        exp2_blood@meta.data$hamster)
colnames(df1) <- c("SCoV2_load", "cluster", "treatment", "hamster")

a = df1 %>% group_by(treatment, hamster) %>% tally(name="tot") 
b = df1 %>% group_by(cluster, treatment, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('treatment', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
#write.table(c, "celltypefractions.tsv", sep="\t", quote=FALSE, row.names = FALSE)

tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("cluster", "treatment"))
ggplot(tgc, aes(x=cluster, y=fraction, fill=paste(cluster, treatment, sep="_")))+
  geom_bar(aes(colour=treatment), position=position_dodge(.8), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=cluster, y=fraction, fill=paste(cluster, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("celltypes per timepoint", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))

ggsave("exp2_blood_integrated_clusterpercentage_pertimepoint.pdf", useDingbats=FALSE)


Idents(exp2_blood) <- exp2_blood@meta.data$seurat_clusters

exp2_blood <- RenameIdents(exp2_blood, 
                          '0'='Neutrophil',
                          '1'='Immature neutrophil',
                          '2'='B',
                          '3'='T',
                          '4'='Platelet',
                          '5'='Classical monocyte',
                          '6'='B',
                          '7'='Neutrophil',
                          '8'='Classical monocyte',
                          '9'='B',
                          '10'='T',
                          '11'='Non-classical monocyte',
                          '12'='B',
                          '13'='mDC',
                          '14'='NK')
exp2_blood@meta.data$celltype <- Idents(exp2_blood)

saveRDS(exp2_blood, "./exp2_blood_combined_integrated_annotated.rds")


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
ggsave("exp2_blood_integrated_celltypes.pdf", useDingbats=FALSE, width=14, height=7)

#percent cell type per treatment
df1 <- cbind.data.frame(exp2_blood@meta.data$SCoV2_load,
                        exp2_blood@meta.data$celltype,
                        exp2_blood@meta.data$treatment,
                        exp2_blood@meta.data$hamster)
colnames(df1) <- c("SCoV2_load", "celltype", "treatment", "hamster")

#if not displaying all celltypes
#df1 <- subset(df1, celltype %in% c("Neutrophil", "Immature neutrophil", "T", "B"))

a = df1 %>% group_by(treatment, hamster) %>% tally(name="tot") 
b = df1 %>% group_by(celltype, treatment, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('treatment', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
write.table(c, "exp2_blood_celltypefractions.tsv", sep="\t", quote=FALSE, row.names = FALSE)

tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("celltype", "treatment"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))

ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")))+
  geom_bar(aes(colour=treatment), position=position_dodge(.9), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.9), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("celltypes per timepoint", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_colors)+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))

ggsave("exp2_blood_integrated_celltypepercentage_pertimepoint.pdf", useDingbats=FALSE)
ggsave("exp2_blood_integrated_celltypepercentage_subs_pertimepoint.pdf", useDingbats=FALSE)


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



###################
#Density plots
#smooth dimplots
cond1="aaUntr"
cond2="adeno2x"
smooth_DimPlot(subset(exp2_blood, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("exp2_blood_density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)
cond1="aaUntr"
cond2="att2x"
smooth_DimPlot(subset(exp2_blood, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("exp2_blood_density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)
cond1="aaUntr"
cond2="mRNA2x"
smooth_DimPlot(subset(exp2_blood, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("exp2_blood_density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)
cond1="aaUntr"
cond2="mRNAatt"
smooth_DimPlot(subset(exp2_blood, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("exp2_blood_density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)





