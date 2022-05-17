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

#Define cell type colors

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
                  "mixed10",
                  "mixed11",
                  "mixed12"
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
  "mixed10" = "#CAD2C5",
  "mixed11" = "#CAD2C5",
  "mixed12" = "#CAD2C5"
)


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



#read integrated object
exp2_lung <- readRDS("./seu_exp2_lung_new_combined_integrated.rds")
exp2_lung@meta.data = cbind(exp2_lung@meta.data, exp2_lung@reductions$umap@cell.embeddings)
exp2_lung@meta.data$treatment <- gsub("_[1-4].*","",exp2_lung@meta.data$orig.ident)
exp2_lung@meta.data$hamster <- gsub(".*_([1-4].*)","Ha\\1",exp2_lung@meta.data$orig.ident)

DefaultAssay(exp2_lung) <- 'RNA'
SCoV2_rawcounts <- FetchData(exp2_lung, grep("SCoV2", exp2_lung@assays$RNA@counts@Dimnames[[1]], value="TRUE"), slot="counts")
exp2_lung@meta.data$SCoV2_sum <- rowSums(SCoV2_rawcounts)
SCoV2_rawcounts$SCoV2sum <- rowSums(SCoV2_rawcounts)
exp2_lung@meta.data$SCoV2_load <- SCoV2_rawcounts$SCoV2sum/exp2_lung@meta.data$nCount_RNA*100
DefaultAssay(exp2_lung) <- 'SCT'

ggplot()+
  geom_point(data=exp2_lung@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=treatment), size=0.5, shape=16, alpha=0.5)+
  ggtitle("treatment")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("exp2_lung_integrated_treatment.pdf", useDingbats=FALSE)
ggplot()+
  geom_point(data=exp2_lung@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=hamster), size=0.2, shape=16, alpha=0.5)+
  ggtitle("hamster")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("exp2_lung_integrated_hamsters.pdf", useDingbats=FALSE)
UMAPPlot(exp2_lung, label=TRUE)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("Clusters")
ggsave("exp2_lung_integrated_clusters.pdf", useDingbats=FALSE)

#Adjust granularity which is needed to separate mDC from pDC
DefaultAssay(exp2_lung) <- 'integrated'
exp2_lung <- FindClusters(exp2_lung, resolution = 0.7)
DefaultAssay(exp2_lung) <- 'SCT'



ggplot()+geom_point(data=subset(exp2_lung@meta.data, SCoV2_load<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", size=0.5, shape = 16)+
  geom_point(data=subset(exp2_lung@meta.data, SCoV2_load>0), aes(x=UMAP_1, y=UMAP_2, colour=log10(SCoV2_load)), size=0.5, shape = 16)+
  scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("log10 SCoV2 load")
ggsave(paste("exp2_lung_integrated_log10SCoV2_load", "pdf", sep="."), useDingbats=FALSE)

unfilt_metadata <- read.table("./exp2_lung_combined_250_metadata.txt")
exp2_lung@meta.data$cell_id <- rownames(exp2_lung@meta.data)
df <- dplyr::left_join(exp2_lung@meta.data, unfilt_metadata, by="cell_id")
exp2_lung@meta.data$DoubSing <- df$DoubSing
exp2_lung@meta.data$endocytotox <- df$endocytotox

#cell type annotation
#LOC101833790 = Ccr5

MELC <- c("")
markers <- c("Ace2", "Ackr1", "Acta2", "Adgre1", "Arg1", "Aspn", "C1qb", "Camp", "Ccl21", "Ccr2", "LOC101833790", "Cd3e", "Cd4", "Cd79b", "Cd8a", "Cldn5", "Col1a2", "Cx3cr1", "Cxcr2", "Dcn", "Ednrb", "Fcgr4", "Flt3", "Foxj1", "Gzma", "Il7r", "Irf8", "Lamp3", "Marco", "Mrc1", "Ms4a1", "Nkg7", "Pdpn", "Ppbp", "Plvap", "Retn", "Rtkn2", "S100a8", "Siglecf", "Tagln", "Tcf4", "Treml4")
for (gene in markers) {
  df <- FetchData(exp2_lung, gene)
  gene <- str_replace_all(gene, "-", "_")
  colnames(df) <- gene
  df = cbind(df, exp2_lung@reductions$umap@cell.embeddings)
  ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
    geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
    scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    coord_fixed(ratio=1)+
    ggtitle(gene)
  ggsave(paste("exp2_lung_integrated", gene, "pdf", sep="."), useDingbats=FALSE)
}

avg <- AverageExpression(exp2_lung, assays=c("SCT"), features = markers, return.seurat = T, slot="data") 
avg2 <- AverageExpression(exp2_lung, assays=c("SCT"), features = markers, return.seurat = F, slot="data")
avg2 <- as.matrix(avg2$SCT)
hc <- hclust(dist(avg2))
gene.order <- row.names(avg2)[order.hclust(hc)]
DoHeatmap(avg, size=5, features=gene.order)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("exp2_lung_integrated_cluster_heatmap.pdf", width=12, height=7)

#Myofibroblasts, fibroblasts, smooth muscle
#LOC101844074 is Cox4i2 (smooth muscle marker)
musclemarkers <- c("Aspn", "Dcn", "Acta2", "Tagln", "Wif1", "Fgf18", "Col1a2", "Bsg", "LOC101844074", "Cnn1", "Myh11", "Actg2", "Gpc3", "Apoe", "Serpinf1", "Gpc3")
avg <- AverageExpression(exp2_lung, assays=c("SCT"), features = musclemarkers, return.seurat = T, slot="data") 
avg2 <- AverageExpression(exp2_lung, assays=c("SCT"), features = musclemarkers, return.seurat = F, slot="data")
avg2 <- as.matrix(avg2$SCT)
hc <- hclust(dist(avg2))
gene.order <- row.names(avg2)[order.hclust(hc)]
DoHeatmap(avg, size=5, features=gene.order)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("exp2_lung_integrated_musclemarkers_heatmap.pdf")

Idents(exp2_lung) <- exp2_lung@meta.data$seurat_clusters

exp2_lung <- RenameIdents(exp2_lung, 
                       '0'='Bcells',
                       '1'='AlveolarMacrophages',
                       '2'='Endothelial',
                       '3'='AT2',
                       '4'='Neutrophils',
                       '5'='TNKcells',
                       '6'='MonocyticMacrophages',
                       '7'='AT1',
                       '8'='InterstitialMacrophages',
                       '9'='MonocyticMacrophages',
                       '10'='InterstitialMacrophages',
                       '11'='Treml4+Macrophages',
                       '12'='Endothelial',
                       '13'='Endothelial',
                       '14'='TNKcells',
                       '15'='MyeloidDendritic',
                       '16'='TNKcells',
                       '17'='Treml4+Macrophages',
                       '18'='Myofibroblast',
                       '19'='mixed1',
                       '20'='Platelets',
                       '21'='mixed2',
                       '22'='TNKcells',
                       '23'='Endothelial',
                       '24'='Endothelial',
                       '25'='SmoothMuscle',
                       '26'='AlveolarMacrophages',
                       '27'='Neutrophils',
                       '28'='mixed3',
                       '29'='Fibroblasts',
                       '30'='pDC',
                       '31'='Bcells',
                       '32'='InterstitialMacrophages',
                       '33'='Ciliated',
                       '34'='mixed4',
                       '35'='mixed5')
exp2_lung@meta.data$celltype <- Idents(exp2_lung)

saveRDS(exp2_lung, "./exp2_lung_combined_integrated_annotated.rds")


#UMAP with cell types
means <- exp2_lung@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(mean_U1 = mean(UMAP_1), mean_U2 = mean(UMAP_2)) 

ggplot()+
  geom_point(data=exp2_lung@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=celltype), size=0.5, shape=16)+
  geom_text(data=means, aes(x=mean_U1, y=mean_U2, label=celltype), size=3)+
  coord_fixed(ratio=1)+
  ggtitle("celltype")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  scale_colour_manual(values = celltypecolors)
ggsave("exp2_lung_integrated_celltypes.pdf", useDingbats=FALSE, width=14, height=7)

#percent cell type per treatment
df1 <- cbind.data.frame(exp2_lung@meta.data$SCoV2_load,
                        exp2_lung@meta.data$celltype,
                        exp2_lung@meta.data$treatment,
                        exp2_lung@meta.data$hamster)
colnames(df1) <- c("SCoV2_load", "celltype", "treatment", "hamster")
a = df1 %>% group_by(treatment, hamster) %>% tally(name="tot") 
b = df1 %>% group_by(celltype, treatment, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('treatment', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
write.table(c, "celltypefractions.tsv", sep="\t", quote=FALSE, row.names = FALSE)
c <- subset(c, celltype %in% c("MonocyticMacrophages", "Neutrophils"))

tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("celltype", "treatment"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")))+
  geom_bar(aes(colour=treatment), position=position_dodge(0.9), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.9), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("celltypes per timepoint", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_colors)+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))

ggsave("exp2_lung_integrated_celltypepercentage_subs_pertimepoint.pdf", useDingbats=FALSE)

#percent virus positive by treatment
a = df1 %>% group_by(celltype, treatment, hamster) %>% tally(name="tot") 
b = df1 %>% filter(SCoV2_load>0) %>% group_by(celltype, treatment, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('celltype', 'treatment', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
#write.table(c, paste("days", gene, "tsv", sep="."), sep="\t", quote=FALSE, row.names = FALSE)
tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("celltype", "treatment"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")))+
  geom_bar(aes(colour=treatment), position=position_dodge(0.9), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.9), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("virus positive per treatment", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_colors)+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))
ggsave("barplot_viruspositivepercelltype_pertimepoint.pdf", useDingbats=FALSE)

#doublets by celltype and treatment
df1 <- cbind.data.frame(exp2_lung@meta.data$DoubSing,
                        exp2_lung@meta.data$endocytotox,
                        exp2_lung@meta.data$celltype,
                        exp2_lung@meta.data$treatment,
                        exp2_lung@meta.data$hamster)
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
ggsave("exp2_lung_combined_250_DoubSing.pdf", useDingbats=FALSE, width=14, height=7)

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
ggsave("exp2_lung_combined_250_endocytotox.pdf", useDingbats=FALSE, width=14, height=7)



###################
#Density plots
#smooth dimplots
cond1="aaUntr"
cond2="adeno2x"
smooth_DimPlot(subset(exp2_lung, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("exp2_lung_density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)
cond1="aaUntr"
cond2="att2x"
smooth_DimPlot(subset(exp2_lung, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("exp2_lung_density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)
cond1="aaUntr"
cond2="mRNA2x"
smooth_DimPlot(subset(exp2_lung, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("exp2_lung_density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)
cond1="aaUntr"
cond2="mRNAatt"
smooth_DimPlot(subset(exp2_lung, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("exp2_lung_density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)




