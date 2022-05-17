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




#read integrated object
exp1_lung <- readRDS("./seu_exp1_lung_new_combined_integrated.rds")
exp1_lung@meta.data = cbind(exp1_lung@meta.data, exp1_lung@reductions$umap@cell.embeddings)
exp1_lung@meta.data$treatment <- gsub("_[1-4].*","",exp1_lung@meta.data$orig.ident)
exp1_lung@meta.data$hamster <- gsub(".*_([1-4].*)","Ha\\1",exp1_lung@meta.data$orig.ident)

DefaultAssay(exp1_lung) <- 'RNA'
SCoV2_rawcounts <- FetchData(exp1_lung, grep("SCoV2", exp1_lung@assays$RNA@counts@Dimnames[[1]], value="TRUE"), slot="counts")
exp1_lung@meta.data$SCoV2_sum <- rowSums(SCoV2_rawcounts)
SCoV2_rawcounts$SCoV2sum <- rowSums(SCoV2_rawcounts)
exp1_lung@meta.data$SCoV2_load <- SCoV2_rawcounts$SCoV2sum/exp1_lung@meta.data$nCount_RNA*100
DefaultAssay(exp1_lung) <- 'SCT'

ggplot()+
  geom_point(data=exp1_lung@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=treatment), size=0.5, shape=16, alpha=0.5)+
  ggtitle("treatment")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("exp1_lung_integrated_treatment.pdf", useDingbats=FALSE)
ggplot()+
  geom_point(data=exp1_lung@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=hamster), size=0.2, shape=16, alpha=0.5)+
  ggtitle("hamster")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("exp1_lung_integrated_hamsters.pdf", useDingbats=FALSE)
UMAPPlot(exp1_lung, label=TRUE)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("Clusters")
ggsave("exp1_lung_integrated_clusters.pdf", useDingbats=FALSE)



ggplot()+geom_point(data=subset(exp1_lung@meta.data, SCoV2_load<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", size=0.5, shape = 16)+
  geom_point(data=subset(exp1_lung@meta.data, SCoV2_load>0), aes(x=UMAP_1, y=UMAP_2, colour=log10(SCoV2_load)), size=0.5, shape = 16)+
  scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("log10 SCoV2 load")
ggsave(paste("exp1_lung_integrated_log10SCoV2_load", "pdf", sep="."), useDingbats=FALSE)

unfilt_metadata <- read.table("./exp1_lung_combined_250_metadata.txt")
exp1_lung@meta.data$cell_id <- rownames(exp1_lung@meta.data)
df <- dplyr::left_join(exp1_lung@meta.data, unfilt_metadata, by="cell_id")
exp1_lung@meta.data$DoubSing <- df$DoubSing
exp1_lung@meta.data$endocytotox <- df$endocytotox

#cell type annotation
#LOC101833790 = Ccr5

markers <- c("Ace2", "Ackr1", "Acta2", "Adgre1", "Arg1", "Aspn", "C1qb", "Camp", "Ccl21", "Ccr2", "LOC101833790", "Cd3e", "Cd4", "Cd79b", "Cd8a", "Cldn5", "Col1a2", "Cx3cr1", "Cxcr2", "Dcn", "Ednrb", "Fcgr4", "Flt3", "Foxj1", "Gzma", "Il7r", "Irf8", "Lamp3", "Marco", "Mrc1", "Ms4a1", "Nkg7", "Pdpn", "Ppbp", "Plvap", "Retn", "Rtkn2", "S100a8", "Siglecf", "Tagln", "Tcf4", "Treml4")
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
  ggsave(paste("exp1_lung_integrated", gene, "pdf", sep="."), useDingbats=FALSE)
}

avg <- AverageExpression(exp1_lung, assays=c("SCT"), features = markers, return.seurat = T, slot="data") 
avg2 <- AverageExpression(exp1_lung, assays=c("SCT"), features = markers, return.seurat = F, slot="data")
avg2 <- as.matrix(avg2$SCT)
hc <- hclust(dist(avg2))
gene.order <- row.names(avg2)[order.hclust(hc)]
DoHeatmap(avg, size=5, features=gene.order)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("exp1_lung_integrated_cluster_heatmap.pdf", width=12, height=7)

#Myofibroblasts, fibroblasts, smooth muscle
#LOC101844074 is Cox4i2 (smooth muscle marker)
musclemarkers <- c("Aspn", "Dcn", "Acta2", "Tagln", "Wif1", "Fgf18", "Col1a2", "Bsg", "LOC101844074", "Cnn1", "Myh11", "Actg2", "Gpc3", "Apoe", "Serpinf1", "Gpc3")
avg <- AverageExpression(exp1_lung, assays=c("SCT"), features = musclemarkers, return.seurat = T, slot="data") 
avg2 <- AverageExpression(exp1_lung, assays=c("SCT"), features = musclemarkers, return.seurat = F, slot="data")
avg2 <- as.matrix(avg2$SCT)
hc <- hclust(dist(avg2))
gene.order <- row.names(avg2)[order.hclust(hc)]
DoHeatmap(avg, size=5, features=gene.order)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("exp1_lung_integrated_musclemarkers_heatmap.pdf")

Idents(exp1_lung) <- exp1_lung@meta.data$seurat_clusters

exp1_lung <- RenameIdents(exp1_lung, 
                          '0'='AlveolarMacrophages',
                          '1'='Neutrophils',
                          '2'='MonocyticMacrophages',
                          '3'='InterstitialMacrophages',
                          '4'='Endothelial',
                          '5'='AT2',
                          '6'='Neutrophils',
                          '7'='Bcells',
                          '8'='mixed1',
                          '9'='Treml4+Macrophages',
                          '10'='Treml4+Macrophages',
                          '11'='AT1',
                          '12'='mixed2',
                          '13'='InterstitialMacrophages',
                          '14'='TNKcells',
                          '15'='TNKcells',
                          '16'='MonocyticMacrophages',
                          '17'='Platelets',
                          '18'='Endothelial',
                          '19'='Bcells',
                          '20'='Treml4+Macrophages',
                          '21'='Neutrophils',
                          '22'='TNKcells',
                          '23'='InterstitialMacrophages',
                          '24'='AT2',
                          '25'='mixed3',
                          '26'='Myofibroblast',
                          '27'='MyeloidDendritic',
                          '28'='Endothelial',
                          '29'='pDC',
                          '30'='SmoothMuscle',
                          '31'='Fibroblasts',
                          '32'='TNKcells',
                          '33'='mixed4',
                          '34'='Bcells',
                          '35'='mixed5',
                          '36'='AlveolarMacrophages',
                          '37'='Ciliated',
                          '38'='mixed6')
exp1_lung@meta.data$celltype <- Idents(exp1_lung)

saveRDS(exp1_lung, "./exp1_lung_combined_integrated_annotated.rds")


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
ggsave("exp1_lung_integrated_celltypes.pdf", useDingbats=FALSE, width=14, height=7)

#percent cell type per treatment
df1 <- cbind.data.frame(exp1_lung@meta.data$SCoV2_load,
                        exp1_lung@meta.data$celltype,
                        exp1_lung@meta.data$treatment,
                        exp1_lung@meta.data$hamster)
colnames(df1) <- c("SCoV2_load", "celltype", "treatment", "hamster")

a = df1 %>% group_by(treatment, hamster) %>% tally(name="tot") 
b = df1 %>% group_by(celltype, treatment, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('treatment', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
write.table(c, "exp1_lung_celltypefractions.tsv", sep="\t", quote=FALSE, row.names = FALSE)

tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("celltype", "treatment"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")))+
  geom_bar(aes(colour=treatment), position=position_dodge(.8), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("celltypes per timepoint", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_colors)+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))

ggsave("exp1_lung_integrated_celltypepercentage_pertimepoint.pdf", useDingbats=FALSE)

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
  geom_bar(aes(colour=treatment), position=position_dodge(), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("virus positive per treatment", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_colors)+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))
ggsave("barplot_viruspositivepercelltype_pertimepoint.pdf", useDingbats=FALSE)

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


###################
#Density plots
#smooth dimplots
cond1="aaUntr"
cond2="adeno"
smooth_DimPlot(subset(exp1_lung, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("exp1_lung_density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)
cond1="aaUntr"
cond2="attenuated"
smooth_DimPlot(subset(exp1_lung, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("exp1_lung_density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)
cond1="aaUntr"
cond2="mRNA"
smooth_DimPlot(subset(exp1_lung, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("exp1_lung_density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)


