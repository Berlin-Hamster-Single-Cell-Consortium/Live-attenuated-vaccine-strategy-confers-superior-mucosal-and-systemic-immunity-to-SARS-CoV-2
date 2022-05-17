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

#Fig. 2, S5 
#pseudobulk_exp2_lung

exp2_lung <- readRDS("./exp2_lung_combined_integrated_annotated.rds")



expr <- list()
for (cluster in unique(Idents(exp2_lung))) {
  for (sample in unique(exp2_lung@meta.data$orig.ident)) {
    cells <- Cells(exp2_lung)[(exp2_lung@meta.data$orig.ident==sample) & (Idents(exp2_lung)==cluster)]
    if (length(cells) > 5) { 
      expr[[paste0(cluster,'_',sample)]] <- rowSums(exp2_lung@assays$RNA@counts[,cells])
    }
  }
}
for (sample in unique(exp2_lung@meta.data$orig.ident)) {
  cells <- Cells(exp2_lung)[(exp2_lung@meta.data$orig.ident==sample)]
  expr[[paste0('all_',sample)]] <- rowSums(exp2_lung@assays$RNA@counts[,cells])
}

colData <- data.frame(cell.type=factor(gsub('_.*','',names(expr))),
                      donor=gsub('[^_]*_([^_]*)_([0-9])','hamster_\\2',names(expr)),
                      treatment=gsub('[^_]*_([^_]*)_([0-9])','\\1',names(expr)),
                      group=names(expr),
                      row.names=names(expr))

counts <- do.call(cbind,expr)

clusters_to_check=c('TNKcells', 'Bcells', 'AlveolarMacrophages', 'Endothelial', 'Myofibroblasts', 'MonocyticMacrophages', 'AT2', 'InterstitialMacrophages', 'Treml4+Macrophages', 'Neutrophils', 'SmoothMuscle', 'MyeloidDendritic', 'Platelets', 'AT1')

res <- list()
for (cluster in c('all',clusters_to_check)) {
  take_row <- rowSums(counts) > 0
  take_col <- (colSums(counts) > 0) & (colData[,'cell.type']==cluster)
  try({
    dds <- DESeqDataSetFromMatrix(countData=counts[take_row,take_col],
                                  colData=colData[take_col,,drop=FALSE],
                                  design=~treatment)
    if (cluster!='all')
      dds <- estimateSizeFactors(dds, type='poscounts')
    dds <- DESeq(dds)
    res[[paste0(cluster,'_adeno2xvsuntr')]] <- lfcShrink(dds,
                                                         contrast=c('treatment','adeno2x','aaUntr'),
                                                         type='normal',
                                                         format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene_name') %>%
      dplyr::mutate(cluster=cluster,contrast='adeno2xvsuntr')
    res[[paste0(cluster,'_att2xvsuntr')]] <- lfcShrink(dds,
                                                       contrast=c('treatment','att2x','aaUntr'),
                                                       type='normal',
                                                       format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene_name') %>%
      dplyr::mutate(cluster=cluster,contrast='att2xvsuntr')
    res[[paste0(cluster,'_mRNA2xvsuntr')]] <- lfcShrink(dds,
                                                        contrast=c('treatment','mRNA2x','aaUntr'),
                                                        type='normal',
                                                        format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene_name') %>%
      dplyr::mutate(cluster=cluster,contrast='mRNA2xvsuntr')
    res[[paste0(cluster,'_mRNAattvsuntr')]] <- lfcShrink(dds,
                                                         contrast=c('treatment','mRNAatt','aaUntr'),
                                                         type='normal',
                                                         format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene_name') %>%
      dplyr::mutate(cluster=cluster,contrast='mRNAattvsuntr')
  })
}

pseudobulk_exp2_lung <- do.call(rbind,res)

write.table(pseudobulk_exp2_lung, "./pseudobulk_exp2_lung.txt", row.names = TRUE, sep="\t", quote=FALSE)


clusters=c('AlveolarMacrophages', 'InterstitialMacrophages', 'MonocyticMacrophages', 'Treml4+Macrophages', 'Neutrophils', 'TNKcells', 'Bcells', 'AT1', 'AT2', 'Endothelial')

contrasts_to_use <- c("adeno2xvsuntr", "att2xvsuntr", "mRNA2xvsuntr", "mRNAattvsuntr")

#Option 1: select specific genes
genes <- c("Irf7", "Mx2", "Isg15", "Ifit3", "Lgals3bp", "Cxcl10", "Ccl5", "Tnfsf10")
the_filename="pseudobulk_exp2_lung_inflISG_all_facet"


#Option 2: To select top x genes per contrast, do:
genes <- pseudobulk_exp2_lung %>%
  dplyr::filter(cluster %in% clusters) %>%
  dplyr::filter(contrast %in% contrasts_to_use) %>%
  dplyr::filter(!grepl('SCoV2',gene_name)) %>%
  dplyr::filter(!(is.na(padj)) & (padj < .01)) %>%
  dplyr::group_by(contrast, cluster) %>%
  dplyr::slice_min(order_by=padj,n=5,with_ties=FALSE) %>%
  dplyr::pull(gene_name)
the_filename="pseudobulk_exp2_lung_top5_all_facet"



#Make plot
df <- pseudobulk_exp2_lung %>%
  dplyr::filter(cluster %in% clusters) %>%
  dplyr::filter(contrast %in% contrasts_to_use) %>%
  #dplyr::filter(!(is.na(padj)) & (padj < .01) & (abs(log2FoldChange) >= 0.7)) %>%
  dplyr::filter(gene_name %in% genes) %>%
  dplyr::mutate(group=paste0(cluster,'_',contrast))  %>%
  dplyr::select(group,gene_name,log2FoldChange) %>%
  spread(group,log2FoldChange) %>%
  tibble::column_to_rownames('gene_name')

df[is.na(df)] <- 0
hc <- hclust(dist(df))
gene.order <- row.names(df)[order.hclust(hc)]
hc <- hclust(dist(t(df)))
group.order <- colnames(df)[order.hclust(hc)]
group.order = sub("(.+)_.*", "\\1", group.order)

genes <- row.names(df)

ggplot(pseudobulk_exp2_lung %>%
         dplyr::filter(cluster %in% clusters) %>%
         dplyr::filter(contrast %in% contrasts_to_use) %>%
         dplyr::filter(gene_name %in% genes) %>%
         dplyr::mutate(padj=ifelse(padj < 1E-20, 1E-20, padj)) %>%
         dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
         dplyr::mutate(gene=factor(gene_name,levels=gene.order)) %>%
         dplyr::mutate(cluster=factor(cluster,levels=clusters)),
       aes(y=gene,x=cluster,size=-log10(padj),color=log2FoldChange)) +
  geom_point(shape=16) +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(2,4,6),
                        labels=c("1e-2","1e-4", "<1e-6"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-3,3),oob=scales::squish) +
  facet_wrap(~contrast, nrow=1) +
  coord_cartesian(clip = 'off') +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
  labs(x='',y='')

#save 8x12 inches
ggsave(paste("exp2_lung", the_filename, "pdf", sep="."), width=10, height = length(genes)/6+2, units="in", useDingbats=FALSE)
