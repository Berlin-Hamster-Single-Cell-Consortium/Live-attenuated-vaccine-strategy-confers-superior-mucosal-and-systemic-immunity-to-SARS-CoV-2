library(ggplot2)
library(stringr)
library(Seurat)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
require(gridExtra)
library(dplyr)
library(DoubletFinder)

data <- Read10X(data.dir = c("exp2/Pool2/part1/part1/outs/per_sample_outs/LU_X1/count/sample_feature_bc_matrix/", "exp2/Pool2/part2/part2/outs/per_sample_outs/LU_X1/count/sample_feature_bc_matrix/", "exp2/Pool2/part3/part3/outs/per_sample_outs/LU_X1/count/sample_feature_bc_matrix/"), gene.column = 2)
aaUntr_1 = CreateSeuratObject(counts = data$`Gene Expression`, project="aaUntr_1", min.cells=5, min.features=250)
aaUntr_1@meta.data$orig.ident <- "aaUntr_1"

data <- Read10X(data.dir = c("exp2/Pool2/part1/part1/outs/per_sample_outs/LU_X2/count/sample_feature_bc_matrix/", "exp2/Pool2/part2/part2/outs/per_sample_outs/LU_X2/count/sample_feature_bc_matrix/", "exp2/Pool2/part3/part3/outs/per_sample_outs/LU_X2/count/sample_feature_bc_matrix/"), gene.column = 2)
aaUntr_2 = CreateSeuratObject(counts = data$`Gene Expression`, project="aaUntr_2", min.cells=5, min.features=250)
aaUntr_2@meta.data$orig.ident <- "aaUntr_2"

data <- Read10X(data.dir = c("exp2/Pool5/part1/part1/outs/per_sample_outs/LU_X3/count/sample_feature_bc_matrix/", "exp2/Pool5/part2/part2/outs/per_sample_outs/LU_X3/count/sample_feature_bc_matrix/", "exp2/Pool5/part3/part3/outs/per_sample_outs/LU_X3/count/sample_feature_bc_matrix/"), gene.column = 2)
aaUntr_3 = CreateSeuratObject(counts = data$`Gene Expression`, project="aaUntr_3", min.cells=5, min.features=250)
aaUntr_3@meta.data$orig.ident <- "aaUntr_3"

data <- Read10X(data.dir = c("exp2/Pool5/part1/part1/outs/per_sample_outs/LU_X4/count/sample_feature_bc_matrix/", "exp2/Pool5/part2/part2/outs/per_sample_outs/LU_X4/count/sample_feature_bc_matrix/", "exp2/Pool5/part3/part3/outs/per_sample_outs/LU_X4/count/sample_feature_bc_matrix/"), gene.column = 2)
aaUntr_4 = CreateSeuratObject(counts = data$`Gene Expression`, project="aaUntr_4", min.cells=5, min.features=250)
aaUntr_4@meta.data$orig.ident <- "aaUntr_4"

data <- Read10X(data.dir = c("exp2/Pool2/part1/part1/outs/per_sample_outs/LU_A1/count/sample_feature_bc_matrix/", "exp2/Pool2/part2/part2/outs/per_sample_outs/LU_A1/count/sample_feature_bc_matrix/", "exp2/Pool2/part3/part3/outs/per_sample_outs/LU_A1/count/sample_feature_bc_matrix/"), gene.column = 2)
att2x_1 = CreateSeuratObject(counts = data$`Gene Expression`, project="att2x_1", min.cells=5, min.features=250)
att2x_1@meta.data$orig.ident <- "att2x_1"

data <- Read10X(data.dir = c("exp2/Pool2/part1/part1/outs/per_sample_outs/LU_A2/count/sample_feature_bc_matrix/", "exp2/Pool2/part2/part2/outs/per_sample_outs/LU_A2/count/sample_feature_bc_matrix/", "exp2/Pool2/part3/part3/outs/per_sample_outs/LU_A2/count/sample_feature_bc_matrix/"), gene.column = 2)
att2x_2 = CreateSeuratObject(counts = data$`Gene Expression`, project="att2x_2", min.cells=5, min.features=250)
att2x_2@meta.data$orig.ident <- "att2x_2"

data <- Read10X(data.dir = c("exp2/Pool2/part1/part1/outs/per_sample_outs/LU_A3/count/sample_feature_bc_matrix/", "exp2/Pool2/part2/part2/outs/per_sample_outs/LU_A3/count/sample_feature_bc_matrix/", "exp2/Pool2/part3/part3/outs/per_sample_outs/LU_A3/count/sample_feature_bc_matrix/"), gene.column = 2)
att2x_3 = CreateSeuratObject(counts = data$`Gene Expression`, project="att2x_3", min.cells=5, min.features=250)
att2x_3@meta.data$orig.ident <- "att2x_3"

data <- Read10X(data.dir = c("exp2/Pool2/part1/part1/outs/per_sample_outs/LU_B1/count/sample_feature_bc_matrix/", "exp2/Pool2/part2/part2/outs/per_sample_outs/LU_B1/count/sample_feature_bc_matrix/", "exp2/Pool2/part3/part3/outs/per_sample_outs/LU_B1/count/sample_feature_bc_matrix/"), gene.column = 2)
adeno2x_1 = CreateSeuratObject(counts = data$`Gene Expression`, project="adeno2x_1", min.cells=5, min.features=250)
adeno2x_1@meta.data$orig.ident <- "adeno2x_1"

data <- Read10X(data.dir = c("exp2/Pool2/part1/part1/outs/per_sample_outs/LU_B2/count/sample_feature_bc_matrix/", "exp2/Pool2/part2/part2/outs/per_sample_outs/LU_B2/count/sample_feature_bc_matrix/", "exp2/Pool2/part3/part3/outs/per_sample_outs/LU_B2/count/sample_feature_bc_matrix/"), gene.column = 2)
adeno2x_2 = CreateSeuratObject(counts = data$`Gene Expression`, project="adeno2x_2", min.cells=5, min.features=250)
adeno2x_2@meta.data$orig.ident <- "adeno2x_2"

data <- Read10X(data.dir = c("exp2/Pool2/part1/part1/outs/per_sample_outs/LU_B3/count/sample_feature_bc_matrix/", "exp2/Pool2/part2/part2/outs/per_sample_outs/LU_B3/count/sample_feature_bc_matrix/", "exp2/Pool2/part3/part3/outs/per_sample_outs/LU_B3/count/sample_feature_bc_matrix/"), gene.column = 2)
adeno2x_3 = CreateSeuratObject(counts = data$`Gene Expression`, project="adeno2x_3", min.cells=5, min.features=250)
adeno2x_3@meta.data$orig.ident <- "adeno2x_3"

data <- Read10X(data.dir = c("exp2/Pool5/part1/part1/outs/per_sample_outs/LU_C1/count/sample_feature_bc_matrix/", "exp2/Pool5/part2/part2/outs/per_sample_outs/LU_C1/count/sample_feature_bc_matrix/", "exp2/Pool5/part3/part3/outs/per_sample_outs/LU_C1/count/sample_feature_bc_matrix/"), gene.column = 2)
mRNA2x_1 = CreateSeuratObject(counts = data$`Gene Expression`, project="mRNA2x_1", min.cells=5, min.features=250)
mRNA2x_1@meta.data$orig.ident <- "mRNA2x_1"

data <- Read10X(data.dir = c("exp2/Pool5/part1/part1/outs/per_sample_outs/LU_C2/count/sample_feature_bc_matrix/", "exp2/Pool5/part2/part2/outs/per_sample_outs/LU_C2/count/sample_feature_bc_matrix/", "exp2/Pool5/part3/part3/outs/per_sample_outs/LU_C2/count/sample_feature_bc_matrix/"), gene.column = 2)
mRNA2x_2 = CreateSeuratObject(counts = data$`Gene Expression`, project="mRNA2x_2", min.cells=5, min.features=250)
mRNA2x_2@meta.data$orig.ident <- "mRNA2x_2"

data <- Read10X(data.dir = c("exp2/Pool5/part1/part1/outs/per_sample_outs/LU_C3/count/sample_feature_bc_matrix/", "exp2/Pool5/part2/part2/outs/per_sample_outs/LU_C3/count/sample_feature_bc_matrix/", "exp2/Pool5/part3/part3/outs/per_sample_outs/LU_C3/count/sample_feature_bc_matrix/"), gene.column = 2)
mRNA2x_3 = CreateSeuratObject(counts = data$`Gene Expression`, project="mRNA2x_3", min.cells=5, min.features=250)
mRNA2x_3@meta.data$orig.ident <- "mRNA2x_3"

data <- Read10X(data.dir = c("exp2/Pool5/part1/part1/outs/per_sample_outs/LU_D1/count/sample_feature_bc_matrix/", "exp2/Pool5/part2/part2/outs/per_sample_outs/LU_D1/count/sample_feature_bc_matrix/", "exp2/Pool5/part3/part3/outs/per_sample_outs/LU_D1/count/sample_feature_bc_matrix/"), gene.column = 2)
mRNAatt_1 = CreateSeuratObject(counts = data$`Gene Expression`, project="mRNAatt_1", min.cells=5, min.features=250)
mRNAatt_1@meta.data$orig.ident <- "mRNAatt_1"

data <- Read10X(data.dir = c("exp2/Pool5/part1/part1/outs/per_sample_outs/LU_D2/count/sample_feature_bc_matrix/", "exp2/Pool5/part2/part2/outs/per_sample_outs/LU_D2/count/sample_feature_bc_matrix/", "exp2/Pool5/part3/part3/outs/per_sample_outs/LU_D2/count/sample_feature_bc_matrix/"), gene.column = 2)
mRNAatt_2 = CreateSeuratObject(counts = data$`Gene Expression`, project="mRNAatt_2", min.cells=5, min.features=250)
mRNAatt_2@meta.data$orig.ident <- "mRNAatt_2"

data <- Read10X(data.dir = c("exp2/Pool5/part1/part1/outs/per_sample_outs/LU_D3/count/sample_feature_bc_matrix/", "exp2/Pool5/part2/part2/outs/per_sample_outs/LU_D3/count/sample_feature_bc_matrix/", "exp2/Pool5/part3/part3/outs/per_sample_outs/LU_D3/count/sample_feature_bc_matrix/"), gene.column = 2)
mRNAatt_3 = CreateSeuratObject(counts = data$`Gene Expression`, project="mRNAatt_3", min.cells=5, min.features=250)
mRNAatt_3@meta.data$orig.ident <- "mRNAatt_3"


obj_list <- list(aaUntr_1, aaUntr_2, aaUntr_3, aaUntr_4, att2x_1, att2x_2, att2x_3, adeno2x_1, adeno2x_2, adeno2x_3, mRNA2x_1, mRNA2x_2, mRNA2x_3, mRNAatt_1, mRNAatt_2, mRNAatt_3)
names(obj_list) <- c("aaUntr_1", "aaUntr_2", "aaUntr_3", "aaUntr_4", "att2x_1", "att2x_2", "att2x_3", "adeno2x_1", "adeno2x_2", "adeno2x_3", "mRNA2x_1", "mRNA2x_2", "mRNA2x_3", "mRNAatt_1", "mRNAatt_2", "mRNAatt_3") 

expr <- list()

for (the_file in names(obj_list)) {
  seu <- obj_list[[the_file]]
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  seu <- RunUMAP(seu, dims = 1:18)
  seu <- FindNeighbors(seu, dims = 1:18)
  seu <- FindClusters(seu, resolution = 0.9)
  
  seu@meta.data = cbind(seu@meta.data, seu@reductions$umap@cell.embeddings)
  sweep.res.list_blood <- paramSweep_v3(seu, PCs = 1:18, sct = FALSE)
  sweep.stats_blood <- summarizeSweep(sweep.res.list_blood, GT = FALSE)
  bcmvn_blood <- find.pK(sweep.stats_blood)
  homotypic.prop <- modelHomotypic(seu@meta.data$seurat_clusters)
  nExp_poi <- round(0.05*nrow(seu@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(seu@meta.data) <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "RNA_snn_res.0.9", "seurat_clusters", "UMAP_1", "UMAP_2", "pAnn", "DoubSing")
  ggplot()+geom_point(data=seu@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=DoubSing), size=0.25)+theme_bw()
  ggsave(paste0(the_file, "_DoubSing.pdf"))
  expr[[the_file]] <- seu@meta.data
}
DoubSing <- do.call(rbind,expr)
DoubSing$cell_id <- rownames(DoubSing)
DoubSing <- DoubSing %>% select(c("cell_id", "DoubSing")) %>% mutate (cell_id = gsub("\\.", "\\_", cell_id))

write.table(DoubSing, "./DoubSing_lung_exp2.txt", sep="\t", quote=FALSE)



hamster_all <- merge(aaUntr_1, y = c(aaUntr_2, aaUntr_3, aaUntr_4, att2x_1, att2x_2, att2x_3, adeno2x_1, adeno2x_2, adeno2x_3, mRNA2x_1, mRNA2x_2, mRNA2x_3, mRNAatt_1, mRNAatt_2, mRNAatt_3), add.cell.ids = c("aaUntr_1", "aaUntr_2", "aaUntr_3", "aaUntr_4", "att2x_1", "att2x_2", "att2x_3", "adeno2x_1", "adeno2x_2", "adeno2x_3", "mRNA2x_1", "mRNA2x_2", "mRNA2x_3", "mRNAatt_1", "mRNAatt_2", "mRNAatt_3"), project = "impf_exp2")

hamster_all@meta.data %>% group_by(orig.ident) %>% tally()


saveRDS(hamster_all, "./seu_lung_exp2_combined_250.rds")

