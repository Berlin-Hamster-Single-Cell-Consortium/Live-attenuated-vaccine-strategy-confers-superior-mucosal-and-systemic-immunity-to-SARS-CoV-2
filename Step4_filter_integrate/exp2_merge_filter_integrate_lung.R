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





hamster_all <- merge(aaUntr_1, y = c(aaUntr_2, aaUntr_3, aaUntr_4, att2x_1, att2x_2, att2x_3, adeno2x_1, adeno2x_2, adeno2x_3, mRNA2x_1, mRNA2x_2, mRNA2x_3, mRNAatt_1, mRNAatt_2, mRNAatt_3), add.cell.ids = c("aaUntr_1", "aaUntr_2", "aaUntr_3", "aaUntr_4", "att2x_1", "att2x_2", "att2x_3", "adeno2x_1", "adeno2x_2", "adeno2x_3", "mRNA2x_1", "mRNA2x_2", "mRNA2x_3", "mRNAatt_1", "mRNAatt_2", "mRNAatt_3"), project = "impf_exp2")

hamster_all@meta.data %>% group_by(orig.ident) %>% tally()



cells_to_keep <- read.table("./exp2_lung_cells_to_keep.txt") 
hamster_all <- subset(hamster_all, cells=cells_to_keep$cell_id)

hamster_all@meta.data %>% group_by(orig.ident) %>% tally()


DefaultAssay(hamster_all) <- "RNA"


## integrate by hamster to remove batch effects

hamster.list <- SplitObject(hamster_all, split.by = "orig.ident")
for (i in 1:length(hamster.list)){
  hamster.list[[i]] <- SCTransform(hamster.list[[i]], verbose =T)
}

hamster.features <- SelectIntegrationFeatures(object.list = hamster.list, nfeatures = 3000)
hamster.list <- PrepSCTIntegration(object.list = hamster.list, anchor.features = hamster.features, 
                                   verbose = T)

hamster.anchors <- FindIntegrationAnchors(object.list = hamster.list, normalization.method = "SCT", 
                                          anchor.features = hamster.features, verbose = T)
hamster.integrated <- IntegrateData(anchorset = hamster.anchors, normalization.method = "SCT", 
                                    verbose = T)

## run dimensional reductions
#   PCA
hamster.integrated<- RunPCA(hamster.integrated, verbose = FALSE)
#   UMAP
hamster.integrated<- RunUMAP(hamster.integrated, dims = 1:30, verbose = FALSE)

hamster.integrated <- FindNeighbors(hamster.integrated, dims = 1:30)
hamster.integrated <- FindClusters(hamster.integrated, resolution = 0.5)


saveRDS(hamster.integrated, "./seu_exp2_lung_new_combined_integrated.rds")
