library(ggplot2)
library(stringr)
library(Seurat)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
require(gridExtra)
library(dplyr)
library(DoubletFinder)


data <- Read10X(data.dir = c("Pool3/part1/part1/outs/per_sample_outs/NM_D1/count/sample_feature_bc_matrix/", "Pool3/part2/part2/outs/per_sample_outs/NM_D1/count/sample_feature_bc_matrix/", "Pool3/part3/part3/outs/per_sample_outs/NM_D1/count/sample_feature_bc_matrix/", "Pool3/part4/part4/outs/per_sample_outs/NM_D1/count/sample_feature_bc_matrix/"), gene.column = 2)
aaUntr_1 = CreateSeuratObject(counts = data$`Gene Expression`, project="aaUntr_1", min.cells=5, min.features=250)
aaUntr_1@meta.data$orig.ident <- "aaUntr_1"

data <- Read10X(data.dir = c("Pool3/part2/part2/outs/per_sample_outs/NM_D2/count/sample_feature_bc_matrix/", "Pool3/part3/part3/outs/per_sample_outs/NM_D2/count/sample_feature_bc_matrix/", "Pool3/part4/part4/outs/per_sample_outs/NM_D2/count/sample_feature_bc_matrix/"), gene.column = 2)
aaUntr_2 = CreateSeuratObject(counts = data$`Gene Expression`, project="aaUntr_2", min.cells=5, min.features=250)
aaUntr_2@meta.data$orig.ident <- "aaUntr_2"

data <- Read10X(data.dir = c("Pool4/part1/part1/outs/per_sample_outs/NM_D3/count/sample_feature_bc_matrix/", "Pool4/part2/part2/outs/per_sample_outs/NM_D3/count/sample_feature_bc_matrix/", "Pool4/part3/part3/outs/per_sample_outs/NM_D3/count/sample_feature_bc_matrix/", "Pool4/part4/part4/outs/per_sample_outs/NM_D3/count/sample_feature_bc_matrix/"), gene.column = 2)
aaUntr_3 = CreateSeuratObject(counts = data$`Gene Expression`, project="aaUntr_3", min.cells=5, min.features=250)
aaUntr_3@meta.data$orig.ident <- "aaUntr_3"

data <- Read10X(data.dir = c("Pool4/part1/part1/outs/per_sample_outs/NM_D4/count/sample_feature_bc_matrix/", "Pool4/part2/part2/outs/per_sample_outs/NM_D4/count/sample_feature_bc_matrix/", "Pool4/part3/part3/outs/per_sample_outs/NM_D4/count/sample_feature_bc_matrix/", "Pool4/part4/part4/outs/per_sample_outs/NM_D4/count/sample_feature_bc_matrix/"), gene.column = 2)
aaUntr_4 = CreateSeuratObject(counts = data$`Gene Expression`, project="aaUntr_4", min.cells=5, min.features=250)
aaUntr_4@meta.data$orig.ident <- "aaUntr_4"

data <- Read10X(data.dir = c("Pool1/part1/part1/outs/per_sample_outs/NM_A1/count/sample_feature_bc_matrix/", "Pool1/part2/part2/outs/per_sample_outs/NM_A1/count/sample_feature_bc_matrix/", "Pool1/part3/part3/outs/per_sample_outs/NM_A1/count/sample_feature_bc_matrix/", "Pool1/part4/part4/outs/per_sample_outs/NM_A1/count/sample_feature_bc_matrix/"), gene.column = 2)
attenuated_1 = CreateSeuratObject(counts = data$`Gene Expression`, project="attenuated_1", min.cells=5, min.features=250)
attenuated_1@meta.data$orig.ident <- "attenuated_1"

data <- Read10X(data.dir = c("Pool1/part1/part1/outs/per_sample_outs/NM_A2/count/sample_feature_bc_matrix/", "Pool1/part2/part2/outs/per_sample_outs/NM_A2/count/sample_feature_bc_matrix/", "Pool1/part3/part3/outs/per_sample_outs/NM_A2/count/sample_feature_bc_matrix/", "Pool1/part4/part4/outs/per_sample_outs/NM_A2/count/sample_feature_bc_matrix/"), gene.column = 2)
attenuated_2 = CreateSeuratObject(counts = data$`Gene Expression`, project="attenuated_2", min.cells=5, min.features=250)
attenuated_2@meta.data$orig.ident <- "attenuated_2"

data <- Read10X(data.dir = c("Pool2/part1/part1/outs/per_sample_outs/NM_A3/count/sample_feature_bc_matrix/", "Pool2/part2/part2/outs/per_sample_outs/NM_A3/count/sample_feature_bc_matrix/", "Pool2/part3/part3/outs/per_sample_outs/NM_A3/count/sample_feature_bc_matrix/", "Pool2/part4/part4/outs/per_sample_outs/NM_A3/count/sample_feature_bc_matrix/"), gene.column = 2)
attenuated_3 = CreateSeuratObject(counts = data$`Gene Expression`, project="attenuated_3", min.cells=5, min.features=250)
attenuated_3@meta.data$orig.ident <- "attenuated_3"

data <- Read10X(data.dir = c("Pool2/part1/part1/outs/per_sample_outs/NM_A4/count/sample_feature_bc_matrix/", "Pool2/part3/part3/outs/per_sample_outs/NM_A4/count/sample_feature_bc_matrix/", "Pool2/part4/part4/outs/per_sample_outs/NM_A4/count/sample_feature_bc_matrix/"), gene.column = 2)
attenuated_4 = CreateSeuratObject(counts = data$`Gene Expression`, project="attenuated_4", min.cells=5, min.features=250)
attenuated_4@meta.data$orig.ident <- "attenuated_4"

data <- Read10X(data.dir = c("Pool1/part1/part1/outs/per_sample_outs/NM_B1/count/sample_feature_bc_matrix/", "Pool1/part2/part2/outs/per_sample_outs/NM_B1/count/sample_feature_bc_matrix/", "Pool1/part3/part3/outs/per_sample_outs/NM_B1/count/sample_feature_bc_matrix/", "Pool1/part4/part4/outs/per_sample_outs/NM_B1/count/sample_feature_bc_matrix/"), gene.column = 2)
adeno_1 = CreateSeuratObject(counts = data$`Gene Expression`, project="adeno_1", min.cells=5, min.features=250)
adeno_1@meta.data$orig.ident <- "adeno_1"

data <- Read10X(data.dir = c("Pool1/part1/part1/outs/per_sample_outs/NM_B2/count/sample_feature_bc_matrix/", "Pool1/part2/part2/outs/per_sample_outs/NM_B2/count/sample_feature_bc_matrix/", "Pool1/part3/part3/outs/per_sample_outs/NM_B2/count/sample_feature_bc_matrix/", "Pool1/part4/part4/outs/per_sample_outs/NM_B2/count/sample_feature_bc_matrix/"), gene.column = 2)
adeno_2 = CreateSeuratObject(counts = data$`Gene Expression`, project="adeno_2", min.cells=5, min.features=250)
adeno_2@meta.data$orig.ident <- "adeno_2"

data <- Read10X(data.dir = c("Pool2/part1/part1/outs/per_sample_outs/NM_B4/count/sample_feature_bc_matrix/", "Pool2/part3/part3/outs/per_sample_outs/NM_B4/count/sample_feature_bc_matrix/", "Pool2/part4/part4/outs/per_sample_outs/NM_B4/count/sample_feature_bc_matrix/"), gene.column = 2)
adeno_4 = CreateSeuratObject(counts = data$`Gene Expression`, project="adeno_4", min.cells=5, min.features=250)
adeno_4@meta.data$orig.ident <- "adeno_4"

data <- Read10X(data.dir = c("Pool3/part1/part1/outs/per_sample_outs/NM_C1/count/sample_feature_bc_matrix/", "Pool3/part2/part2/outs/per_sample_outs/NM_C1/count/sample_feature_bc_matrix/", "Pool3/part3/part3/outs/per_sample_outs/NM_C1/count/sample_feature_bc_matrix/", "Pool3/part4/part4/outs/per_sample_outs/NM_C1/count/sample_feature_bc_matrix/"), gene.column = 2)
mRNA_1 = CreateSeuratObject(counts = data$`Gene Expression`, project="mRNA_1", min.cells=5, min.features=250)
mRNA_1@meta.data$orig.ident <- "mRNA_1"

data <- Read10X(data.dir = c("Pool3/part1/part1/outs/per_sample_outs/NM_C2/count/sample_feature_bc_matrix/", "Pool3/part2/part2/outs/per_sample_outs/NM_C2/count/sample_feature_bc_matrix/", "Pool3/part3/part3/outs/per_sample_outs/NM_C2/count/sample_feature_bc_matrix/", "Pool3/part4/part4/outs/per_sample_outs/NM_C2/count/sample_feature_bc_matrix/"), gene.column = 2)
mRNA_2 = CreateSeuratObject(counts = data$`Gene Expression`, project="mRNA_2", min.cells=5, min.features=250)
mRNA_2@meta.data$orig.ident <- "mRNA_2"

data <- Read10X(data.dir = c("Pool4/part1/part1/outs/per_sample_outs/NM_C3/count/sample_feature_bc_matrix/", "Pool4/part2/part2/outs/per_sample_outs/NM_C3/count/sample_feature_bc_matrix/", "Pool4/part3/part3/outs/per_sample_outs/NM_C3/count/sample_feature_bc_matrix/", "Pool4/part4/part4/outs/per_sample_outs/NM_C3/count/sample_feature_bc_matrix/"), gene.column = 2)
mRNA_3 = CreateSeuratObject(counts = data$`Gene Expression`, project="mRNA_3", min.cells=5, min.features=250)
mRNA_3@meta.data$orig.ident <- "mRNA_3"




#hamster_all <- merge(aaUntr_1, y = c(aaUntr_2, aaUntr_3, aaUntr_4, attenuated_1, attenuated_2, attenuated_3, attenuated_4, adeno_1, adeno_2, adeno_4, mRNA_1, mRNA_2, mRNA_3), add.cell.ids = c("aaUntr_1", "aaUntr_2", "aaUntr_3", "aaUntr_4", "attenuated_1", "attenuated_2", "attenuated_3", "attenuated_4", "adeno_1", "adeno_2", "adeno_4", "mRNA_1", "mRNA_2", "mRNA_3"), project = "impf_N")

hamster_all <- merge(aaUntr_1, y = c(aaUntr_2, aaUntr_3, aaUntr_4, attenuated_1, attenuated_2, attenuated_3, attenuated_4, adeno_1, adeno_2, mRNA_1, mRNA_2, mRNA_3), add.cell.ids = c("aaUntr_1", "aaUntr_2", "aaUntr_3", "aaUntr_4", "attenuated_1", "attenuated_2", "attenuated_3", "attenuated_4", "adeno_1", "adeno_2", "mRNA_1", "mRNA_2", "mRNA_3"), project = "impf_N")


hamster_all@meta.data %>% group_by(orig.ident) %>% tally()


cells_to_keep <- read.table("./exp1_nasal_cells_to_keep.txt", sep="\t") 
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


saveRDS(hamster.integrated, "./seu_exp1_nasal_new_combined_integrated.rds")
