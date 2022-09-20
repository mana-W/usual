library(Seurat)

#scmap
library(SingleCellExperiment)
library(scmap)

args1 <- commandArgs()
args <- args1[c(6:length(args1))]
print(args)
#args1 refrds ; args2 colname of annotation ; args3 userds
rds1 <- readRDS(args[1])
rds2 <- readRDS(args[3])

ann1 <- data.frame(cell_type1=rds1@meta.data[args[2]])
names(ann1) <- "cell_type1"
row.names(ann1) <- row.names(rds1@meta.data)
exp <- data.frame(rds1@assays$RNA@counts)
sce1 <- SingleCellExperiment(assays = list(normcounts = as.matrix(exp)),
                                        colData = ann1)
logcounts(sce1) <- log2(normcounts(sce1) + 1)
rowData(sce1)$feature_symbol <- rownames(sce1)
sce1 <- sce1[!duplicated(rownames(sce1)), ]
ann1 <- data.frame(cell_type1=rds2@meta.data$orig.ident)
row.names(ann1) <- row.names(rds2@meta.data)
exp <- data.frame(rds2@assays$RNA@counts)
sce2 <- SingleCellExperiment(assays = list(normcounts = as.matrix(exp)),
                                                                                colData = ann1)
logcounts(sce2) <- log2(normcounts(sce2) + 1)
rowData(sce2)$feature_symbol <- rownames(sce2)
sce2 <- sce2[!duplicated(rownames(sce2)), ]

sce1 <- selectFeatures(sce1, suppress_plot = FALSE,n_features = 1000)
sce2 <- selectFeatures(sce2, suppress_plot = FALSE,n_features = 1000)
sce1 <- indexCell(sce1)
sce2 <- indexCell(sce2)

sce2_to_sce1 <- scmapCell(
  projection = sce2,
  index_list = list(
    seger = metadata(sce1)$scmap_cell_index
  ),
  w = 50
)

cell_type_scmap <- colData(sce1)$cell_type1[sce2_to_sce1$seger[[1]][1,]]
anndata <- data.frame(cellbarcode= colnames(sce2),
                                                cell_type_scmap=cell_type_scmap)

save.image("image.rda")

library(tidyverse)
library(SummarizedExperiment)
library(scuttle)
library(SingleR)

p_count=rds1[["RNA"]]@counts
rds1@meta.data$Index <- row.names(rds1@meta.data)
pdata=rds1@meta.data[,c("Index",args[2])]
rownames(pdata)=pdata$Index
pdata$Index=NULL
colnames(pdata)="ref_label"

pdata_SE <- SummarizedExperiment(assays=list(counts=p_count),colData = pdata)
pdata_SE <- logNormCounts(pdata_SE)

my_count=rds2[["RNA"]]@counts
common_gene <- intersect(rownames(my_count), rownames(pdata_SE))
pdata_SE <- pdata_SE[common_gene,]
my_count <- my_count[common_gene,]

my_SE <- SummarizedExperiment(assays=list(counts=my_count))
my_SE <- logNormCounts(my_SE)

singleR_res <- SingleR(test = my_SE,
                       ref = pdata_SE,
                       labels = pdata_SE$ref_label)

anno_df <- as.data.frame(singleR_res$labels)
anno_df$cellbarcode <- rownames(singleR_res)
colnames(anno_df)[1] <- 'cell_type_singleR'

anno <- merge(anndata,anno_df,by="cellbarcode")
save.image("image.rda")
write.csv(anno,"annores.csv",row.names=F,quote=F)
rds2@meta.data$cell_type_singleR <- anno$cell_type_singleR[match(row.names(rds2@meta.data),anno$cellbarcode)]
rds2@meta.data$cell_type_scmap <- anno$cell_type_scmap[match(row.names(rds2@meta.data),anno$cellbarcode)]
saveRDS(rds2,"anno_rds.rds")



