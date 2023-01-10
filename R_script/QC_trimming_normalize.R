#####################
# Library Loading
#####################
#suppressPackageStartupMessages(library(tidyverse))
library(ggplot2)
library(dplyr)
library(Seurat)
library(cowplot)
library(scales)


# ++++++++++++++++++++++++
# I/O
# -----------------------
args = commandArgs(trailingOnly=TRUE)
input_file_Path = args[1]
output_dir_path  = args[2]
Sample_name  = args[3]
# min_nFeature = args[4]
# max_nFeature = args[5]
# max_Mito = args[6]
# min_express_cell = args[7]


dir.create(output_dir_path, showWarnings = FALSE)


# ++++++++++++++++++++++++
# LOAD DATA
# ------------------------
Seurat_Object<- readRDS(input_file_Path)


# ++++++++++++++++++++++++
# calculate mitochondrial and Ribosomal 
# ------------------------
Seurat_Object = PercentageFeatureSet(Seurat_Object, "^MT-", col.name = "percent_mito")
Seurat_Object = PercentageFeatureSet(Seurat_Object, "^RP[SL]", col.name = "percent_ribo")


# ++++++++++++++++++++++++
# before trimming
# ------------------------

feats <- c("nFeature_RNA","nCount_RNA","percent_mito","percent_ribo")

png(
    filename= paste(output_dir_path, Sample_name, "_Seurat_Object_QC_before_trimming.png", sep=""),
    width = 15, 
    height = 5,
    units="in",res=300
    )

    VlnPlot(Seurat_Object, group.by= "orig.ident", features = feats, pt.size =0, ncol = 4) + NoLegend()
dev.off()


# ++++++++++++++++++++++++
# the trimming
# ------------------------

Seurat_Object <- subset(Seurat_Object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mito < 5)

# ++++++++++++++++++++++++
# after trimming
# ------------------------

feats <- c("nFeature_RNA","nCount_RNA","percent_mito","percent_ribo")

png(
    filename= paste(output_dir_path, Sample_name, "_Seurat_Object_QC_after_trimming.png", sep=""),
    width = 15, 
    height = 5,
    units="in",res=300
    )

    VlnPlot(Seurat_Object, group.by= "orig.ident", features = feats, pt.size =0, ncol = 4) + NoLegend()
dev.off()



# ++++++++++++++++++++++++
# log normlaize
# ------------------------
Seurat_Object <- NormalizeData(
    object = Seurat_Object,
    normalization.method = "LogNormalize",
    scale.factor = 10000
    )

Seurat_Object <- ScaleData(Seurat_Object, verbose = FALSE)


suppressWarnings(
    suppressMessages(
        Seurat_Object <- FindVariableFeatures(
            Seurat_Object, mean.function = ExpMean, dispersion.function = LogVMR, 
                x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,verbose = FALSE,assay = "RNA"
                )
            )
        )


top20 <- head(VariableFeatures(Seurat_Object), 20)
png(
    filename= paste(output_dir_path, Sample_name, "_Seurat_Object_QC_variable_genes.png", sep=""),
    width = 15, 
    height = 5,
    units="in", res=300
    )
    LabelPoints(plot = VariableFeaturePlot(Seurat_Object), points = top20, repel = TRUE)

dev.off() 



all.genes <- rownames(Seurat_Object)
Seurat_Object <- ScaleData(Seurat_Object, features = all.genes, vars.to.regress = "percent_mito")

# ++++++++++++++++++++++++
# Save RDS object
# ------------------------
saveRDS(Seurat_Object, 
    paste(output_dir_path, Sample_name, "_Seurat_Object_processed.rds", sep="")
   )

# ++++++++++++++++++++++++
# Finito message
# ------------------------
print(
    paste("QC_trimming_normalize for ", Sample_name, " completed Successfully.!!!")
    )