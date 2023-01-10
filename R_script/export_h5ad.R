#####################
# Library Loading
#####################
suppressPackageStartupMessages(library(sceasy))
suppressPackageStartupMessages(library(Seurat))

# ++++++++++++++++++++++++
# I/O
# -----------------------
args = commandArgs(trailingOnly=TRUE)
input_file_Path = args[1]
output_dir_path  = args[2]
Sample_name  = args[3]


dir.create(output_dir_path, showWarnings = FALSE)



# ++++++++++++++++++++++++
# LOAD DATA
# ------------------------
Seurat_Object<- readRDS(input_file_Path)


# ++++++++++++++++++++++++
# Save RDS object
# ------------------------

sceasy::convertFormat(Seurat_Object, from="seurat", to="anndata", outFile=paste(output_dir_path, Sample_name, "_Seurat_Object_processed_UMAP.h5ad", sep=""))
# ++++++++++++++++++++++++
# Finito message
# ------------------------
print(
    paste("export h5ad for ", Sample_name, " completed Successfully.!!!")
    )