#####################
# Library Loading
#####################
#suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))


# ++++++++++++++++++++++++
# I/O
# -----------------------
args = commandArgs(trailingOnly=TRUE)
input_file_Path = args[1]
output_dir_path  = args[2]
Sample_name  = args[3]


dir.create(output_dir_path, showWarnings = FALSE)



# +++++++++++++++++++++++
# Build seurat object
# -----------------------
h5_Object = Read10X_h5(input_file_Path)
Seurat_Object <- CreateSeuratObject(h5_Object, project = Sample_name)


# ++++++++++++++++++++++++
# Save RDS object
# ------------------------
saveRDS(Seurat_Object, 
    paste(output_dir_path, Sample_name, "_Seurat_Object.rds", sep="")
   )

# ++++++++++++++++++++++++
# Finito message
# ------------------------
print(
    paste("build seurat object for ", Sample_name, " completed Successfully.!!!")
    )