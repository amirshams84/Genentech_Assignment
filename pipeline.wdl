# WORKFLOW DEFINITION 
workflow pipeline {

  Array[File] h5_files
   
  scatter(s in h5_files) {
    call build_seurat_object { input: h5_file=s }
    call QC_trimming_normalize { input: Seurat_object = build_seurat_object.Seurat_Object }
    call PCA_UMAP { input: Seurat_object_processed = QC_trimming_normalize.Seurat_object_processed }
    call cluster_biomarker { input: Seurat_object_processed_UMAP = PCA_UMAP.Seurat_object_processed_UMAP }

  }

}


task build_seurat_object {
  File h5_file
  String prefix = sub(h5_file, ".h5", "")

  command <<<
    set +eu && PS1=dummy && . /opt/conda/etc/profile.d
    conda activate Seurat_conda
    
    Rscript /script/Genentech_Assignment/R_script/build_seurat_object.R ${h5_file} /result ${prefix}
  >>>
  output {
    File Seurat_Object  = "/result/${prefix}_Seurat_Object.rds"
  }
  runtime {}
}


task QC_trimming_normalize {
  File Seurat_object
  String prefix = sub(Seurat_object, "_Seurat_Object.rds", "") 

  command <<<
    set +eu && PS1=dummy && . /opt/conda/etc/profile.d
    conda activate Seurat_conda
    Rscript /script/Genentech_Assignment/R_script/QC_trimming_normalize.R ${Seurat_object} /result ${prefix}
  >>>
  output {
    File Seurat_object_processed  = "/result/${prefix}_Seurat_Object_processed.rds"
  }
  runtime {}
}


task PCA_UMAP {
  File Seurat_object_processed
  String prefix = sub(Seurat_object_processed, "_Seurat_Object_processed.rds", "")  

  command <<<
    set +eu && PS1=dummy && . /opt/conda/etc/profile.d
    conda activate Seurat_conda
    Rscript /script/Genentech_Assignment/R_script/PCA_UMAP.R ${Seurat_object_processed} /result ${prefix}
  >>>
  output {
    File Seurat_object_processed_UMAP  = "/result/${prefix}_Seurat_Object_processed_UMAP.rds"
  }
  runtime {}
}


task cluster_biomarker {
  File Seurat_object_processed_UMAP
  String prefix = sub(Seurat_object_processed_UMAP, "_Seurat_Object_processed_UMAP.rds", "") 

  command <<<
    set +eu && PS1=dummy && . /opt/conda/etc/profile.d
    conda activate Seurat_conda
    Rscript /script/Genentech_Assignment/R_script/cluster_biomarker.R ${Seurat_object_processed_UMAP} /result ${prefix}
  >>>
  output {
    File Seurat_object_processed_UMAP_cluster = "/result/${prefix}_Seurat_Object_cluster_biomarkers_dotplot.png"
  }
  runtime {}
}