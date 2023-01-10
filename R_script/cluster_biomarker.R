#####################
# Library Loading
#####################
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
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
# Find Markers
# ------------------------
Seurat_Object <- SetIdent(Seurat_Object, value = "cluster_label")
markers_genes <- FindAllMarkers(Seurat_Object,
                               test.use = "wilcox",
                               only.pos = FALSE,
                               assay = "RNA",
                               verbose=TRUE)


#
# Process gene ranke table
# 
rand_pvalue<- function(){
    
    runif(n=1000,min=1,max=9)
    }
processed_celltype_biomarker_genes = markers_genes %>%
    mutate(p_val_adj = ifelse(p_val_adj==0, rand_pvalue()* 1e-100 ,p_val_adj )) %>%
    mutate(p_val = ifelse(p_val==0, rand_pvalue()* 1e-300,p_val )) %>%
    rownames_to_column(var="Gene_ID")
#head(processed_celltype_biomarker_genes, n=10L)

top_limit=10
top_biomarkers = processed_celltype_biomarker_genes %>% 
                    group_by(cluster) %>% 
                    filter(!str_detect(Gene_ID, "^RP[SL]")) %>%
                    filter(!str_detect(Gene_ID, "^MT-")) %>%
                    arrange(p_val_adj, -abs(avg_log2FC)) %>%
                    slice(1:top_limit) %>% 
                    arrange(-pct.1)
#head(top_biomarkers)
top_biomarkers_List = top_biomarkers %>% pull(gene) %>% unique()
#length(top_biomarkers_List)








png(
    filename= paste(output_dir_path, Sample_name, "_Seurat_Object_cluster_biomarkers_heatmap.png", sep=""),
    width = 15, 
    height = 5,
    units="in",res=300
    )
    DoHeatmap(features = top_biomarkers_List, Seurat_Object,group.by = "cluster_label", assay="RNA", size=3)

dev.off() 


# #################
#
# #################


cluster_cluster_expression_dotplot = DotPlot(Seurat_Object, 
        features = top_biomarkers_List,
        group.by = "cluster_label",
        cols = c( "grey", "blue"),
        assay = "RNA",
        scale.by = "radius"
       ) +
scale_x_discrete(position = "top") + 
#scale_y_discrete(breaks = Feature_Order_List, labels=Feature_Label_Dict, limits=Feature_Order_List, values=Feature_Color_Dict) +
#scale_y_continuous(breaks = Feature_Order_List, labels=Feature_Label_Dict, limits=Feature_Order_List, values=Feature_Color_Dict) + 
labs(title= NULL, x =NULL, y=NULL, size="% expressed \n in group", colour ="Average\nexpression") +
theme(
            text = element_text(size=8, face="plain"),
            plot.title = element_text(size=10, color = "black", face = "plain", hjust = 0),
            plot.subtitle = element_text(size=10, color = "black", face = "plain", hjust = 0),
            plot.caption = element_text(size=10, color = "black", face = "plain"),
        axis.text.x = element_text(size=8, color = "black", face = "bold", angle=45, hjust=0),
        axis.text.y = element_text(size=8, color = "black", face = "bold"),
        axis.line = element_line(colour = "black", size=0.0),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
        ) +
guides(
    colour = guide_colorbar(reverse = FALSE, order = 1), 
    size = guide_legend(reverse = TRUE, order = 2), 
    shape = guide_legend(reverse = TRUE)
) 





png(
    filename= paste(output_dir_path, Sample_name, "_Seurat_Object_cluster_biomarkers_dotplot.png", sep=""),
    width = 15, 
    height = 5,
    units="in",res=300
    )
    cluster_cluster_expression_dotplot

dev.off() 



# ++++++++++++++++++++++++
# Finito message
# ------------------------
print(
    paste("Cluster_biomarker for ", Sample_name, " completed Successfully.!!!")
    )