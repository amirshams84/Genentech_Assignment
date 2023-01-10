#####################
# Library Loading
#####################
library(ggplot2)
library(dplyr)
library(tibble)
library(Seurat)
library(cowplot)
library(scales)
library(ggrepel)

##################
#Define constant
##################
COLOR_Pallett = c('#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf')
PLOTLY_Pallett = COLOR_Pallett
FEATURE_PALLETT = PLOTLY_Pallett
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
# PCA
# ------------------------
Seurat_Object <- RunPCA(Seurat_Object, 
                  npcs = 50, 
                  reduction.name = "PCA", 
                  assay = "RNA",  
                  verbose = F
                 )


# ++++++++++++++++++++++++
# NeighborJoining
# ------------------------


Seurat_Object <- FindNeighbors(Seurat_Object, 
    reduction = "PCA", 
    dims = 1:50, 
    assay="RNA", 
    # k.param = 20,
    # compute.SNN = TRUE,
    # prune.SNN = 1/15,
    # nn.method = "rann",
    # annoy.metric = "euclidean",
    # nn.eps = 0,
    # verbose = FALSE,
    # force.recalc = FALSE,
    # do.plot = FALSE
)

##############
#Set cluster resolution
##############
cluster_resolution = 1.3
cluster_label = paste("RNA_snn_res.", cluster_resolution, sep="")

##############
#Clustering
##############
Seurat_Object <- FindClusters(Seurat_Object,
                              resolution =cluster_resolution,
                              #algorithm = "Leiden"
                              verbose = FALSE)
#head(Seurat_Object)
levels(Seurat_Object@meta.data[[cluster_label]])


# ++++++++++++++++++++++++
# UMAP
# ------------------------
Seurat_Object <- RunUMAP(Seurat_Object, 
                         reduction = "PCA", 
                         reduction.name = "UMAP",
                         assay = "RNA",
                         dims=1:30,
                         verbose = F,
                         #features = PCA_PC1_Gene_Loading_List
                        )
# ++++++++++++++++++++++++
# before trimming
# ------------------------


cluster_length = levels(Seurat_Object@meta.data[[cluster_label]])
# length(cluster_length)
new_cluster_name_List = sprintf("cluster_%01d", seq(1,length(cluster_length)))
#new_cluster_name_List
#################
#Relabel clusters
#################
Seurat_Object <- SetIdent(Seurat_Object, value = cluster_label)
#
old_cluster_name_List = levels(Seurat_Object@meta.data[[cluster_label]])
names(new_cluster_name_List) <- old_cluster_name_List
#
Seurat_Object <- RenameIdents(Seurat_Object, new_cluster_name_List)
#head(alldata@meta.data)

cluster_label_DF = as.data.frame(Seurat_Object@active.ident) 
colnames(cluster_label_DF)[1] <- "cluster_label"
# head(cluster_label_DF)
####################
#FEATURE ORDER/COLOR
####################
Feature_Order_List = new_cluster_name_List
Feature_List_Size = length(Feature_Order_List)
Feature_Label_Dict = Feature_Order_List
names(Feature_Label_Dict) = Feature_Order_List
Feature_Order_Factor = factor(Feature_Order_List, levels=Feature_Order_List)
if(Feature_List_Size > 10){
    Feature_Color_List = colorRampPalette(COLOR_Pallett)(Feature_List_Size)
}else{
    Feature_Color_List = COLOR_Pallett[1:Feature_List_Size]
}
Feature_Color_Dict = Feature_Color_List
names(Feature_Color_Dict) = Feature_Order_List
#################
#set specifc Color  of each clusters
#################
color_alldata <- RenameIdents(Seurat_Object, Feature_Color_Dict)
# head(alldata@meta.data)

cluster_color_DF = as.data.frame(color_alldata@active.ident) #%>% 
    #rownames_to_column(var="UMI") %>%
    #rename(cluster_color = `color_alldata@active.ident`)
cluster_color_DF = as.data.frame(color_alldata@active.ident) 
colnames(cluster_color_DF)[1] <- "cluster_color"
# head(cluster_color_DF)
###################
#update metadata
###################
Seurat_Object <- AddMetaData(
  object = Seurat_Object,
  metadata = cluster_label_DF
)
Seurat_Object <- AddMetaData(
  object = Seurat_Object,
  metadata = cluster_color_DF
)
#head(Seurat_Object@meta.data)



######################
# Extract UMAP from Seurat object
######################
UMAP_DF = as.data.frame(Seurat_Object@reductions$UMAP@cell.embeddings)
# head(UMAP_DF)
###################
#update metadata with UMAP column
###################
Seurat_Object <- AddMetaData(
  object = Seurat_Object,
  metadata = UMAP_DF
)
# head(Seurat_Object@meta.data)
#################
#Extracted dataframe
#################
metadata_UMAP_DF = as.data.frame(Seurat_Object@meta.data) %>%
    rownames_to_column(var="UMI") 
######################
# Centroids of each cluster
######################
selected_metadata_UMAP_DF = metadata_UMAP_DF %>%
    dplyr::select(UMI, cluster_label, UMAP_1, UMAP_2)
# head(selected_metadata_tSNE_DF)
cluster_centroids = as.data.frame(selected_metadata_UMAP_DF) %>%
    group_by(cluster_label) %>%
    summarize(UMAP_1_centroids = median(UMAP_1), UMAP_2_centroids = median(UMAP_2))
#head(cluster_centroids)


####################
# Cluster ScatterPlot
####################
Cluster_ScatterPlot_Object = ggplot(metadata_UMAP_DF, aes(x=UMAP_1, y=UMAP_2, color=cluster_label)) +
    geom_point(size=1)  +
    geom_label_repel(data = cluster_centroids, aes(label =cluster_label ,x = UMAP_1_centroids, y = UMAP_2_centroids), 
               show.legend = FALSE, size=3,point.padding = NA, color="black"
              ) +
    #xlim(min_exp, max_exp) +
    #ylim(min_exp, max_exp) +
    labs(title = NULL, subtitle = NULL, caption = NULL, tag = NULL, x="UMAP-1", y="UMAP-2") +
    scale_color_manual(breaks = Feature_Order_List, labels=Feature_Label_Dict, values=Feature_Color_Dict) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    #annotate(geom="segment", x = -Inf, y = -Inf, xend=-Inf, yend=5, arrow = arrow(length = unit(10, "mm"))) +
    theme_bw() + 
    theme(
        ###LEGEND
        #aspect.ratio=1,
        legend.position="right",
        legend.title = element_blank(),
        legend.key=element_blank(),
        #legend.key.size = unit(2,"line"),
        legend.text.align = 0,
        legend.text = element_text(size=10,  family="sans", face="plain", color = "black", angle=0, vjust=0.5),
        legend.background = element_blank(),
        #legend.box.margin=LEGEND_BOX_MARGIN,
        #legend.margin=LEGEND_MARGIN,
        ###
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #plot.margin = SCATTER_MARGIN,
        axis.line = element_line(colour = "black"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(size=10,  family="sans", face="plain", color = "black", angle=0, hjust =0),
        axis.title.y = element_text(size=10,  family="sans", face="plain", color = "black", angle=90, hjust =0),
        axis.line.x.top = element_line(size=0.1),
        axis.line.x.bottom = element_line(size=0.1),
        axis.line.y.left = element_line(size=0.1),
        axis.line.y.right = element_line(size=0.1),
        #
        
        
    )





png(
    filename= paste(output_dir_path, Sample_name, "_Seurat_Object_UMAP.png", sep=""),
    width = 10, 
    height = 5,
    units="in",res=300
    )

    Cluster_ScatterPlot_Object
dev.off()



# ++++++++++++++++++++++++
# Save RDS object
# ------------------------
saveRDS(Seurat_Object, 
    paste(output_dir_path, Sample_name, "_Seurat_Object_processed_UMAP.rds", sep="")
   )

# ++++++++++++++++++++++++
# Finito message
# ------------------------
print(
    paste("PCA_UMAP for ", Sample_name, " completed Successfully.!!!", sep="")
    )