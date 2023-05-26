################################################
#Figure2b; pseudotime analysis
################################################

# Load packages
library(SingCellaR)
library(monocle3)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(plyr)
library(gg3D)

# Load SingCellaR object
load(file = "../../../Sources/AP_HTMPNAML_Harmony.revised.rdata")

# Check clusters on UMAP
plot_umap_label_by_clusters(AP.HTMPNAML,show_method="merged_louvain",mark.clusters = T)

#cl1 - LSC
#cl2 - erythroid
#cl3 - preleukemic

# Prepare input files for Monocle3
    # Expression matrix
    cells.used <- AP.HTMPNAML@sc.clusters$Cell
    umi <- get_umi_count(AP.HTMPNAML)
    used.umi <- umi[,cells.used]
    expression_matrix <- used.umi
    dim(expression_matrix) # check the dimension of object

    # Cell cluster metadata
    cell_metadata <- AP.HTMPNAML@sc.clusters
    rownames(cell_metadata) <- cell_metadata$Cell

    # Gene metadata
    gene_annotation <- as.data.frame(rownames(used.umi))
    colnames(gene_annotation) <- "gene_short_name"
    rownames(gene_annotation) <- gene_annotation$gene_short_name

# Create Monocle3 object.
cds <- new_cell_data_set(expression_data = expression_matrix,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)

# Integrate Monocle3 and SingCellaR results
    # Pre-process Monocle3 object
    cds <- preprocess_cds(cds,num_dim = 100,method = "PCA")
    cds <- align_cds(cds)

    # Substitute Monocle3’s embeddings with SingCellaR’s embeddings
    embeddings <- AP.HTMPNAML@Harmony.embeddings
    cds@int_colData@listData$reducedDims$Aligned <- embeddings

    # Nonlinear dimension reduction
    cds <- reduce_dimension(cds,
                              reduction_method = "UMAP",
                              umap.min_dist = 0.3,
                              preprocess_method = "Aligned")

    # Identify and assign clusters
    cds <- cluster_cells(cds,
                           reduction_method = "UMAP",
                           k = 30,
                           cluster_method = "louvain")

    # Project diffusion map as 2D
        # Retrieve coordinates
        DiffMap <- AP.HTMPNAML@diffusionmap.result
        DiffMap <- DiffMap[,c("Cell", "DC1", "DC2")]

        # Multiply coordinates by scale factor to match UMAP range
        fivenum(AP.HTMPNAML@umap.result$UMAP1)
        fivenum(DiffMap$DC1)
        DiffMap$DC1 <- DiffMap$DC1*100
        DiffMap$DC2 <- DiffMap$DC2*100

        # Annotate cluster IDs
        DiffMap <- join(DiffMap,
                        AP.HTMPNAML@sc.clusters[,c("Cell", "merged_louvain")],
                        by="Cell",
                        type="left"
                        )

        # Assign color
        DiffMap$color <- NA
        DiffMap$color[which(DiffMap$merged_louvain=="cl3")] <- "blue" # preLSC
        DiffMap$color[which(DiffMap$merged_louvain=="cl1")] <- "green" # LSC
        DiffMap$color[which(DiffMap$merged_louvain=="cl2")] <- "red" # erythroid

        # Quick plot
        plot(DiffMap$DC1, DiffMap$DC2, col=DiffMap$color, cex=0.1)
    
        # Update UMAP slot with 2D diffusion map
        UMap <- AP.HTMPNAML@umap.result[,c("Cell","UMAP1","UMAP2")]
        DiffMap <- join(UMap[,"Cell",drop=FALSE], DiffMap, by="Cell", type="left")
        table(DiffMap$Cell==UMap$Cell) # TRUE
        DiffMap$merged_louvain <- NULL
        DiffMap$color <- NULL
        names(DiffMap)[which(names(DiffMap)=="DC1")] <- "UMAP1"
        names(DiffMap)[which(names(DiffMap)=="DC2")] <- "UMAP2"
        AP.HTMPNAML@umap.result <- DiffMap
    
    # Substitute Monocle3’s UMAP embeddings with SingCellaR’s embedding
    newcds<- cds # change monocle3 objects name

    #SingCellaR.umap <-, c("Cell","UMAP1","UMAP2")]
    #AP.HTMPNAML@umap.result[monocle3.umap] <- newcds@int_colData$reducedDims$UMAP
    SingCellaR.umap <- AP.HTMPNAML@umap.result[,c("Cell","UMAP1","UMAP2")]
    monocle3.umap <- newcds@int_colData$reducedDims$UMAP
    umap <- SingCellaR.umap[match(rownames(monocle3.umap),SingCellaR.umap$Cell),]
    rownames(umap) <- umap$Cell
    umap$Cell <- NULL
    newcds@int_colData$reducedDims$UMAP <- umap

    # Substitute Monocle3’s cluster identity with SingCellaR’s cluster identity
    #anno.clusters <- AP.HTMPNAML@sc.clusters$louvain_cluster
    anno.clusters <- AP.HTMPNAML@sc.clusters$merged_louvain
    names(anno.clusters) <- AP.HTMPNAML@sc.clusters$Cell
    table(names(newcds@clusters$UMAP$clusters)==names(anno.clusters)) # TRUE
    newcds@clusters$UMAP$clusters <- anno.clusters

#########################################################################
####################### PSEUDOTIME ANALYSIS #############################
#########################################################################

# Generate trajectory graph and order cells by pseudotime
    # Create graph
    newcds <- learn_graph(newcds)

    # Save Monocle object
    save(newcds, file="../../Sources/AP.HTMPNAML_monocle3.rdata")

    # Load Monocle object
    # load(file="../../Sources/AP.HTMPNAML_monocle3.rdata")

    # Apply function to retrieve root node
    #root.nodes <- get_earliest_principal_node(newcds,cluster = "cl1")

    # Order cells by pseudotime relative to root node
    #newcds <- order_cells(newcds, root_pr_nodes = root.nodes)
    newcds <- order_cells(newcds) # Choose starting point at preLSC cluster
                                  # Please save screenshot

# Diffusion map
    # Retrieve dfm coordinates
    df <- AP.HTMPNAML@diffusionmap.result
    dim(df)

    # Annotate pseudotime
    new_data_cl8 <- data.frame(pseudotime = pseudotime(newcds,reduction_method = "UMAP"))
    new_data_cl8$Cell <- rownames(new_data_cl8)
    dim(new_data_cl8)

    df.cl8 <- merge(df,new_data_cl8,by = "Cell",all.x = TRUE)

    head(df.cl8)

    df.cl8 <- df.cl8[order(df.cl8$pseudotime), ]

    # tiff(file = "monocle_mannually_diffusionmap3D_cl8.tiff",units = "in",width =6, height = 5, res=300,compression = "lzw")
    qplot(x=0, y=0, z=0, geom="blank") +
      axes_3D() +
      theme_void()

    ggplot(df.cl8, aes(x=-DC1, y=-DC2, z=-DC4, color=pseudotime)) +
      theme_void() +
      axes_3D() +
      scale_color_viridis_c(name = 'Pseudotime',option = "C", values=c(0,0.6,0.9,1))+
      # geom_point(size=0.03) +
      stat_3D()

    ggsave(filename="Figure2b.png",width = 8,height = 8)
