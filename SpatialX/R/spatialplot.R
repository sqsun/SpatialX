#####################################################################
## Package: SPATIALX, 
## Package includes: SpatialVG, for single-slice
## Version: 0.1.0
## Date: 2020-2-14
## Modified: 2020-2-22 04:50:39;2020-3-14 19:55:43;2022-6-4 09:26:03;2022-9-17 09:37:51;
## 2024-8-13 21:45:39
## Title : Spatially-aware expression analysis of spatially resolved transcriptomics
## Determining the spatially variable genes using Poisson mixed models and 
## spatially variable genes using logistic mixed models
## Authors: S.Q. Sun
## Contacts: sqsunsph@xjtu.edu.cn
##          Xi'an Jiaotong University, 
##          Center for Single-Cell Omics and Health, Department of Biostatistics
######################################################################



#' Spatial Visualization for Spatial Transcriptomics Data
#'
#' This function creates various types of spatial plots for visualizing spatial
#' transcriptomics data, including gene expression patterns, clustering results,
#' dimensionality reductions, and deconvolution results.
#'
#' @param object The data to visualize. Can be a list of data frames, an Assay
#' object, or a SpatialX object.
#' @param plot.type Type of plot to create. Available options:
#' \itemize{
#'   \item{"dot"}{ - Dot plot for continuous values (e.g., gene expression)}
#'   \item{"rgb"}{ - RGB plot for multi-dimensional data (e.g., t-SNE, UMAP)}
#'   \item{"vln"}{ - Violin plot (not yet implemented)}
#'   \item{"prop"}{ - Proportion/weight plot for deconvolution results}
#' }
#' @param point.size Size of points in the plot (default: 2).
#' @param text.size Size of text elements in the plot (default: 15).
#' @param colors Color palette for the plot. For dot plots, a gradient from low
#' to high values. For RGB plots, ignored.
#' @param title Plot title. If NULL, uses default titles based on data.
#' @param legend Position of legend ("none" to hide).
#' @param verbose Whether to display progress information (default: TRUE).
#' @param slot.use Assay slot to use for data extraction (default: "data").
#' @param plot.data Type of data to plot. Available options:
#' \itemize{
#'   \item{"features"}{ - Plot gene expression for specified features}
#'   \item{"cluster"}{ - Plot clustering results}
#'   \item{"reduction"}{ - Plot dimensionality reduction results (PCA)}
#'   \item{"tsne"}{ - Plot t-SNE results as RGB}
#'   \item{"umap"}{ - Plot UMAP results as RGB}
#'   \item{"deconvolution"}{ - Plot cell type deconvolution results}
#' }
#' @param location Spatial coordinates matrix. If NULL, extracted from object.
#' @param features Vector of feature names to plot when plot.data = "features".
#' @param dims.plot Dimensions to plot when plot.data = "reduction" (default: 1:2).
#' @param assay Name of assay to use for SpatialX objects.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A list of ggplot objects that can be further customized or combined.
#'
#' @details
#' This function provides comprehensive spatial visualization capabilities for
#' spatial transcriptomics data analysis:
#'
#' \itemize{
#'   \item{\strong{Dot Plots}:} Display continuous values (gene expression,
#'   PCA scores) as colored dots on spatial coordinates.
#'   \item{\strong{RGB Plots}:} Visualize multi-dimensional data (t-SNE, UMAP)
#'   by mapping dimensions to RGB color channels.
#'   \item{\strong{Proportion Plots}:} Show cell type proportions from
#'   deconvolution analysis.
#' }
#'
#' @section Visualization Types:
#' \describe{
#'   \item{Gene Expression}{Plot expression levels of specific genes across
#'   spatial locations using dot plots.}
#'   \item{Clustering Results}{Visualize spatial distribution of cell clusters.}
#'   \item{Dimensionality Reduction}{Plot PCA components or embed t-SNE/UMAP
#'   results in spatial context.}
#'   \item{Deconvolution}{Display estimated cell type proportions across
#'   spatial locations.}
#' }
#'
#' @section Color Schemes:
#' \itemize{
#'   \item{\strong{Dot Plots}:} Uses a diverging color gradient (orange-grey-green)
#'   by default, centered at the mean value.
#'   \item{\strong{RGB Plots}:} Maps the first three dimensions of embeddings
#'   to red, green, and blue channels respectively.
#'   \item{\strong{Custom Colors}:} Users can specify custom color palettes
#'   through the \code{colors} parameter.
#' }
#'
#' @examples
#' \dontrun{
#' # Create example spatial data
#' spatial_data <- list(
#'   gene1 = data.frame(
#'     x = runif(100),
#'     y = runif(100),
#'     v = rnorm(100)
#'   )
#' )
#'
#' # Create dot plot
#' plots <- SpatialPlot(object = spatial_data, plot.type = "dot")
#'
#' # Create RGB plot (for multi-dimensional data)
#' rgb_data <- list(
#'   embedding = data.frame(
#'     x = runif(100),
#'     y = runif(100),
#'     r = runif(100),
#'     g = runif(100),
#'     b = runif(100)
#'   )
#' )
#' plots <- SpatialPlot(object = rgb_data, plot.type = "rgb")
#'
#' # With SpatialX object
#' data("pbmc_small")
#'
#' # Plot specific genes
#' plots <- SpatialPlot(object = object,
#'                      plot.data = "features",
#'                      features = c("Gene1", "Gene2"))
#'
#' # Plot clustering results
#' plots <- SpatialPlot(object = object,
#'                      plot.data = "cluster")
#'
#' # Plot PCA results
#' plots <- SpatialPlot(object = object,
#'                      plot.data = "reduction",
#'                      dims.plot = 1:3)
#'
#' # Plot t-SNE as RGB
#' plots <- SpatialPlot(object = object,
#'                      plot.data = "tsne")
#'
#' # Customize appearance
#' plots <- SpatialPlot(object = object,
#'                      point.size = 3,
#'                      text.size = 12,
#'                      colors = c("blue", "white", "red"),
#'                      title = "Custom Plot")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_identity ggtitle
#' @importFrom ggplot2 theme_void theme element_text scale_color_gradient2
#' @importFrom grDevices rgb
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap
#' @rdname SpatialPlot
#' @concept visualization
#' @concept spatial
#' @export
#' @method SpatialPlot default
#' 
SpatialPlot.default <- function(object,
                                plot.type = "dot",
                                point.size = 2,
                                text.size = 15,
                                colors = c("#4E84C4", "#FFDB6D"),
                                title = NULL,
                                legend = "none",
                                verbose = TRUE, ...) {
	ggp <- list()
	x <- y <- r <- g <- b <- v <- NULL
	## 
	if(verbose){cat(paste0("SpatialPlot::draw the ", plot.type, " plot. \n"))}
	## main functions
	if(plot.type == "rgb"){
	##******************##
	##   RGB plots      ##
	##******************##
	  for(ip in 1:length(object)){
	    ggp[[ip]] <- ggplot2::ggplot(data=object[[ip]], ggplot2::aes(x = .data$x, 
	                                                                 y = .data$y, 
	                                                                 col = rgb(.data$r, .data$g, .data$b)))+
	      ggplot2::geom_point(size=point.size)+ggplot2::scale_color_identity()+
	      ggplot2::ggtitle(ifelse(is.null(title), names(object)[ip], title))+
	      ggplot2::theme_void()+
	      ggplot2::theme(plot.title = ggplot2::element_text(size = text.size),text = ggplot2::element_text(size = text.size),legend.position = "bottom")
	  }## end for
	  
	}else if(plot.type == "dot"){
	##*****************##
	##   Dot plots     ##
	##*****************##
	  for(ip in 1:length(object)){
	      mid <- mean(object[[ip]]$v)
	      #colors <- colorRampPalette(c("grey92", "#32a676","#32a676"))(3)
	      p <- ggplot2::ggplot(data=object[[ip]], ggplot2::aes(x = .data$x, 
	                                                           y = .data$y, 
	                                                           color = .data$v))+
	      ggplot2::geom_point(alpha = 1, size=point.size)+
	        ggplot2::scale_color_gradient2(midpoint = mid, low = "#dd8a0b", mid = "grey92", high = "#32a676")+
	        ggplot2::ggtitle(ifelse(is.null(title), names(object)[ip], title))+
	        ggplot2::theme_void()+
	        ggplot2::theme(plot.title = ggplot2::element_text(size = text.size,face = "bold"),
	              text = ggplot2::element_text(size = text.size),legend.position = "bottom")
	    if(FALSE){
	      if(length(unique(object[[ip]]$v)) > 0.3*nrow(object[[ip]])){
  	      p <- p + ggplot2::scale_fill_manual(values = c("grey70", "white", "red"))
  	    }else{
  	      colors <- c("#5CB85C" ,"#9C9EDE" ,"#FFDC91", "#4DBBD5" ,"#FF9896" ,"#FED439", "#E377C2", "#FED439")
  	      p <- p + ggplot2::scale_color_manual(values = colors)
  	    }## end fi
	    }## end fi
	    ggp[[ip]] <- p
	  }## end for
	  
	  
	}else if(plot.type == "vln"){
  	##*********************##
  	##    Violin plots     ##
  	##*********************##
  	##*
  	
	}else if(plot.type == "prop"){
	  ##*********************##
	  ##    Weight plots     ##
	  ##*********************##
	  ##*
	  if(length(colors) < 3){
	    colors <- c("lightblue","lightyellow","red")
	  }## end fi
	  
	  #ggp <- suppressMessages(ggplot2::ggplot(mData, ggplot2::aes(x, y)) + 
	  #                       ggplot2::geom_point(ggplot2::aes(colour = value),size = 3.0) +
	  #                      ggplot2::scale_color_gradientn(colours = colors) + 
	  #                       ggplot2::scale_x_discrete(expand = c(0, 1)) + ggplot2::scale_y_discrete(expand = c(0,1))+ 
	  #                       ggplot2::facet_wrap(~Cell_Type,ncol = NumCols)+ 
	  #                       ggplot2::coord_fixed()+
	  #                       ggplot2::theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
	  #                             #legend.position=c(0.14,0.76),
	  #                             panel.background = ggplot2::element_blank(),
	  #                             plot.background = ggplot2::element_blank(),
	  #                             panel.border = ggplot2::element_rect(colour = "grey89", fill=NA, size=0.5),
	  #                             axis.text =ggplot2::element_blank(),
	  #                             axis.ticks =ggplot2::element_blank(),
	  #                             axis.title =ggplot2::element_blank(),
	  #                             legend.title=ggplot2::element_text(size = 14,face="bold"),
	  #                             legend.text=ggplot2::element_text(size = 11),
	  #                             strip.text = ggplot2::element_text(size = 12,face="bold"),
	  #                             legend.key = ggplot2::element_rect(colour = "transparent", fill = "white"),
	  #                             legend.key.size = unit(0.45, 'cm')))
	}## end fi## end fi
	
	return(ggp)
}## end function 



#' @rdname SpatialPlot
#' @concept visualization
#' @export
#' @method SpatialPlot Assay
#'
SpatialPlot.Assay <- function(object,
                              slot.use = "data",
                              plot.data = "reduction",
                              plot.type = "dot",
                              location = NULL,
                              features = NULL,
                              dims.plot = 1:2,
                              point.size = 2,
                              text.size = 15,
                              colors = c("#4E84C4", "#FFDB6D"),
                              title = NULL,
                              verbose = TRUE, ...) {
	
	##
	if(is.null(location)){
		location <- GetAssayData(object = object, slot = "location");
		if(!is.matrix(location)){
			location <- as.matrix(location)
		}## end fi
	}## end fi
  
  if((slot.use == "data") && length(GetAssayData(object = object, slot = slot.use)) ==0 ){
    stop("SpatialPlot::The data from provided slot do not exist.")
  }##
  
  ## plot data preparation
  pltdata <- list()
  if(plot.data == "features"){
    if(is.null(features)){stop("SpatialPlot::'features' is required to pass.")}
    data.use <- GetAssayData(object = object, slot = slot.use)
    x <- location[,1]
    y <- location[,2]
    for(ig in features){
      v <- data.use[ig,]
      pltdata <- c(pltdata, list(data.frame(x, y, v)))
    }## end for
    names(pltdata) <- features
    plot.type <- "dot"
  }else if(plot.data == "cluster"){
    v <- GetAssayData(object = object, slot = "ident")
    if(length(v) == 0){stop("SpatialPlot::Please run 'SpatialClust()' function first.")}
    x <- location[,1]
    y <- location[,2]
  
    pltdata[[1]] <- data.frame(x, y, v)
    names(pltdata) <- "clusters"
    plot.type <- "dot"
  }else if(plot.data == "reduction"){
    data.use <- GetAssayData(object = object, slot = "dr")
    x <- location[,1]
    y <- location[,2]
    for(id in dims.plot){
      v <- data.use$embeddings[,id]
      pltdata <- c(pltdata, list(data.frame(x, y, v)))
    }## end for 
    plot.type <- "dot"
    names(pltdata) <- paste0("PC", dims.plot)
  }else if(plot.data == "tsne"){
    data.use <- GetAssayData(object = object, slot = "dr")
    tsne <- Rtsne::Rtsne(data.use$embeddings, dims = 3,check_duplicates = FALSE)
    r = (tsne$Y[,1]-min(tsne$Y[,1]))/(max(tsne$Y[,1])-min(tsne$Y[,1]))
    g = (tsne$Y[,2]-min(tsne$Y[,2]))/(max(tsne$Y[,2])-min(tsne$Y[,2]))
    b = (tsne$Y[,3]-min(tsne$Y[,3]))/(max(tsne$Y[,3])-min(tsne$Y[,3]))
    x = location[,1]
    y = location[,2]
    pltdata[[1]] <- data.frame(x,y,r,g,b)
    names(pltdata) <- "tSNE"
    plot.type <- "rgb"
  }else if(plot.data == "umap"){
    data.use <- GetAssayData(object = object, slot = "dr")
    umap <- umap::umap(data.use$embeddings, n_components = 3)
    r = (umap$layout[,1]-min(umap$layout[,1]))/(max(umap$layout[,1])-min(umap$layout[,1]))
    g = (umap$layout[,2]-min(umap$layout[,2]))/(max(umap$layout[,2])-min(umap$layout[,2]))
    b = (umap$layout[,3]-min(umap$layout[,3]))/(max(umap$layout[,3])-min(umap$layout[,3]))
    x = location[,1]
    y = location[,2]
    pltdata[[1]] <- data.frame(x,y,r,g,b)
    names(pltdata) <- "UMAP"
    plot.type <- "rgb"
  }else if(plot.data == "deconvolution"){
    #data.use <- GetAssayData(object = object, slot = "dec")
    #res_CARD = as.data.frame(proportion)
    #res_CARD = res_CARD[,order(colnames(res_CARD))]
    #location = as.data.frame(spatial_location)
    #if(sum(rownames(res_CARD)==rownames(location))!= nrow(res_CARD)){
    #  stop("The rownames of proportion data does not match with the rownames of spatial location data")
    #}
    #ct.select = ct.visualize
    #res_CARD = res_CARD[,ct.select]
    #res_CARD_scale = as.data.frame(apply(res_CARD,2,function(x){
    #  (x - min(x)) / (max(x) - min(x))
    #} ))
    #res_CARD_scale$x = as.numeric(location$x)
    #res_CARD_scale$y = as.numeric(location$y)
    #mData = melt(res_CARD_scale,id.vars = c("x","y"))
    #colnames(mData)[3] <- "Cell_Type"
    #plot.type <- "rgb"
  }## end fi

	
	## run default SpatialPlot
	ggp <- SpatialPlot(object = pltdata,
	                   plot.type = plot.type,
	                   point.size = point.size,
	                   text.size = text.size,
	                   colors = colors,
	                   title = title,
	                   verbose = verbose, ...)
					

	## store the scaled data in the slot
	return(ggp)
}## end func

#' 
#' @rdname SpatialPlot
#' @concept visualization
#' @export
#' @method SpatialPlot SpatialX
#' 
SpatialPlot.SpatialX <- function(object, 
                                 assay = NULL,
                                 slot.use = "data",
                                 plot.data = "reduction",
                                 plot.type = "dot",
                                 location = NULL,
                                 features = NULL,
                                 dims.plot = 1:2,
                                 point.size = 2,
                                 text.size = 15,
                                 colors = c("#4E84C4", "#FFDB6D"),
                                 title = NULL,
                                 verbose = TRUE, ...) {
	
	## assays
	assay <- assay %||% DefaultAssay(object = object)
	assay.data <- GetAssay(object = object, assay = assay)

	ggp <- SpatialPlot(object = assay.data,
	                   slot.use = slot.use,
	                   plot.data = plot.data,
	                   plot.type = plot.type,
	                   location = location,
	                   features = features,
	                   dims.plot = dims.plot,
	                   point.size = point.size,
	                   text.size = text.size,
	                   colors = colors,
	                   title = title,
	                   verbose = verbose, ... )

	## store back the treated data
	return(ggp)
}## end func


#' Plots for SRT
#' @rdname SpatialPlot
#' @export
#'
SpatialPlot <- function(object, ...) {
	UseMethod(generic = "SpatialPlot", object = object)
}## end func


#########################################
#             CODE END                  #
#########################################
