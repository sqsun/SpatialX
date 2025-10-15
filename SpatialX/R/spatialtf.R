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



#' Normalize raw count data using library size factors
#'
#' This function normalizes raw count data by scaling counts per cell using
#' library size factors and then applies a log transformation. This is a
#' standard preprocessing step for single-cell and spatial transcriptomics
#' data to make counts comparable across cells with different sequencing depths.
#'
#' @param data A matrix, data.frame, or sparse matrix containing raw count data
#'            with features (genes) as rows and cells as columns.
#' @param scale.factor The scaling factor used for normalization. Counts are
#'                    scaled by (library size / scale.factor) before log
#'                    transformation. Default is 10000.
#' @param verbose Logical indicating whether to print progress messages and
#'               normalization statistics. Default is TRUE.
#'
#' @return Returns a normalized matrix with the same dimensions as the input
#'         data. The values are log-transformed (log1p) counts normalized by
#'         library size factors.
#'
#' @details
#' This function implements a standard normalization approach for count-based
#' genomics data:
#'
#' \deqn{\text{normalized} = \log\left(1 + \frac{\text{counts} \times \text{scale.factor}}{\text{library size}}\right)}
#'
#' The normalization process:
#' \enumerate{
#'   \item Calculate library size (total counts) for each cell
#'   \item Scale counts by (library size / scale.factor) to account for
#'         sequencing depth differences
#'   \item Apply log(1 + x) transformation to stabilize variance and
#'         make the data more normally distributed
#' }
#'
#' Key features:
#' \itemize{
#'   \item Handles both dense and sparse matrix inputs efficiently
#'   \item Preserves the structure and dimensions of the input data
#'   \item Uses log1p transformation to handle zero counts gracefully
#'   \item Provides informative messages about the normalization process
#' }
#'
#' This normalization method is particularly useful for:
#' \itemize{
#'   \item Making gene expression comparable across cells with different
#'         sequencing depths
#'   \item Preparing data for dimensionality reduction methods like PCA
#'   \item Improving the performance of clustering algorithms
#'   \item Stabilizing variance for downstream statistical tests
#' }
#'
#' @examples
#' \dontrun{# Create example count matrix
#' mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
#' rownames(mat) <- paste0("Gene_", 1:5)
#' colnames(mat) <- paste0("Cell_", 1:5)
#' 
#' # View raw counts
#' print("Raw counts:")
#' print(mat)
#' 
#' # Normalize with default parameters
#' mat_norm <- LogLSFactor(data = mat)
#' print("Normalized counts:")
#' print(mat_norm)
#' 
#' # Normalize with custom scale factor
#' mat_norm_custom <- LogLSFactor(data = mat, scale.factor = 1000, verbose = TRUE)
#' 
#' # Using with sparse matrix
#' library(Matrix)
#' sparse_mat <- Matrix(mat, sparse = TRUE)
#' sparse_norm <- LogLSFactor(data = sparse_mat)
#' 
#' # Compare library sizes before and after normalization
#' lib_sizes_raw <- colSums(mat)
#' lib_sizes_norm <- colSums(expm1(mat_norm))  # Convert back to count scale
#' plot(lib_sizes_raw, lib_sizes_norm, 
#'      xlab = "Raw library size", ylab = "Normalized library size",
#'      main = "Library Size Comparison")
#' abline(0, 1, col = "red")
#'}
#' @seealso
#' \code{\link{NormalizeVST}},
#' \code{\link{log1p}}
#'
#' @importFrom Matrix Matrix
#' @importFrom methods as
#' @export
#' @concept pre-processing
#'
LogLSFactor <- function(data, scale.factor = 1e4, verbose = TRUE) {
  if (is.data.frame(x = data)) {
    data <- as.matrix(x = data)
  }## end fi
  
  if (!inherits(x = data, what = 'dgCMatrix')) {
    data <- as(object = data, Class = "dgCMatrix")
  }## end fi
  
  ## call Rcpp function to normalize
  if (verbose) {
    cat("Performing log-library size factor normalization\n", file = stderr())
  }## end fi
  
  norm.data <- LogLSFactorCpp(data, scale_factor = scale.factor, display_progress = verbose)
  colnames(x = norm.data) <- colnames(x = data)
  rownames(x = norm.data) <- rownames(x = data)
  return(norm.data)
}## end func



#' Normalize count data using deconvolution size factors
#'
#' This function normalizes raw count data using deconvolution size factors
#' as described in Lun et al. (2016), followed by log transformation. This
#' method is particularly effective for single-cell and spatial transcriptomics
#' data as it handles the high proportion of zeros and varying library sizes
#' more robustly than simple library size normalization.
#'
#' @param data A matrix, data.frame, or sparse matrix containing raw count data
#'            with features (genes) as rows and cells as columns.
#' @param scale.factor Precomputed size factors. If NULL (default), size factors
#'                    are calculated using the deconvolution method from the
#'                    \code{scran} package.
#' @param pseudo.count Pseudo-count to add before log transformation to avoid
#'                    taking the log of zero. Default is 1.
#' @param verbose Logical indicating whether to print progress messages and
#'               normalization statistics. Default is TRUE.
#'
#' @return Returns a normalized matrix with the same dimensions as the input
#'         data. The values are log-transformed (log1p) counts normalized by
#'         deconvolution size factors.
#'
#' @details
#' This function implements the deconvolution normalization method described
#' in Lun et al. (2016), which is particularly well-suited for single-cell
#' and spatial transcriptomics data. The method:
#'
#' \itemize{
#'   \item Groups cells into clusters to pool information across cells
#'   \item Calculates pool-based size factors that are more robust to
#'         differential expression between cell types
#'   \item Deconvolves the pool-based factors to obtain cell-specific
#'         size factors
#'   \item Applies log(1 + x) transformation after normalization
#' }
#'
#' The normalization process:
#' \enumerate{
#'   \item Cluster cells using \code{scran::quickCluster}
#'   \item Calculate deconvolution size factors using \code{scran::calculateSumFactors}
#'   \item Normalize counts by dividing by size factors and multiplying by
#'         the mean size factor
#'   \item Apply log(1 + pseudo.count + normalized counts) transformation
#' }
#'
#' Advantages of deconvolution normalization:
#' \itemize{
#'   \item More robust to heterogeneous cell populations
#'   \item Better handles genes expressed in only a subset of cells
#'   \item Reduces technical bias in downstream analyses
#'   \item Particularly effective for data with many zero counts
#' }
#'
#' @examples
#' \dontrun{
#' # Create example count matrix
#' mat <- matrix(data = rbinom(n = 100, size = 5, prob = 0.2), nrow = 10)
#' rownames(mat) <- paste0("Gene_", 1:10)
#' colnames(mat) <- paste0("Cell_", 1:10)
#'
#' # View raw counts
#' print("Raw counts:")
#' print(mat[1:5, 1:5])
#'
#' # Normalize with deconvolution size factors
#' mat_norm <- LogSizeFactor(data = mat)
#' print("Normalized counts:")
#' print(mat_norm[1:5, 1:5])
#'
#' # Normalize with custom pseudo-count
#' mat_norm_pseudo <- LogSizeFactor(data = mat, pseudo.count = 0.5, verbose = TRUE)
#'
#' # Use precomputed size factors
#' library(scran)
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat))
#' clust <- scran::quickCluster(sce)
#' precomputed_factors <- scran::calculateSumFactors(sce, cluster = clust)
#'
#' mat_norm_precomputed <- LogSizeFactor(
#'   data = mat,
#'   scale.factor = precomputed_factors,
#'   verbose = TRUE
#' )
#'
#' # Compare different normalization methods
#' libsize_norm <- LogLSFactor(mat)
#' deconv_norm <- LogSizeFactor(mat)
#'
#' # Compare distributions
#' par(mfrow = c(1, 2))
#' hist(libsize_norm, main = "Library Size Normalization", breaks = 20)
#' hist(deconv_norm, main = "Deconvolution Normalization", breaks = 20)
#' }
#'
#' @references
#' Lun, A. T., Bach, K., & Marioni, J. C. (2016). Pooling across cells to
#' normalize single-cell RNA sequencing data with many zero counts. Genome
#' Biology, 17(1), 75.
#'
#' @seealso
#' \code{\link{LogLSFactor}}, \code{\link{NormalizeVST}},
#' \code{\link[scran]{calculateSumFactors}}, \code{\link[scran]{quickCluster}}
#'
#' @importFrom Matrix Matrix
#' @importFrom methods as
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom scran quickCluster calculateSumFactors
#' @export
#' @concept preprocessing
#' 
LogSizeFactor <- function(data, 
                          scale.factor = NULL, 
                          pseudo.count = 1, 
                          verbose = TRUE) {
  if (is.data.frame(x = data)) {
    data <- as.matrix(x = data)
  }## end fi
  
  if (!inherits(x = data, what = 'dgCMatrix')) {
    data <- as(object = data, Class = "dgCMatrix")
  }## end fi
  
  ## call Rcpp function to normalize
  if (verbose) {
    cat("Performing log-deconvlution size factor normalization\n", file = stderr())
  }## end fi
  
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=data))
  clust <- scran::quickCluster(sce) 
  if(is.null(scale.factor)){
    scale.factor <- scran::calculateSumFactors(sce, cluster=clust)
  }## end fi
  
  
  norm.data <- LogSizeFactorCpp(data, size_factor = scale.factor, pseudo_count = pseudo.count, display_progress = verbose)
  colnames(x = norm.data) <- colnames(x = data)
  rownames(x = norm.data) <- rownames(x = data)
  return(norm.data)
}## end func

#'
#' Normalize Spatial Transcriptomics Data
#'
#' This function performs various normalization methods on spatial transcriptomics data
#' to account for technical variations and make expression values comparable across
#' cells/spots.
#'
#' @param object A matrix-like object containing raw count data where rows represent
#' features (genes) and columns represent cells/spots.
#' @param norm.method Method for normalization. Available options:
#' \itemize{
#'   \item{"log"}{ - Feature counts for each cell are divided by the total counts for 
#'   that cell and multiplied by the scale.factor. This is then natural-log transformed 
#'   using log1p.}
#'   \item{"logsf"}{ - Log normalization with size factors and pseudo count.}
#' }
#' @param scale.factor Sets the scale factor for cell-level normalization (default: 1e4).
#' For "log" method, this becomes the target count per cell after normalization.
#' @param margin If performing CLR normalization, normalize across features (1) or cells (2).
#' @param block.size How many cells should be run in each chunk, will try to split evenly 
#' across threads. If NULL, automatically determined based on number of workers.
#' @param pseudo.count Pseudo count to add before log transformation for "logsf" method 
#' (default: 1).
#' @param verbose Logical indicating whether to display progress bar for normalization 
#' procedure (default: TRUE).
#' @param ... Additional arguments passed to normalization methods.
#'
#' @return A normalized matrix of the same dimensions as the input object, with 
#' normalized expression values.
#'
#' @details
#' The function provides two main normalization approaches:
#' \itemize{
#'   \item{\strong{log method}:} Performs library size normalization followed by 
#'   log transformation. Each cell's counts are divided by its total counts and 
#'   scaled by \code{scale.factor}, then log1p transformed.
#'   \item{\strong{logsf method}:} Performs log normalization with size factors 
#'   and pseudo count for stabilization.
#' }
#' 
#' For parallel processing, the function automatically detects the number of available
#' workers and splits the data into chunks for efficient computation.
#'
#' @examples
#' \dontrun{# Create example spatial data
#' data <- matrix(rpois(1000, 5), nrow = 100, ncol = 10)
#' rownames(data) <- paste0("Gene", 1:100)
#' colnames(data) <- paste0("Spot", 1:10)
#'
#' # Perform log normalization
#' normalized_data <- NormalizedSpData(object = data, norm.method = "log")
#'
#' # Perform log normalization with custom scale factor
#' normalized_data <- NormalizedSpData(object = data, norm.method = "log", 
#'                                    scale.factor = 10000)
#'
#' # Perform logsf normalization with pseudo count
#' normalized_data <- NormalizedSpData(object = data, norm.method = "logsf", 
#'                                    pseudo.count = 1)
#'}
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_lapply
#' @importFrom methods is
#' @importFrom Matrix rowSums colSums
#' @importFrom rlang %||%
#' 
#' 
#' @rdname NormalizedSpData
#' @concept preprocessing
#' @concept normalization
#' @export
#' 
NormalizedSpData.default <- function(object,
                                  norm.method = "log",
                                  scale.factor = 1e4,
                                  margin = 1,
                                  block.size = NULL,
                                  pseudo.count = 1,
                                  verbose = TRUE, ...) {
	

  if (is.null(x = norm.method)) {
	  return(object)
	}## end fi
	
  normalized.data <- if (future::nbrOfWorkers() > 1) {
		norm.function <- switch(
      EXPR = norm.method,
      'log' = LogLSFactor,
      'logsf' = LogSizeFactor,
      stop("Unknown normalization method: ", norm.method))

		# tryCatch(
		#   expr = Parenting(parent.find = 'SpatialX', margin = margin),
		#   error = function(e) {
		# 	invisible(x = NULL)
		# })
		dsize <- switch(
		  EXPR = margin,
		  '1' = nrow(x = object),
		  '2' = ncol(x = object),
		  stop("'margin' must be 1 or 2")	)
		
		chunk.points <- ChunkPoints(dsize = dsize,
		  csize = block.size %||% ceiling(x = dsize / future::nbrOfWorkers())	)
		
		normalized.data <- future.apply::future_lapply(
			  X = 1:ncol(x = chunk.points),
			  FUN = function(i) {
				block <- chunk.points[, i]
				data <- if (margin == 1) {
				  object[block[1]:block[2], , drop = FALSE]
				} else {
				  object[, block[1]:block[2], drop = FALSE]
				}
				clr_function <- function(x) {
				  return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
				}
				args <- list(
				  data = data,
				  scale.factor = scale.factor,
				  verbose = FALSE,
				  custom_function = clr_function, margin = margin
				)
				args <- args[names(x = formals(fun = norm.function))]
				return(do.call(
				  what = norm.function,
				  args = args
				))
			  }
			)
		do.call(
		  what = switch(
			EXPR = margin,
			'1' = 'rbind',
			'2' = 'cbind',
			stop("'margin' must be 1 or 2")
		  ),
		  args = normalized.data)
	} else {## directly call normalize function
		switch(
		  EXPR = norm.method,
		  'log' = LogLSFactor(data = object, scale.factor = scale.factor, verbose = verbose),
		  'logsf' = LogSizeFactor(data = object, scale.factor = scale.factor, pseudo.count = pseudo.count, verbose = verbose),
		  stop("Unkown normalization method: ", norm.method)	)
	}## end fi
	return(normalized.data)
}## end func


#' @rdname NormalizedSpData
#' @export
#' 
NormalizedSpData <- function(object, ...) {
  UseMethod(generic = "NormalizedSpData", object = object)
}## end func



#' Scale Spatial Transcriptomics Data
#'
#' This function centers and scales the expression data to standardize feature 
#' distributions, making them comparable across cells/spots. It supports both 
#' standard scaling and robust scaling methods, with efficient parallel processing 
#' for large datasets.
#'
#' @param object A matrix-like object containing normalized expression data where 
#' rows represent features (genes) and columns represent cells/spots.
#' @param features Vector of feature names to scale/center. Default is all features 
#' in the object.
#' @param scale.method Method for scaling. Available options:
#' \itemize{
#'   \item{"standard"}{ - Standard scaling (Z-score normalization)}
#'   \item{"robust"}{ - Robust scaling using median and MAD}
#' }
#' @param do.scale Whether to scale the data (default: TRUE).
#' @param do.center Whether to center the data (default: TRUE).
#' @param scale.max Maximum absolute value to return for scaled data. Values beyond 
#' this threshold will be clipped. The default is 10. Setting this can help reduce 
#' the effects of extreme outliers.
#' @param block.size Number of features to process in a single computation block. 
#' Increasing block.size may speed up calculations but at an additional memory cost 
#' (default: 1000).
#' @param min.cells.to.block If object contains fewer than this number of cells, 
#' don't use blocking for scaling calculations (default: 3000).
#' @param verbose Whether to display progress bar for scaling procedure (default: TRUE).
#' @param ... Additional arguments passed to scaling methods.
#'
#' @return A scaled matrix of the same dimensions as the input object, with 
#' centered and/or scaled expression values.
#'
#' @details
#' The function provides two scaling approaches:
#' \itemize{
#'   \item{\strong{standard scaling}:} Centers the data (subtract mean) and scales 
#'   (divide by standard deviation) to produce Z-scores.
#'   \item{\strong{robust scaling}:} Centers using median and scales using median 
#'   absolute deviation (MAD), making it more robust to outliers.
#' }
#' 
#' For large datasets, the function automatically uses parallel processing when 
#' multiple workers are available. Data is processed in blocks to manage memory 
#' usage efficiently.
#'
#' @section Parallel Processing:
#' When multiple workers are available (detected via \code{future::nbrOfWorkers()}),
#' the function processes data in parallel using \code{future.apply::future_lapply}.
#' The data is split into blocks of features for efficient computation.
#'
#' @section Memory Management:
#' For sparse matrices (dgCMatrix, dgTMatrix), the function uses optimized sparse 
#' matrix operations. For dense matrices, it converts to standard matrix format 
#' and uses efficient row-wise operations.
#'
#' @examples
#' \dontrun{# Create example normalized spatial data
#' data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(data) <- paste0("Gene", 1:100)
#' colnames(data) <- paste0("Spot", 1:10)
#'
#' # Perform standard scaling (center and scale)
#' scaled_data <- ScaledSpData(object = data, 
#'                            scale.method = "standard",
#'                            do.scale = TRUE,
#'                            do.center = TRUE)
#'
#' # Perform only centering (no scaling)
#' centered_data <- ScaledSpData(object = data,
#'                              do.scale = FALSE,
#'                              do.center = TRUE)
#'
#' # Perform robust scaling
#' robust_scaled_data <- ScaledSpData(object = data,
#'                                   scale.method = "robust")
#'
#' # Scale only specific features
#' variable_genes <- paste0("Gene", 1:50)
#' subset_scaled <- ScaledSpData(object = data,
#'                              features = variable_genes)
#'}
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_lapply
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats na.omit
#' @rdname ScaledSpData
#' @concept preprocessing
#' @concept scaling
#' @export
#' 
ScaledSpData.default <- function(object,
                              features = NULL,
                              scale.method = "standard",
                              do.scale = TRUE,
                              do.center = TRUE,
                              scale.max = 10,
                              block.size = 1000,
                              min.cells.to.block = 3000,
                              verbose = TRUE, ... ) {
  ## which features to be normalized
  features <- features %||% rownames(x = object)
  features <- as.vector(x = intersect(x = features, y = rownames(x = object)))
  
  object <- object[features, , drop = FALSE]
  object.names <- dimnames(x = object)
  min.cells.to.block <- min(min.cells.to.block, ncol(x = object))
  
  #split.by <- split.by %||% TRUE
  split.cells <- split(x = colnames(x = object), f = TRUE)
  
  #############
  if (verbose && (do.scale || do.center)) {
    msg <- paste(
      na.omit(object = c(
        ifelse(test = do.center, yes = 'centering', no = NA_character_),
        ifelse(test = do.scale, yes = 'scaling', no = NA_character_)
      )),
      collapse = ' and '
    )
    msg <- paste0(
      toupper(x = substr(x = msg, start = 1, stop = 1)),
      substr(x = msg, start = 2, stop = nchar(x = msg)),
      ' data matrix'
    )
    message(msg)
  }## end fi
  
  if (inherits(x = object, what = c('dgCMatrix', 'dgTMatrix'))) {
    scale.function <- FastSparseRowScaleData
  } else {
    object <- as.matrix(x = object)
    scale.function <- FastRowScaleData
  }## end fi
  
  ## nbrOfWorkers
  if(future::nbrOfWorkers() > 1) {
    blocks <- ChunkPoints(dsize = length(x = features), csize = block.size)
    chunks <- expand.grid(
      names(x = split.cells),
      1:ncol(x = blocks),
      stringsAsFactors = FALSE)
    
    scaled.data <- future.apply::future_lapply(X = 1:nrow(x = chunks),
                                 FUN = function(index) {
                                   row <- chunks[index, ]
                                   group <- row[[1]]
                                   block <- as.vector(x = blocks[, as.numeric(x = row[[2]])])
                                   arg.list <- list(
                                     mat = object[features[block[1]:block[2]], split.cells[[group]], drop = FALSE],
                                     scale = do.scale,
                                     center = do.center,
                                     scale_max = scale.max,
                                     display_progress = FALSE )
        ## 
        arg.list <- arg.list[intersect(x = names(x = arg.list), y = names(x = formals(fun = scale.function)))]
        data.scale <- do.call(what = scale.function, args = arg.list)
        dimnames(x = data.scale) <- dimnames(x = object[features[block[1]:block[2]], split.cells[[group]]])
        suppressWarnings(expr = data.scale[is.na(x = data.scale)] <- 0)
        return(data.scale)
      })
    
    if (length(x = split.cells) > 1) {
      merge.indices <- lapply(X = 1:length(x = split.cells),
        FUN = seq.int,
        to = length(x = scaled.data),
        by = length(x = split.cells) )
      scaled.data <- lapply(X = merge.indices,
        FUN = function(x) {
          return(suppressWarnings(expr = do.call(what = 'rbind', args = scaled.data[x])))
        })
      
      scaled.data <- suppressWarnings(expr = do.call(what = 'cbind', args = scaled.data))
    } else {
      suppressWarnings(expr = scaled.data <- do.call(what = 'rbind', args = scaled.data))
    }
  } else {
    ## initialize data
    scaled.data <- matrix(data = NA_real_, nrow = nrow(x = object), ncol = ncol(x = object), dimnames = object.names )
    
    max.block <- ceiling(x = length(x = features) / block.size)
    for (x in names(x = split.cells)) {
      if (verbose) {
        if (length(x = split.cells) > 1 && (do.scale || do.center)) {
          message(gsub(pattern = 'matrix', replacement = 'from split ', x = msg), x)
        }
        pb <- txtProgressBar(min = 0, max = max.block, style = 3, file = stderr())
      }## end fi
      
      for (i in 1:max.block) {
        my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
        my.inds <- my.inds[my.inds <= length(x = features)]
        arg.list <- list(
          mat = object[features[my.inds], split.cells[[x]], drop = FALSE],
          scale = do.scale,
          center = do.center,
          scale_max = scale.max,
          display_progress = FALSE)
        
        arg.list <- arg.list[intersect(x = names(x = arg.list), y = names(x = formals(fun = scale.function)))]
        data.scale <- do.call(what = scale.function, args = arg.list)
        dimnames(x = data.scale) <- dimnames(x = object[features[my.inds], split.cells[[x]]])
        scaled.data[features[my.inds], split.cells[[x]]] <- data.scale
        rm(data.scale)
     
        if(verbose){
          setTxtProgressBar(pb = pb, value = i)
        }## end fi
      }## end for
      
      if(verbose){
        close(con = pb)
      }## end fi
      
    }## end for
  }## end fi
  
  dimnames(x = scaled.data) <- object.names
  scaled.data[is.na(x = scaled.data)] <- 0

  return(scaled.data)
}## end func


#' 
#' @rdname ScaledSpData
#' @export
#'
ScaledSpData <- function(object, ...) {
  UseMethod(generic = "ScaledSpData", object = object)
}## end func


#' Spatial Transcriptomics Data Transformation
#'
#' This function performs various transformations on spatial transcriptomics data
#' including normalization, variance stabilization, and scaling. It supports multiple
#' transformation methods and can handle both Assay and SpatialX objects.
#'
#' @param object An object of class \code{Assay} or \code{SpatialX} containing
#' spatial transcriptomics data.
#' @param features Vector of feature names to transform. Default is all features
#' in the object.
#' @param tf.method Transformation method to apply. Available options:
#' \itemize{
#'   \item{"log"}{ - Log normalization with library size factors}
#'   \item{"logsf"}{ - Log normalization with size factors and pseudo count}
#'   \item{"vst"}{ - Variance Stabilizing Transformation using sctransform}
#'   \item{"standard"}{ - Standard scaling (Z-score normalization)}
#' }
#' @param slot.use Name of the assay slot to use as input (default: "counts").
#' @param do.scale Whether to scale the data when using "standard" method (default: TRUE).
#' @param do.center Whether to center the data when using "standard" method (default: TRUE).
#' @param scale.factor Scale factor for normalization methods (default: 1e4).
#' @param pseudo.count Pseudo count to add before log transformation for "logsf"
#' method (default: 1).
#' @param margin Margin for normalization (1 for features, 2 for cells) (default: 1).
#' @param scale.max Maximum absolute value for scaled data (default: 10).
#' @param block.size Number of features to process in a single computation block
#' (default: 1000).
#' @param min.cells.to.block Minimum number of cells to use blocking (default: 3000).
#' @param verbose Whether to display progress information (default: TRUE).
#' @param assay Name of assay to use for SpatialX objects.
#' @param ... Additional arguments passed to transformation methods.
#'
#' @return For \code{Assay} objects: Returns a modified Assay object with transformed
#' data stored in the appropriate slot. For \code{SpatialX} objects: Returns a
#' modified SpatialX object with the specified assay transformed.
#'
#' @details
#' This function provides a unified interface for multiple data transformation
#' techniques commonly used in spatial transcriptomics analysis:
#'
#' \itemize{
#'   \item{\strong{log}:} Performs log normalization with library size factors.
#'   Transformed data is stored in the 'data' slot.
#'   \item{\strong{logsf}:} Performs log normalization with size factors and
#'   pseudo count for stabilization. Transformed data is stored in the 'data' slot.
#'   \item{\strong{vst}:} Applies Variance Stabilizing Transformation using the
#'   \code{sctransform} package. This method models the technical noise and
#'   returns residuals that are approximately normally distributed. Transformed
#'   data is stored in the 'scale.data' slot.
#'   \item{\strong{standard}:} Applies standard scaling (centering and scaling
#'   to Z-scores). Transformed data is stored in the 'data' slot.
#' }
#'
#' @section Slot Usage:
#' Different transformation methods store results in different slots:
#' \itemize{
#'   \item{log, logsf, standard:} Store results in the 'data' slot
#'   \item{vst:} Stores results in the 'scale.data' slot
#' }
#'
#' @section Method-Specific Details:
#' \describe{
#'   \item{VST Transformation}{When using the 'vst' method, the function calls
#'   \code{sctransform::vst} with latent variables including "log_umi". This
#'   method is particularly effective for removing technical variability while
#'   preserving biological signal.}
#'   \item{Log Transformations}{Both 'log' and 'logsf' methods use the
#'   \code{NormalizedSpData} function internally, with different normalization
#'   strategies.}
#'   \item{Standard Scaling}{Uses the \code{ScaledSpData} function with the
#'   "standard" method for Z-score normalization.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example with SpatialX object
#' data("pbmc_small")
#' 
#' # Apply log transformation
#' object <- SpatialTF(object = object, tf.method = "log")
#' 
#' # Apply variance stabilizing transformation
#' object <- SpatialTF(object = object, tf.method = "vst")
#' 
#' # Apply standard scaling to specific features
#' variable_features <- c("Gene1", "Gene2", "Gene3")
#' object <- SpatialTF(object = object, 
#'                    tf.method = "standard",
#'                    features = variable_features)
#' 
#' # Apply log transformation with custom scale factor
#' object <- SpatialTF(object = object,
#'                    tf.method = "log",
#'                    scale.factor = 10000)
#' }
#'
#' @importFrom sctransform vst
#' @rdname SpatialTF
#' @concept preprocessing
#' @concept transformation
#' @export
#' @method SpatialTF Assay
#'
SpatialTF.Assay <- function(object,
                            features = NULL,
                            tf.method = c("log","logsf","vst","standard"),
                            slot.use = "counts",
                            do.scale = TRUE,
                            do.center = TRUE,
                            scale.factor = 1e4,
                            pseudo.count = 1,
                            margin = 1,
                            scale.max = 10,
                            block.size = 1000,
                            min.cells.to.block = 3000,
                            verbose = TRUE, ... ) {
  
  ## only match predefined arg
  tf.method <- match.arg(tf.method)
	##
  features <- features %||% rownames(x = GetAssayData(object = object, slot = slot.use))
	##
  switch(EXPR = tolower(tf.method),
    'log' = {new.data <- NormalizedSpData(object = GetAssayData(object = object, slot = slot.use),
                                          norm.method = "log",
                                          scale.factor = scale.factor,
                                          pseudo.count = pseudo.count,
                                          verbose = verbose,
                                          margin = margin, ...)
              ## store the scaled data in the slot
              object <- SetAssayData(object = object, slot = 'data', new.data = new.data)},
    'vst' = {## Variance Stabilizing Transformation
            vst_out <- sctransform::vst(umi = GetAssayData(object = object, slot = slot.use), 
                                  latent_var = c("log_umi"), 
                                  return_gene_attr = TRUE,
                                  return_cell_attr = TRUE, 
                                  verbosity = 1)
            new.data <- vst_out$y;
            ## store the scaled data in the slot
            object <- SetAssayData(object = object, slot = 'scale.data', new.data = new.data)
      },
    'standard' = {new.data <- ScaledSpData(object = GetAssayData(object = object, slot = "data"), 
                                           features = features, 
                                           scale.method = "standard", 
                                           do.scale = do.scale, 
                                           do.center = do.center, 
                                           scale.max = scale.max, 
                                           block.size = block.size, 
                                           min.cells.to.block = min.cells.to.block, 
                                           verbose = verbose, ...)
                ## store the scaled data in the slot
                object <- SetAssayData(object = object, slot = 'data', new.data = new.data)},
    'logsf' = {new.data <- NormalizedSpData(object = GetAssayData(object = object, slot = slot.use),
                                            norm.method = "logsf",
                                            scale.factor = scale.factor,
                                            pseudo.count = pseudo.count,
                                            verbose = verbose,
                                            margin = margin, ...)
              ## store the scaled data in the slot
              object <- SetAssayData(object = object, slot = 'data', new.data = new.data)},
    stop("'tf.method' must be 'log', 'logsf', 'standard', or 'logsf'.")	)
  
	## return the results
	return(object)
}## end func

#' SRT data transformation for SpatialX object
#' @rdname SpatialTF
#' @concept preprocessing
#' @export
#' @method SpatialTF SpatialX
#'
SpatialTF.SpatialX <- function(object,
                            assay = NULL,
                            features = NULL,
                            tf.method = "log",
                            slot.use = "counts",
                            do.scale = TRUE,
                            do.center = TRUE,
                            scale.factor = 1e4,
                            pseudo.count = 1,
                            margin = 1,
                            scale.max = 10,
                            block.size = 1000,
                            min.cells.to.block = 3000,
                            verbose = TRUE, ...) {
  
	## create a folder to store the results if does not exist
	assay <- assay %||% DefaultAssay(object = object)
	assay.data <- GetAssay(object = object, assay = assay)
  
	##
	features <- features %||% rownames(x = object)
	
	## normalize data using different methods
    assay.data <- SpatialTF(object = assay.data,
                            tf.method = tf.method,
                            scale.factor = scale.factor,
                            features = features,
                            pseduo.count = pseudo.count,
                            margin = margin,
                            do.scale = do.scale,
                            do.center = do.center,
                            scale.max = scale.max,
                            block.size = block.size,
                            min.cells.to.block = min.cells.to.block,
                            verbose = verbose,...)

  object@assays[[assay]] <- assay.data
  return(object)
}## end func


#' Transformation for SRT data
#' @rdname SpatialTF
#' @export
#' 
SpatialTF <- function(object, ...) {
	UseMethod(generic = "SpatialTF", object = object)
}## end func



#########################################
#             CODE END                  #
#########################################
