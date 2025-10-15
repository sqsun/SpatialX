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



#' Compute single kernel matrix
#'
#' This function computes various types of kernel matrices for spatial analysis,
#' machine learning, and statistical modeling. It supports multiple kernel types
#' including Gaussian (RBF), polynomial, linear, CAR (Conditional Autoregressive),
#' cosine, Matern, and spline kernels, with options for sparse matrix computation
#' and parallel processing.
#'
#' @param x1 Primary data matrix for kernel computation. Rows represent samples
#'          and columns represent features. Required unless ED is provided.
#' @param x2 Optional secondary data matrix for cross-kernel computation.
#'          If NA (default), computes kernel within x1 (self-similarity).
#' @param ED Precomputed Euclidean distance matrix. If provided, skips distance
#'          computation and uses this matrix directly.
#' @param kernel.type Type of kernel to compute. Options: "gaussian" (RBF),
#'          "poly" (polynomial), "linear", "car" (Conditional Autoregressive),
#'          "cosine", "matern", or "spline". Default is "gaussian".
#' @param kernel.params Parameters for the kernel function. For Gaussian kernel,
#'          this is the bandwidth parameter; for polynomial, the degree; etc.
#'          If NA, automatically computed using \code{ComputeKernelParamLessMem}.
#' @param sparse Logical indicating whether to compute sparse kernel matrix.
#'          Useful for large datasets to save memory. Default is FALSE.
#' @param num.core Number of CPU cores to use for parallel computation when
#'          \code{sparse = TRUE}. Default is 5.
#' @param check.positive Logical indicating whether to check and ensure the
#'          kernel matrix is positive definite. Default is FALSE.
#' @param tol Tolerance threshold for sparse matrix computation. Elements below
#'          this value are treated as zero. Default is 1e-3.
#'
#' @return Returns a kernel matrix of the specified type. The matrix can be
#'         dense or sparse depending on the \code{sparse} parameter. For
#'         self-similarity kernels, returns a square symmetric matrix.
#'
#' @details
#' This function provides a unified interface for computing various kernel
#' matrices commonly used in spatial statistics, machine learning, and
#' genomic analysis. Key features:
#'
#' \itemize{
#'   \item \strong{Multiple Kernel Types}:
#'     \itemize{
#'       \item \code{gaussian}: Radial Basis Function (RBF) kernel
#'       \item \code{poly}: Polynomial kernel
#'       \item \code{linear}: Linear kernel (dot product)
#'       \item \code{car}: Conditional Autoregressive kernel for spatial data
#'       \item \code{cosine}: Cosine similarity kernel
#'       \item \code{matern}: Matern covariance kernel for spatial data
#'       \item \code{spline}: Thin plate spline kernel
#'     }
#'   \item \strong{Automatic Parameter Tuning}: When \code{kernel.params} is NA,
#'         automatically computes appropriate parameters based on the data.
#'   \item \strong{Sparse Computation}: Option to compute sparse kernels for
#'         memory efficiency with large datasets.
#'   \item \strong{Parallel Processing}: Supports parallel computation for
#'         sparse kernel matrices.
#'   \item \strong{Data Preprocessing}: Automatically scales and centers input data.
#' }
#'
#' The function automatically handles data preprocessing by scaling and centering
#' the input matrices. For spatial kernels (CAR, Matern, spline), the function
#' assumes the input represents spatial coordinates.
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' data_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' spatial_coords <- matrix(runif(20), ncol = 2)
#'
#' # Compute Gaussian (RBF) kernel
#' gaussian_kernel <- ComputeSingleKernelMat(
#'   x1 = data_matrix,
#'   kernel.type = "gaussian",
#'   kernel.params = 1.0
#' )
#'
#' # Compute sparse Gaussian kernel with parallel processing
#' sparse_gaussian <- ComputeSingleKernelMat(
#'   x1 = data_matrix,
#'   kernel.type = "gaussian",
#'   sparse = TRUE,
#'   num.core = 4,
#'   tol = 1e-4
#' )
#'
#' # Compute polynomial kernel
#' poly_kernel <- ComputeSingleKernelMat(
#'   x1 = data_matrix,
#'   kernel.type = "poly",
#'   kernel.params = 3
#' )
#'
#' # Compute CAR kernel for spatial data
#' car_kernel <- ComputeSingleKernelMat(
#'   x1 = spatial_coords,
#'   kernel.type = "car",
#'   kernel.params = 0.5
#' )
#'
#' # Compute Matern kernel
#' matern_kernel <- ComputeSingleKernelMat(
#'   x1 = spatial_coords,
#'   kernel.type = "matern",
#'   kernel.params = 1.5
#' )
#'
#' # Use precomputed distance matrix
#' dist_matrix <- as.matrix(dist(data_matrix))
#' kernel_from_dist <- ComputeSingleKernelMat(
#'   ED = dist_matrix,
#'   kernel.type = "gaussian",
#'   kernel.params = 1.0
#' )
#' }
#'
#' @seealso
#' \code{\link{ComputeKernelParamLessMem}}, \code{\link{SysMatEigen}},
#' \code{\link{ComputeDmW}}, \code{\link{MatrixInverse}}
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom pdist pdist
#' @importFrom parallel mclapply
#' @importFrom tidyr unnest_wider
#' @importFrom fields Matern
#' @importFrom dplyr %>%
#' @importFrom mgcv smoothCon
#' @importFrom stats dist sd
#' @export
#' 
ComputeSingleKernelMat <- function(x1 = NULL, 
                                   x2 = NA, 
                                   ED = NULL, 
                                   kernel.type = NULL, 
                                   kernel.params = NA, 
                                   sparse = FALSE, 
                                   num.core = 5, 
                                   check.positive = FALSE,
                                   tol = 1e-3){
    
    ## set the parameters for x1
    if(is.null(ED)){
		n1 = nrow(x1)
		d1 = ncol(x1)
		
		##normalize the data 
		#x1 <- apply(x1, MARGIN = 1, FUN = function(X) (X - mean(X)) / sd(X)) 
		#x1 <- t(x1)
		x1 <- scale(x1)
		## set the parameters for x2
		if(any(is.na(x2))){
			n2 = n1
			d2 = d1
			flag = 0
		}else{
			n2 = nrow(x2)
			d2 = ncol(x2)
			x2 <- apply(x2, MARGIN = 1, FUN = function(X) (X - mean(X)) / sd(X)) 
			if(d1!=d2) {stop("Error in the data.")}
			flag = 1
		}## end fi
	}## end fi null ED
	
	if(any(is.na(kernel.params))){
		# kernel.params <- 2.0
	  kernel.params <- ComputeKernelParamLessMem(x1)$kparam[5]
	}## end fi
    
	if(is.null(kernel.type)){kernel.type <- "gaussian"}
   
	## consider any kernel type
    kernel.type <- tolower(kernel.type)
	if(kernel.type == "linear"){
        kernel_mat <- x1 %*% t(x2)
    }else if(kernel.type == "poly"){
        kernel_mat <- (x1 %*% t(x2))^kernel.params
	}else if(kernel.type == "gaussian"){
        #ED <- EucDist2(x1,x2)
		    if(sparse){
		      fx_gaussian <- function(i){
		        line_i <- rep(0, dim(x1)[1])
		        line_i[i] <- 1
		        line_i[-i] <- exp(-(pdist::pdist(x1[i, ], x1[-i, ])@dist^2)/(kernel.params))
		        ind_i <- which(line_i >= tol)
		        return(list("ind_i" = ind_i, "ind_j" = rep(i, length(ind_i)), "val_i" = line_i[ind_i]))
		      }
		      ## Aggregate the sparse matrix
		      results <- parallel::mclapply(1:dim(x1)[1], fx_gaussian, mc.cores = num.core)
		      tib <- tidyr::tibble(results) %>% tidyr::unnest_wider(results)
		      kernel_mat <- Matrix::sparseMatrix(i = unlist(tib[[1]]), j = unlist(tib[[2]]), 
		                                       x = unlist(tib[[3]]), dims = c(dim(x1)[1],
		                                                                      dim(x1)[1]))
		    }else{
		      if(is.null(ED)){ED <- dist(x1)}
		      kernel_mat <- exp(-ED/(2*kernel.params*kernel.params) )
		      kernel_mat <- as.matrix(kernel_mat)
		    }
	  }else if(kernel.type == "cosine"){
        #ED <- EucDist2(x1,x2)
		    if(sparse){
		      fx_cosine <- function(i){
		        line_i <- rep(0, dim(x1)[1])
		        line_i[i] <- 1
		        line_i[-i] <- cos(2*pi*(pdist::pdist(x1[i, ], x1[-i, ])@dist)/(kernel.params))
		        ind_i <- which(line_i >= tol)
		        return(list("ind_i" = ind_i, "ind_j" = rep(i, length(ind_i)), "val_i" = line_i[ind_i]))
		      }
		      ## Aggregate the sparse matrix
		      results <- parallel::mclapply(1:dim(x1)[1], fx_cosine, mc.cores = num.core)
		      tib <- tidyr::tibble(results) %>% tidyr::unnest_wider(results)
		      kernel_mat <- Matrix::sparseMatrix(i = unlist(tib[[1]]), j = unlist(tib[[2]]), 
		                                         x = unlist(tib[[3]]), dims = c(dim(x1)[1],
		                                                                        dim(x1)[1]))
		    }else{
		      if(is.null(ED)){ED <- dist(x1)}
		      #kernel_mat <- exp(-sqrt(Dist)/(2*kernel.params^2))
		      kernel_mat <- cos(2*pi*ED/kernel.params)
		      kernel_mat <- as.matrix(kernel_mat)
		    }
	  }else if(kernel.type == "car"){
        #ED <- EucDist2(x1,x2)
    		if(is.null(ED)){ED <- as.matrix(dist(x1))}
    		kernel2 <- ComputeDmW(ED, kernel.params, 1.0, 0.5);
    		kernel_mat <- MatrixInverse(kernel2)$inv_mat
	  }else if(kernel.type == "matern"){## edited by sun, 2025-04-11
    	  if(is.null(ED)){ED <- as.matrix(dist(x1))}
    	  kernel_mat <- fields::Matern(d = ED, range = 1, alpha = kernel.params, smoothness = 1.5);
    	  # kernel_mat <- MatrixInverse(kernel_mat)$inv_mat
	  }else if(kernel.type == "spline"){## edited by sun, 2025-04-11
    	  center_coords <- sweep(x1, 2, apply(x1, 2, mean), '-')
    	  center_coords <- as.data.frame(center_coords / sd(as.matrix(center_coords)))
    	  x <- y <- NULL
    	  colnames(center_coords) <- c("x", "y")
    	  sm <- mgcv::smoothCon(mgcv::s(x, y, k = 4, fx = T, bs = 'tp'), data = center_coords)[[1]]
    	  mm <- as.matrix(data.frame(sm$X))
    	  mm_inv <- solve(crossprod(mm, mm))
    	  kernel_mat <- list(mm %*% mm_inv %*% t(mm))
    	  # kernel_mat <- MatrixInverse(kernel_mat)$inv_mat
	  }else{
	      stop("ComputeSingleKernel::Provide a propriate kernel type.")
	  }## end fi
	
	## kernel matrix should be positive definition matrix
	if(check.positive){
		cat(paste("#### chekcing kernel matrix that should be positive definition \n"))
		kernel_mat <- SysMatEigen(kernel_mat)$kernel_mat
		gc()
	}## end fi check
    return(kernel_mat)
}## end function



#' Combine multiple p-values using Cauchy combination rule
#'
#' This function combines p-values from multiple tests or kernels using the
#' Cauchy Combination Test (also known as ACAT). It is designed to efficiently
#' handle matrix inputs where each row represents a feature (e.g., gene) and
#' each column represents a different test or kernel.
#'
#' @param pvalues A matrix of p-values where rows represent features (e.g., genes)
#'               and columns represent different tests or kernels. Typically a
#'               p x k matrix where k is the number of kernels (e.g., 10).
#' @param weights An optional matrix of weights with the same dimensions as
#'               \code{pvalues}. If NULL, equal weights are used for all tests.
#'               Weights should be non-negative and will be normalized for each
#'               feature (row-wise normalization).
#'
#' @return A numeric vector of combined p-values, one for each row (feature)
#'         in the input matrix. The vector has the same length as the number
#'         of rows in \code{pvalues}.
#'
#' @details
#' The Cauchy Combination Test (ACAT) is a powerful method for combining
#' dependent p-values from multiple tests. This implementation is optimized
#' for the common use case in spatial transcriptomics where multiple spatial
#' kernels are tested for each gene.
#'
#' For each feature (row), the function:
#' \enumerate{
#'   \item Normalizes weights to sum to 1 for that feature
#'   \item Applies the Cauchy transformation: \eqn{w_i \cdot \tan((0.5 - p_i) \cdot \pi)}
#'   \item Sums the transformed values across tests/kernels
#'   \item Converts back to a p-value using the Cauchy CDF
#' }
#'
#' Key advantages of this method:
#' \itemize{
#'   \item Robust to dependencies between tests
#'   \item Powerful under sparse alternatives
#'   \item Computationally efficient for large matrices
#'   \item Handles both weighted and unweighted combination
#' }
#'
#' The function automatically handles edge cases:
#' \itemize{
#'   \item If any p-value in a row is 0, the combined p-value is 0
#'   \item P-values of 1 are adjusted to avoid numerical issues
#'   \item NA values are removed with appropriate warnings
#'   \item Invalid p-values (<0 or >1) trigger errors
#' }
#'
#' @examples
#' \dontrun{
#' # Example with 100 genes and 10 kernels
#' pvals_matrix <- matrix(runif(1000), nrow = 100, ncol = 10)
#' 
#' # Combine with equal weights
#' combined_equal <- CombinePValues(pvals_matrix)
#' 
#' # Create custom weights (e.g., based on kernel importance)
#' weights_matrix <- matrix(rep(1:10, each = 100), nrow = 100, ncol = 10)
#' combined_weighted <- CombinePValues(pvals_matrix, weights_matrix)
#' 
#' # Use in spatial transcriptomics analysis
#' spatial_pvals <- GetSpatialPValues(spatial_object, kernels = 10)
#' significant_genes <- CombinePValues(spatial_pvals)
#' 
#' # Identify significant genes at FDR 0.05
#' fdr_adjusted <- p.adjust(significant_genes, method = "fdr")
#' sig_genes <- which(fdr_adjusted < 0.05)
#' }
#'
#' @references
#' Liu Y., et al. (2019). ACAT: A Fast and Powerful p Value Combination
#' Method for Rare-Variant Analysis in Sequencing Studies.
#' The American Journal of Human Genetics, 104(3), 410-421.
#'
#' @seealso
#' \code{\link{ComputeACAT}}, \code{\link[stats]{p.adjust}}
#'
#' @importFrom stats pcauchy
#' @export
#' 
CombinePValues <- function(pvalues, weights=NULL){
	if(!is.matrix(pvalues)){pvalues <- as.matrix(pvalues)}
	## to avoid extremely values
	pvalues[which(pvalues==0)] <- 5.55e-17
	pvalues[which((1-pvalues)<1e-3)] <- 0.99
	
	num_pval <- ncol(pvalues)
	num_gene <- nrow(pvalues)
	if(is.null(weights)){
		weights <- matrix(rep(1.0/num_pval, num_pval*num_gene), ncol=num_pval )
	}# end fi
	if( (nrow(weights) != num_gene) || (ncol(weights) != num_pval)){
		stop("the dimensions of weights does not match that of combined pvalues")
	}# end fi
	
	Cstat <- tan((0.5 - pvalues)*pi)

	wCstat <- weights*Cstat
	Cbar <- apply(wCstat, 1, sum)
	#combined_pval <- 1.0/2.0 - atan(Cbar)/pi
	combined_pval <- 1.0 - pcauchy(Cbar)	
	combined_pval[which(combined_pval <= 0)] <- 5.55e-17
	return(combined_pval)
}# end func



#' Combine p-values using Aggregated Cauchy Association Test (ACAT)
#'
#' This function implements the ACAT method for combining p-values from multiple
#' tests or kernels. ACAT is particularly effective for combining dependent
#' p-values and is widely used in genetic association studies and spatial
#' analysis where test statistics may be correlated.
#'
#' @param Pvals A numeric vector of p-values to be combined. P-values should be
#'              between 0 and 1. NA values will be removed with a warning.
#' @param Weights An optional numeric vector of weights for each p-value.
#'               If NULL (default), equal weights are used. Weights should be
#'               non-negative and will be normalized to sum to 1.
#'
#' @return A single numeric value representing the combined p-value using the
#'         ACAT method. The returned value is between 0 and 1.
#'
#' @details
#' The Aggregated Cauchy Association Test (ACAT) combines p-values by:
#' \enumerate{
#'   \item Transforming each p-value using the Cauchy transformation:
#'         \eqn{w_i \cdot \tan((0.5 - p_i) \cdot \pi)}
#'   \item Summing the transformed values
#'   \item Converting the sum back to a p-value using the Cauchy cumulative
#'         distribution function
#' }
#'
#' ACAT has several advantages:
#' \itemize{
#'   \item Robust to dependencies between tests
#'   \item Powerful under various alternative hypotheses
#'   \item Computationally efficient
#'   \item Works well with both small and large numbers of tests
#' }
#'
#' The method automatically handles extreme p-values:
#' \itemize{
#'   \item If any p-value is 0, the combined p-value is 0
#'   \item P-values of 1 are adjusted to 1 - ε to avoid numerical issues
#'   \item NA values are removed with a warning
#' }
#'
#' @examples
#' \dontrun{
#' # Combine p-values with equal weights
#' pvals <- c(0.01, 0.05, 0.1, 0.2)
#' combined_p <- ComputeACAT(pvals)
#' 
#' # Combine p-values with custom weights
#' weights <- c(2, 1, 1, 2)
#' combined_p_weighted <- ComputeACAT(pvals, weights)
#' 
#' # Handle extreme values
#' extreme_pvals <- c(0, 0.001, 0.5, 1)
#' combined_extreme <- ComputeACAT(extreme_pvals)
#' 
#' # Use in spatial analysis with multiple kernel p-values
#' kernel_pvals <- c(0.03, 0.07, 0.12, 0.25, 0.08)
#' spatial_combined <- ComputeACAT(kernel_pvals)
#' }
#'
#' @references
#' Liu Y., et al. (2019). ACAT: A Fast and Powerful p Value Combination
#' Method for Rare-Variant Analysis in Sequencing Studies.
#' The American Journal of Human Genetics, 104(3), 410-421.
#'
#' @seealso
#' \code{\link[stats]{p.adjust}} for other p-value adjustment methods
#' @importFrom stats pcauchy
#' 
#' @export
#' 
ComputeACAT <- function(Pvals, Weights=NULL){
	#### check if there is NA
	if(sum(is.na(Pvals)) > 0){
		stop("Cannot have NAs in the p-values!")
	}## end fi
	
	#### check if Pvals are between 0 and 1
	if((sum(Pvals<0) + sum(Pvals>1)) > 0){
		stop("P-values must be between 0 and 1!")
	}## end fi
	
	#### check if there are pvals that are either exactly 0 or 1.
	# is.zero <- (sum(Pvals==0) >= 1)
	# is.one <- (sum(Pvals==1) >= 1)
	Pvals[which(Pvals == 0)] <- 5.55e-17
	Pvals[which((1.0 - Pvals)<1e-3)] <- 0.99
	
	#if(is.zero && is.one){stop("Cannot have both 0 and 1 p-values!")}## end fi
	# if(is.zero && is.one){return(NA)}## end fi
	
	# if(is.zero){return(1e-300)}## end fi

	# if(is.one){
		##warning("There are p-values that are exactly 1!")
		# return(0.9999)
	# }## end fi

	#### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
	if(is.null(Weights)){
		Weights <- rep(1/length(Pvals), length(Pvals))
	}else if (length(Weights) != length(Pvals)){
		stop("The length of weights should be the same as that of the p-values")
	}else if (sum(Weights<0) > 0){
		stop("All the weights must be positive!")
	}else{
		Weights <- Weights/sum(Weights)
	}## end fi


	#### check if there are very small non-zero p values
	is.small <- (Pvals < 1e-17)
	if(sum(is.small) == 0){
		cct.stat <- sum(Weights*tan((0.5 - Pvals)*pi))
	}else{
		cct.stat <- sum((Weights[is.small]/Pvals[is.small])/pi)
		cct.stat <- cct.stat + sum(Weights[!is.small]*tan((0.5 - Pvals[!is.small])*pi))
	}## end fi

	#### check if the test statistic is very large.
	if(cct.stat > 1e+15){
		pval <- (1/cct.stat)/pi
	}else{
		pval <- 1 - pcauchy(cct.stat)
	}## end fi
	return(pval)
}## end func

#' Anscombe variance stabilizing transformation for negative binomial data
#'
#' These functions perform Anscombe's variance stabilizing transformation on
#' count data following a negative binomial distribution. The transformation
#' stabilizes the variance across different expression levels, making the data
#' more suitable for downstream statistical analyses that assume homoscedasticity.
#'
#' @param counts A gene expression count matrix with features (genes) as rows
#'              and samples (cells/spots) as columns. Should contain
#'              non-negative integer counts.
#' @param sv Normalization parameter used as starting value for estimating
#'          the dispersion parameter phi. Default is 1.
#'
#' @return Returns a matrix of variance-stabilized and normalized expression
#'         values with the same dimensions as the input counts matrix. The
#'         values are on a continuous scale with stabilized variance.
#'
#' @details
#' The Anscombe transformation for negative binomial data is defined as:
#' \deqn{\log\left(\frac{counts + \frac{1}{2\phi}}{total}\right)}
#' where \eqn{\phi} is the dispersion parameter estimated from the data.
#'
#' The transformation consists of two main steps:
#' \enumerate{
#'   \item Estimate the dispersion parameter \eqn{\phi} by fitting the
#'         relationship: \eqn{Var(X) = \mu + \phi \mu^2}
#'   \item Apply the variance stabilizing transformation and regress out
#'         library size effects
#' }
#'
#' This approach is particularly useful for:
#' \itemize{
#'   \item Preparing count data for methods assuming homoscedastic variance
#'   \item Reducing the influence of extreme values in highly variable genes
#'   \item Improving performance of clustering and dimensionality reduction
#'   \item Enabling the use of Euclidean distance-based methods
#' }
#'
#' @section Functions:
#' \describe{
#'   \item{\code{NormalizeVST}}{Standard implementation using base R functions
#'         for variance and mean calculation. Suitable for small to medium datasets.}
#'   \item{\code{NormalizeVSTfast}}{Optimized implementation using Rcpp functions
#'         for efficient computation of row means, variances, and sums. Recommended
#'         for large datasets.}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate example count data
#' set.seed(123)
#' counts <- matrix(rnbinom(1000, size = 2, mu = 10), nrow = 100, ncol = 10)
#' 
#' # Apply variance stabilizing transformation
#' vst_data <- NormalizeVST(counts)
#' 
#' # Check variance stabilization
#' original_vars <- apply(counts, 1, var)
#' transformed_vars <- apply(vst_data, 1, var)
#' 
#' # Compare variance distributions
#' summary(original_vars)
#' summary(transformed_vars)
#' 
#' # Use fast version for large datasets
#' vst_data_fast <- NormalizeVSTfast(counts)
#' 
#' # Use in downstream analysis
#' pca_result <- prcomp(t(vst_data), scale. = TRUE)
#' }
#'
#' @references
#' Anscombe, F. J. (1948). The transformation of Poisson, binomial and
#' negative-binomial data. Biometrika, 35(3/4), 246-254.
#'
#' @importFrom stats lm resid nls coef
#'
#' @name NormalizeVST
#' @rdname NormalizeVST
#' @export
#' 
NormalizeVST <- function(counts, sv = 1) {
  varx = apply(counts, 1, var)
  meanx = apply(counts, 1, mean)
  phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = sv)))
  
  ## regress out log total counts
  norm_counts <- log(counts + 1/(2 * phi))
  total_counts <- apply(counts, 2, sum)
  res_norm_counts <- t(apply(norm_counts, 1, function(x){resid(lm(x ~ log(total_counts)))} ))
  
  return(res_norm_counts)
}# end func


#' @rdname NormalizeVST
#' @export
#' 
NormalizeVSTfast <- function(counts, sv = 1) {
  varx = as.vector(sp_vars_Rcpp(counts, rowVars=TRUE))
  meanx =  as.vector(sp_means_Rcpp(counts, rowMeans=TRUE))
  phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = sv)))
  
  ## regress out log total counts
  norm_counts <- log(counts + 1/(2 * phi))
  total_counts <- as.vector(sp_sums_Rcpp(counts, rowSums=FALSE))
  
  res_norm_counts <- t(apply(norm_counts, 1, function(x){resid(lm(x ~ log(total_counts)))} ))
  
  return(res_norm_counts)
}# end func



#' Generate chunk points for data processing
#'
#' This internal function generates start and end points for chunking large
#' datasets into manageable pieces for parallel processing or memory-efficient
#' computation. It is particularly useful for processing large matrices or
#' data frames in batches.
#'
#' @param dsize Total size of the data to be chunked (number of elements)
#' @param csize Size of each chunk. If \code{NA}, assumes a single chunk
#'              covering the entire dataset
#'
#' @return A matrix with two columns where each row represents a chunk:
#'   \itemize{
#'     \item Column 1: Start index of the chunk
#'     \item Column 2: End index of the chunk
#'   }
#'   The matrix has as many rows as there are chunks.
#'
#' @details
#' This function is used internally to break down large computational tasks
#' into smaller chunks for:
#' \itemize{
#'   \item Parallel processing across multiple cores
#'   \item Memory-efficient batch processing of large datasets
#'   \item Progress tracking in iterative operations
#'   \item Error handling in partial computations
#' }
#'
#' The function ensures that:
#' \itemize{
#'   \item All data points are covered without gaps
#'   \item Chunks are as evenly sized as possible
#'   \item The last chunk handles any remainder
#'   \item Single-chunk case is handled efficiently
#' }
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Chunk 100 elements into groups of 25
#' chunks <- ChunkPoints(100, 25)
#' print(chunks)
#' #      start end
#' # [1,]     1  25
#' # [2,]    26  50
#' # [3,]    51  75
#' # [4,]    76 100
#'
#' # Single chunk case
#' single_chunk <- ChunkPoints(50, NA)
#' print(single_chunk)
#' #      start end
#' # [1,]     1  50
#'
#' # Chunk with remainder
#' remainder_chunks <- ChunkPoints(23, 10)
#' print(remainder_chunks)
#' #      start end
#' # [1,]     1  10
#' # [2,]    11  20
#' # [3,]    21  23
#'
#' # Use in parallel processing
#' chunks <- ChunkPoints(nrow(large_matrix), 1000)
#' results <- parallel::mclapply(1:nrow(chunks), function(i) {
#'   chunk <- large_matrix[chunks[i, 1]:chunks[i, 2], ]
#'   process_chunk(chunk)
#' }, mc.cores = 4)
#' }
#'
#'
ChunkPoints <- function(dsize, csize) {
  if (is.na(x = csize)) {
    return(matrix(
      data = c(1, dsize),
      ncol = 2,
      dimnames = list(NULL, c('start', 'end'))
    ))
  }
  return(t(x = vapply(
    X = seq.default(from = 1L, to = ceiling(dsize / csize)),
    FUN = function(i) {
      return(c(
        start = (csize * (i - 1L)) + 1L,
        end = min(csize * i, dsize)
      ))
    },
    FUN.VALUE = numeric(length = 2L)
  )))
}## end func


#' Aggregate spots using square binning
#'
#' This function decreases the resolution of spatial transcriptomics data by
#' aggregating neighboring spots into square bins. This is particularly useful
#' for data with very high spatial resolution (e.g., HDST, Slide-seq) where
#' binning can improve signal-to-noise ratio and computational efficiency.
#'
#' @param counts Gene expression counts matrix with features as rows and spots
#'              as columns. Column names should be in the format "location_x x location_y"
#'              or compatible with the spatial coordinates.
#' @param location Spatial location information for each spot. Should be a
#'                data.frame or matrix with rows corresponding to spots and
#'                columns for spatial coordinates (x and y).
#' @param cell.type Cell type information for each spot, particularly useful
#'                 for high-resolution technologies like HDST. Can be a vector
#'                 or matrix with cell type assignments or proportions.
#' @param bin.fold The fold for zooming in (aggregation factor). Must be an
#'                integer (e.g., 2, 3, 4). A value of 2 will combine 2x2 spots
#'                into one bin, 3 will combine 3x3 spots, etc.
#' @param num.core Number of CPU cores to use for parallel processing. Default is 1.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#'
#' @return Returns a list with the following components:
#' \itemize{
#'   \item \code{counts}: Aggregated gene expression counts matrix
#'   \item \code{location}: Updated spatial coordinates for each bin (center points)
#'   \item \code{cell.type}: Aggregated cell type composition (if provided)
#'   \item \code{bin.assignment}: Mapping of original spots to bin identifiers
#' }
#'
#' @details
#' Square binning aggregates neighboring spots into larger square regions,
#' effectively reducing spatial resolution while increasing counts per bin.
#' This approach is beneficial for:
#' \itemize{
#'   \item Reducing technical noise in high-resolution spatial data
#'   \item Improving statistical power for lowly expressed genes
#'   \item Decreasing computational requirements for downstream analysis
#'   \item Matching resolutions across different spatial technologies
#' }
#'
#' The function works by:
#' \enumerate{
#'   \item Dividing the spatial coordinates into a grid of square bins
#'   \item Aggregating counts within each bin (summation)
#'   \item Calculating bin centroids as new spatial coordinates
#'   \item Optionally aggregating cell type information (proportions or counts)
#' }
#'
#' @examples
#' \dontrun{
#' # Load high-resolution spatial data
#' data(hdst_data)
#' 
#' # Aggregate 2x2 spots into bins
#' binned_data <- AggrSquareBin(
#'   counts = hdst_counts,
#'   location = hdst_locations,
#'   bin.fold = 2,
#'   num.core = 4
#' )
#' 
#' # Aggregate with cell type information
#' binned_data_with_ct <- AggrSquareBin(
#'   counts = hdst_counts,
#'   location = hdst_locations,
#'   cell.type = hdst_celltypes,
#'   bin.fold = 3,
#'   verbose = TRUE
#' )
#' 
#' # Access the binned data
#' binned_counts <- binned_data$counts
#' binned_locations <- binned_data$location
#' bin_assignments <- binned_data$bin.assignment
#' 
#' # Compare dimensions before and after binning
#' dim(hdst_counts)      # Original dimensions
#' dim(binned_counts)    # Binned dimensions
#' }
#'
#' @seealso
#' \code{\link{ComputeBandWidth}}
#'
#' @importFrom parallel mclapply
#' 
#' @export
#' 
AggrSquareBin <- function(counts,
                      location,
                      cell.type = NULL, 
                      bin.fold = 2, 
                      num.core = 1,
                      verbose = TRUE){
  
  num_spot <- ncol(counts)
  num_gene <- nrow(counts)

  gene_name <- rownames(counts)
  spot_name <- colnames(counts)
  
  if(!is.null(cell.type) && length(cell.type) != num_spot){
    stop("The 'cell.type' is a character; length of 'cell.type' must be the same as the column of counts. ")
  }## end fi
  if(nrow(location) != num_spot){
    stop("The number of column of 'counts' must be equal to the number of row of 'location'.")
  }## end fi
  if(ncol(location) != 2){
    stop("The number of column of 'location' must be 2!")
  }##
  
  ## convert location to integer, location1
  loc1 <- location[,1]
  names(loc1) <- rownames(location)
  loc1 <- loc1[!duplicated(loc1)]
  #loc1 <- unique(data.frame(location[,1],row.names = rownames(location)))
  order_loc1 <- loc1[order(loc1)]
  #xcoord <- data.frame('xcoord' = c(1:length(order_loc1)), row.names = names(order_loc1))
  xcoord <- data.frame('xcoord' = c(1:length(order_loc1)), row.names = order_loc1)
  
  ## convert location to integer, location2
  loc2 <- location[,2]
  names(loc2) <- rownames(location)
  loc2 <- loc2[!duplicated(loc2)]
  order_loc2 <- loc2[order(loc2)]
  #ycoord <- data.frame('ycoord' = c(1:length(order_loc2)), row.names = names(order_loc2))
  ycoord <- data.frame('ycoord' = c(1:length(order_loc2)), row.names = order_loc2)
  location2 <- data.frame('xcoord' = xcoord[as.character(location[,1]),], 'ycoord' = ycoord[as.character(location[,2]),])
  rownames(location2) <- rownames(location)
  location <- location2
  
  grid_cor <- expand.grid(seq(min(location[,1]), max(location[,1]), by = bin.fold), 
                          seq(min(location[,2]), max(location[,2]), by = bin.fold) )
  
  
  if(!requireNamespace("data.table",quietly = TRUE)){
    stop("Please install package 'data.table'.")
  }
  "%&%" <- function(x, y) paste0(x, y)
  ## assign new coordinates the bottom left corner coordinates
  index <- parallel::mclapply(1:nrow(grid_cor), mc.cores = num.core, function(i){ 
    x <- grid_cor[i, 1]
    y <- grid_cor[i, 2]
    xs <- seq(x, x + bin.fold - 1)
    ys<- seq(y, y + bin.fold - 1)
    
    xy <- expand.grid(xs, ys)
    xy <- xy[, 1] %&% "x" %&% xy[, 2]
    if(length(which(spot_name %in% xy))>0){
      return(which(spot_name %in% xy))
    }else{
      return(NULL)
    } })
  
  ## remove NULL values
  index <- Filter(Negate(is.null), index)
  
  if(verbose){
    cat(paste0("The number of expected samples is: ", length(index), ".\n"))
  }## end fi
  
  ## aggregate the counts
  AggCounts <- lapply(index, function(x){
    if(length(x)>1){
      aggc <- Matrix::rowSums(counts[, x])
    }else{
      aggc <- counts[, x]
    }## end for
    return(aggc)
  })
  aggcounts <- do.call(cbind, AggCounts)
  
  
  ## aggregate locations
  AggLocation <- lapply(index, function(x){
    if(length(x)>1){
      aggl <- as.numeric(Matrix::colMeans(location[x,]))
    }else{
      aggl <- as.numeric(location[x,])
    }## end for
    return(aggl)
  })
  agglocation <- do.call(rbind, AggLocation)
  
  ## aggregate the number of spots
  AggNumSpot <- lapply(index, function(x){
    return(length(x))
  })
  aggnumspot <- do.call(rbind, AggNumSpot)
  
  ## aggregate the percentage the cell types
  if(!is.null(cell.type)) {
    AggFreq <- lapply(index, function(x){
      aggfreq <- data.frame(matrix(0,nrow=length(unique(cell.type))), row.names = as.character(unique(cell.type)))
      tmp <- plyr::count(cell.type[x]);
      aggfreq[as.character(tmp$x),] <- tmp$freq/length(x)
      return(aggfreq)
    })
    aggfreq <- t(do.call(cbind, AggFreq))
  }## end fi
  
  ## column and row names of aggregate counts
  colnames(aggcounts) <- agglocation[, 1] %&% "x" %&% agglocation[, 2]
  counts <- as(aggcounts, "dgCMatrix")
  
  location <- apply(agglocation, 2, as.numeric)
  location <- as.data.frame(location)
  colnames(location) <- c("x","y")
  rownames(location) <- colnames(aggcounts)

  ## cell type information
  if(is.null(cell.type)) {
    meta.data <- data.frame(location, 'aggrnum' = aggnumspot)
  }else{
    rownames(aggfreq) <- colnames(aggcounts)
    meta.data <- data.frame(location, 'aggrnum' = aggnumspot, aggfreq)
  }## end fi
  
  ## return the results
  return(list(counts = counts, meta.data=meta.data))
}## end funcs



#' Select bandwidth for Gaussian kernel
#'
#' This function selects an appropriate bandwidth for Gaussian kernel smoothing
#' in spatial transcriptomics analysis. The bandwidth parameter controls the
#' spatial scale of smoothing and is crucial for proper spatial pattern
#' detection and noise reduction.
#'
#' @param norm_counts A normalized gene expression matrix with genes as rows and
#'                   spatial locations as columns. Should be normalized and
#'                   scaled appropriately for spatial analysis.
#' @param location Optional spatial coordinates data.frame or matrix with rows
#'                corresponding to spatial locations and columns for coordinates
#'                (typically x and y). If provided, uses spatial distances for
#'                bandwidth selection.
#' @param method Method for bandwidth selection. Options include:
#'   \itemize{
#'     \item "Silverman": Silverman's rule of thumb, suitable for large sample sizes
#'     \item "SJ": Sheather-Jones method, more accurate for small sample sizes
#'     \item "Scott": Scott's rule for multivariate data
#'     \item "UCV": Unbiased cross-validation
#'     \item "BCV": Biased cross-validation
#'   }
#'   Default is "Silverman".
#' @param num.core Number of CPU cores to use for parallel computation when
#'                applicable. Default is 2.
#'
#' @return A numeric value representing the calculated bandwidth for Gaussian
#'         kernel smoothing. This value can be used directly in spatial
#'         smoothing functions.
#'
#' @details
#' Bandwidth selection is critical for Gaussian kernel smoothing in spatial
#' transcriptomics:
#' \itemize{
#'   \item Too small bandwidth: Overfitting, noise amplification
#'   \item Too large bandwidth: Over-smoothing, loss of spatial patterns
#' }
#'
#' The available methods have different characteristics:
#' \describe{
#'   \item{Silverman}{Fast rule-of-thumb based on standard deviation and sample size.
#'         Best for normally distributed data and large samples.}
#'   \item{SJ (Sheather-Jones)}{Plug-in method that estimates the optimal bandwidth
#'         by estimating the density's curvature. More accurate for small samples
#'         and non-normal distributions.}
#'   \item{Scott}{Multivariate extension of Silverman's rule. Suitable for
#'         spatial coordinates when provided.}
#'   \item{Cross-validation methods}{Computationally intensive but data-driven
#'         approaches that minimize prediction error.}
#' }
#'
#' When spatial coordinates are provided, the function calculates bandwidth
#' based on spatial distances. Otherwise, it uses expression distances between
#' spatial locations.
#'
#' @examples
#' \dontrun{
#' # Using normalized expression data only
#' normalized_data <- normalize_counts(raw_counts)
#' bandwidth1 <- ComputeBandWidth(normalized_data, method = "Silverman")
#'
#' # Using spatial coordinates for more accurate bandwidth
#' spatial_coords <- GetSpatialCoordinates(spatial_object)
#' bandwidth2 <- ComputeBandWidth(normalized_data, 
#'                               location = spatial_coords, 
#'                               method = "SJ")
#'
#' # Using cross-validation for optimal bandwidth
#' bandwidth3 <- ComputeBandWidth(normalized_data, 
#'                               location = spatial_coords,
#'                               method = "UCV",
#'                               num.core = 4)
#'
#' }
#'
#' @importFrom stats bw.SJ bw.nrd0
#'
#' @export
#' 
ComputeBandWidth <- function(norm_counts, 
                             location = NULL, 
                             method = "Silverman", 
                             num.core = 2){
    
    num_cell <- ncol(norm_counts)
    
    if (method == "SJ") {
      bw_SJ <- parallel::mclapply(1:nrow(norm_counts), FUN = function(x){
        res <- tryCatch({
          bw <- bw.SJ(norm_counts[x, ], method = "dpi")
        }, error=function(e){cat("Gene",x," :",conditionMessage(e), "\n")})
        return(res)
      }, mc.cores = getOption("mc.cores", num.core))
      beta <- unlist(bw_SJ);
      names(beta) <- rownames(norm_counts)
      
        if(FALSE){bw_SJ = c()
        for (i in 1:nrow(norm_counts)) {
            tryCatch({
            bw_SJ[i] = bw.SJ(norm_counts[i, ], method = "dpi")
             }, error=function(e){cat("Gene",i," :",conditionMessage(e), "\n")})
        }## end for
        # beta = median(na.omit(bw_SJ))
        beta <- bw_SJ;
        }
    }else if (method == "Silverman") {
        
        bw_Silverman <- parallel::mclapply(1:nrow(norm_counts), FUN = function(x){
          res <- tryCatch({
            bw <- bw.nrd0(norm_counts[x, ])
          }, error=function(e){cat("Gene",x," :",conditionMessage(e), "\n")})
          return(res)
        }, mc.cores = getOption("mc.cores", num.core))
        beta <- unlist(bw_Silverman);
        names(beta) <- rownames(norm_counts)
        if(FALSE){bw_Silverman = c()
        for (i in 1:dim(norm_counts)[1]) {
            tryCatch({
            bw_Silverman[i] = bw.nrd0(norm_counts[i, ])
            }, error=function(e){cat("Gene",i," :",conditionMessage(e), "\n")})
        }## end for
        # beta = median(na.omit(bw_Silverman))
        beta <- bw_Silverman;
        }
    }else if (method == "Variability") {
      if(is(norm_counts, "dgCMatrix")){
        VarFunc <- FastSparseRowVar;
      }else if(is(norm_counts, "matrix")){
        VarFunc <- FastRowVar;
      }else{
        norm_counts <- Matrix::Matrix(norm_counts, sparse = TRUE)
        VarFunc <- FastSparseRowVar;
      }
      
      var_gene <- VarFunc(norm_counts)
      var_gene <- 1.06*sqrt(var_gene)*num_cell^(-1.0/5)
      if(FALSE){var_gene = c()
        for (i in 1:dim(norm_counts)[1]) {
            tryCatch({
              var_gene[i] = 1.06*sd(norm_counts[i,])*num_cell^(-1.0/5)
            }, error=function(e){cat("Gene",i," :",conditionMessage(e), "\n")})
        }## end fi
      }
      total_gene <- Matrix::rowSums(norm_counts)
      total_gene <- total_gene/sum(total_gene)
      beta <- var_gene*total_gene;
      # beta = median(na.omit(var_gene))
    }else {
      beta <- ComputeKernelParamLessMem(as.matrix(location))$kparam
    }
	return(beta)
}## end func



#-------------------------------------------------------------------------------
#' Wrapper around mclapply to track progress
#' 
#' Based on http://stackoverflow.com/questions/10984556
#' 
#' @param X         a vector (atomic or list) or an expressions vector. Other
#'                  objects (including classed objects) will be coerced by
#'                  ‘as.list’
#' @param FUN       the function to be applied to
#' @param ...       optional arguments to ‘FUN’
#' @param mc.preschedule see mclapply
#' @param mc.set.seed see mclapply
#' @param mc.silent see mclapply
#' @param mc.cores see mclapply
#' @param mc.cleanup see mclapply
#' @param mc.allow.recursive see mclapply
#' @param mc.progress track progress?
#' @param mc.style    style of progress bar (see txtProgressBar)
#'
#' @export
#' 
#' @examples
#' x <- mclapply2(1:1000, function(i, y) Sys.sleep(0.01))
#' x <- mclapply2(1:3, function(i, y) Sys.sleep(1), mc.cores=1)
#' 
#' dat <- lapply(1:10, function(x) rnorm(100)) 
#' func <- function(x, arg1) mean(x)/arg1 
#' mclapply2(dat, func, arg1=10, mc.cores=2)
#'
mclapply2 <- function(X, FUN, ..., 
                      mc.preschedule = TRUE, mc.set.seed = TRUE,
                      mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
                      mc.cleanup = TRUE, mc.allow.recursive = TRUE,
                      mc.progress=TRUE, mc.style=3) {
  if (!is.vector(X) || is.object(X)) X <- as.list(X)
  
  if (mc.progress) {
    f <- fifo(tempfile(), open="w+b", blocking=T)
    p <- parallel:::mcfork()
    pb <- txtProgressBar(0, length(X), style=mc.style)
    setTxtProgressBar(pb, 0) 
    progress <- 0
    if (inherits(p, "masterProcess")) {
      while (progress < length(X)) {
        readBin(f, "double")
        progress <- progress + 1
        setTxtProgressBar(pb, progress) 
      }
      cat("\n")
      parallel:::mcexit()
    }
  }
  tryCatch({
    result <- parallel::mclapply(X, ..., function(...) {
      res <- FUN(...)
      if (mc.progress) writeBin(1, f)
      res
    }, 
    mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,
    mc.silent = mc.silent, mc.cores = mc.cores,
    mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive
    )
    
  }, finally = {
    if (mc.progress) close(f)
  })
  
  return(result)
}## end func


######################
#' Log-likelihood functions for count distribution models
#'
#' These internal functions calculate log-likelihoods for various count
#' distribution models including negative binomial, Poisson, and their
#' zero-inflated variants. They are used internally for maximum likelihood
#' estimation and model comparison in spatial transcriptomics data analysis.
#'
#' @param param,th,mu,lam Distribution parameters (dispersion, mean, etc.)
#' @param y A numeric vector of observed count values
#' @param x A numeric vector of count values to be fitted
#' @param maxiter Maximum number of iterations for optimization algorithms
#'
#' @return For log-likelihood functions: Returns the log-likelihood value (scalar)
#'
#' @return For fitting functions: Returns a named vector with parameter estimates,
#'         log-likelihood, and convergence information
#'
#' @details
#' These functions provide the statistical foundation for modeling count data
#' in spatial transcriptomics:
#'
#' \describe{
#'   \item{Negative Binomial}{Models overdispersed count data with two parameters:
#'         mean (mu) and dispersion (theta)}
#'   \item{Poisson}{Models equidispersed count data with one parameter (lambda)}
#'   \item{Zero-inflated variants}{Account for excess zeros beyond standard distributions}
#' }
#'
#' The negative binomial parameterization uses the following probability mass function:
#' \deqn{P(Y = y) = \frac{\Gamma(y + \theta)}{y! \Gamma(\theta)} \left(\frac{\theta}{\theta + \mu}\right)^\theta \left(\frac{\mu}{\theta + \mu}\right)^y}
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom MASS glm.nb
#' @importFrom stats optim var logLik dpois
NULL

#' @rdname internal-count-models
#' @keywords internal

nb_loglik <- function(param,y){
  ## this version is two parameter estimation,slower than the fix mu estimation.
  ## but the result is same as the fitdistr
  
  # th = param[1]
  # mu = param[2]
  th = max(param[1],1e-10)
  mu = max(param[2],1e-10)
  
  return(-sum((lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(ifelse(th>0, th,1e-10)) + y * 
                 log(mu + (y == 0)) - (th + y) * log(th + mu))))
}## end func


#' log-likelihood for two-parameter negative binomial but fix the mu
#' @rdname internal-count-models
#' @keywords internal

nb_loglik_mom <- function(th, mu, y){
  th = max(th,1e-10)
  mu = max(mu,1e-10)
  return(-sum((lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(ifelse(th>0, th,1e-10)) + y * 
                 log(mu + (y == 0)) - (th + y) * log(th + mu))))
}## end func


#' fitting data with nb through optim function
#' @rdname internal-count-models
#' @keywords internal

fit_nb_optim <- function(x,maxiter= 500){
  m = mean(x)
  v = stats::var(x)
  
  if(v>m){
    size = m^2/(v-m)
  }else{
    size = 100
  }
  
  
  if(length(x)<1000){
    param2est <- c(size,m)
    ## the BFGS method is more stable than Nelder-Mead? need to figure out which is fast and which is better
    ## optim may give error in some case of 10X, similar problem observed in the MASS (fitdistr)
    fitres <- tryCatch({
      opt_out <- optim(param2est,nb_loglik,y=x,control = list(maxit = maxiter),method="BFGS")
      c(opt_out$par,-opt_out$value,1-opt_out$convergence)
    },
    error = function(cond){
      # library(MASS)
      glmfit <- MASS::glm.nb(x~1)
      c(glmfit$theta, exp(glmfit$coefficients),as.numeric(logLik(glmfit)),min(glmfit$converged,is.null(glmfit$th.warn))) 
    })
  }else{
    ## fix the mu through the moment matching to facilitate the estimation
    fitres <- tryCatch({
      opt_out <- optim(size,nb_loglik_mom,mu=m,y=x,control = list(maxit = maxiter), method="Brent",lower=0,upper=1000)
      c(opt_out$par,m,-opt_out$value,1-opt_out$convergence)
    },
    error = function(cond){
      # library(MASS)
      glmfit <- MASS::glm.nb(x~1)
      c(glmfit$theta, exp(glmfit$coefficients),as.numeric(logLik(glmfit)),min(glmfit$converged,is.null(glmfit$th.warn))) 
    })
  }
  
  names(fitres) <- c("theta","mu","llk","convergence")
  return(fitres)
}## end funcs



#' log-likelihood for poisson distribution
#' @rdname internal-count-models
#' @keywords internal

pos_loglik <- function(mu,y){
  # th = param[1]
  # mu = param[2]
  mu = max(mu,1e-10)
  return(-sum(y*log(mu) - mu - lgamma(y + 1)))
}## end func

#' fitting data with poisson through optim function
#' @rdname internal-count-models
#' @keywords internal

fit_pos_optim <- function(x,maxiter= 500){
  ## need to improve for the sparsematrix
  m = mean(x)
  
  opt_out <- optim(m,pos_loglik,y=x,control = list(maxit = maxiter), method="Brent",lower=0,upper=max(x))
  fitres <- c(opt_out$par,-opt_out$value,1-opt_out$convergence)
  names(fitres) <- c("lambda","llk","convergence")
  return(fitres)
}## end funcs


#' Log-likelihood for zero-inflated Poisson distribution
#'
#' This internal function calculates the log-likelihood for a zero-inflated
#' Poisson (ZIP) distribution given the Poisson mean parameter and observed
#' count data. The zero-inflation probability is typically estimated separately
#' or fixed at a predetermined value.
#'
#' @param lam The mean parameter (lambda) for the Poisson distribution component
#' @param y A numeric vector of observed count values to be fitted
#'
#' @return Returns the log-likelihood value (scalar) for the given Poisson
#'         mean parameter and observed data. Higher values indicate better fit.
#'
#' @details
#' The zero-inflated Poisson distribution models count data with excess zeros.
#' However, this particular function calculates the likelihood for the Poisson
#' component only, assuming the zero-inflation probability is handled externally.
#'
#' The standard zero-inflated Poisson probability mass function is:
#' \deqn{P(Y = y) = \pi \cdot I(y=0) + (1-\pi) \cdot \frac{e^{-\lambda} \lambda^y}{y!}}
#'
#' where:
#' \itemize{
#'   \item \eqn{\pi} is the zero-inflation probability (typically estimated separately)
#'   \item \eqn{\lambda} is the Poisson mean parameter (lam)
#'   \item \eqn{I(y=0)} is an indicator function for zero counts
#' }
#'
#' This function focuses on the Poisson component and is often used in
#' conjunction with external zero-inflation probability estimation in
#' two-step fitting procedures.
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#' \dontrun{
#' # Example count data with zero-inflation
#' counts <- c(rep(0, 25), rpois(75, lambda = 3))
#' 
#' # Calculate log-likelihood for different lambda values
#' ll1 <- zip_loglik(lam = 2, y = counts)
#' ll2 <- zip_loglik(lam = 3, y = counts)
#' ll3 <- zip_loglik(lam = 4, y = counts)
#' 
#' # Find maximum likelihood estimate for lambda
#' optimize(
#'   f = function(lam) zip_loglik(lam, counts),
#'   interval = c(0.1, 10),
#'   maximum = TRUE
#' )
#' 
#' # Use in EM algorithm or two-step ZIP fitting
#' # Step 1: Estimate zero-inflation probability
#' p0_est <- mean(counts == 0) - dpois(0, lambda = mean(counts[counts > 0]))
#' p0_est <- max(0, p0_est)  # Ensure non-negative
#' 
#' # Step 2: Estimate Poisson mean using this function
#' lambda_est <- optimize(
#'   f = function(lam) zip_loglik(lam, counts),
#'   interval = c(0.1, 10),
#'   maximum = TRUE
#' )$maximum
#' }

zip_loglik <- function(lam,y){
  zeroidx <- which(y==0)
  nzero   <- length(zeroidx)
  numSam  <- length(y)
  
  if(nzero>0){
    nz_y <- y[-zeroidx]
  }else{
    nz_y <- y
  }
  
  pp <- (nzero/numSam-exp(-lam))/(1-exp(-lam))
  
  ll <- nzero*log(pp + (1-pp)*exp(-lam)) + sum(log(1-pp) + nz_y*log(lam) - lam -lgamma(nz_y + 1))
  
  return(-(ll))
}## end funcs



#' Fit zero-inflated Poisson distribution using optimization
#'
#' This internal function fits a zero-inflated Poisson (ZIP) distribution
#' to count data using maximum likelihood estimation via the \code{optim} function.
#' It handles zero-inflation commonly found in count data from spatial 
#' transcriptomics and single-cell RNA sequencing when overdispersion is not present.
#'
#' @param x A numeric vector of count values to be fitted
#' @param maxiter Maximum number of iterations for the optimization algorithm
#'               (default: 500)
#'
#' @return Returns a numeric vector with the following components:
#' \itemize{
#'   \item \code{p0}: Zero-inflation proportion (probability of structural zeros)
#'   \item \code{mu}: Mean parameter of the Poisson component
#'   \item \code{llk}: Log-likelihood value at the maximum
#'   \item \code{convergence}: Convergence code from \code{optim} (0 indicates success)
#'   \item \code{zip}: Logical indicating if zero-inflated model was fitted (TRUE)
#' }
#'
#' @details
#' The zero-inflated Poisson model accounts for two types of zeros:
#' \itemize{
#'   \item Structural zeros (due to zero-inflation)
#'   \item Sampling zeros (from the Poisson distribution)
#' }
#'
#' The probability mass function is:
#' \deqn{P(X = k) = \pi \cdot I(k=0) + (1-\pi) \cdot f_{P}(k; \mu)}
#' where:
#' \itemize{
#'   \item \eqn{\pi} is the zero-inflation probability (p0)
#'   \item \eqn{f_{P}} is the Poisson probability mass function
#'   \item \eqn{\mu} is the mean of the Poisson component
#' }
#'
#' The function uses the \code{optim} function with L-BFGS-B method for bounded
#' optimization, ensuring parameters stay within valid ranges (p0 in [0,1], mu > 0).
#'
#' Compared to the zero-inflated negative binomial (ZINB) model, ZIP is more
#' appropriate when there is zero-inflation but no overdispersion beyond what
#' the Poisson distribution can accommodate.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats optim
#'
#' @examples
#' \dontrun{
#' # Generate example count data with zero-inflation (no overdispersion)
#' set.seed(123)
#' counts <- c(rep(0, 40), rpois(160, lambda = 3))
#' 
#' # Fit ZIP model
#' fit <- fit_zip_optim(counts, maxiter = 1000)
#' 
#' # Extract parameters
#' p0 <- fit["p0"]        # Zero-inflation proportion
#' mu <- fit["mu"]        # Mean parameter
#' llk <- fit["llk"]      # Log-likelihood
#' 
#' # Check convergence
#' if (fit["convergence"] == 0) {
#'   print("Model converged successfully")
#' }
#' 
#' # Compare with regular Poisson (no zero-inflation)
#' poisson_fit <- fit_poisson_optim(counts)
#' }

fit_zip_optim <- function(x,maxiter= 500){
  
  ## need to improve for the sparsematrix
  zeroidx <- which(x==0)
  n0    <- length(zeroidx)
  n     <- length(x)
  m     <- mean(x)
  
  opt_out <- optim(m,zip_loglik,y=x,control = list(maxit = maxiter), method="Brent",lower=min(-log(n0/n),0),upper=max(x))
  
  pp_est  <- (n0/n-exp(-opt_out$par))/(1-exp(-opt_out$par))
  
  if(pp_est<=0|pp_est==1|opt_out$convergence==1){
    opt_out   <- fit_pos_optim(x,maxiter= maxiter)
    fitres    <- c(0.0,opt_out,0)
  }else{
    fitres  <- c(pp_est,opt_out$par,-opt_out$value,1-opt_out$convergence,1)
  }
  names(fitres) <- c("p0","mu","llk","convergence","zip")
  return(fitres)
}## end func



#' Log-likelihood for zero-inflated negative binomial distribution
#'
#' This internal function calculates the log-likelihood for a zero-inflated
#' negative binomial (ZINB) distribution given parameter values and observed
#' count data. It is used internally by optimization routines for maximum
#' likelihood estimation of ZINB parameters.
#'
#' @param par_init A numeric vector of three parameters for the ZINB distribution:
#'   \itemize{
#'     \item \code{p0}: Zero-inflation probability (probability of structural zeros)
#'     \item \code{theta}: Dispersion parameter (inverse of overdispersion)
#'     \item \code{mu}: Mean parameter of the negative binomial component
#'   }
#' @param y A numeric vector of observed count values to be fitted
#'
#' @return Returns the log-likelihood value (scalar) for the given parameters
#'         and observed data. Higher values indicate better fit.
#'
#' @details
#' The zero-inflated negative binomial distribution models count data with
#' excess zeros and overdispersion. The probability mass function is:
#'
#' \deqn{P(Y = y) = \pi \cdot I(y=0) + (1-\pi) \cdot f_{NB}(y; \mu, \theta)}
#'
#' where:
#' \itemize{
#'   \item \eqn{\pi} is the zero-inflation probability (p0)
#'   \item \eqn{I(y=0)} is an indicator function for zero counts
#'   \item \eqn{f_{NB}} is the negative binomial probability mass function
#'   \item \eqn{\mu} is the mean of the negative binomial component
#'   \item \eqn{\theta} is the dispersion parameter
#' }
#'
#' The log-likelihood is calculated as the sum of log probabilities across
#' all observations. The function handles edge cases and ensures numerical
#' stability through careful implementation.
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#' \dontrun{
#' # Example count data
#' counts <- c(rep(0, 30), rnbinom(70, size = 2, mu = 5))
#' 
#' # Test parameters
#' test_params <- c(p0 = 0.3, theta = 2, mu = 5)
#' 
#' # Calculate log-likelihood
#' ll <- zinb_loglik(test_params, counts)
#' print(ll)
#' 
#' # Use in optimization (minimize negative log-likelihood)
#' optim_result <- optim(
#'   par = c(0.2, 1.5, 4),
#'   fn = function(par) -zinb_loglik(par, counts),
#'   method = "L-BFGS-B",
#'   lower = c(0, 0.1, 0.1),
#'   upper = c(1, 100, 100)
#' )
#' }
zinb_loglik <- function(par_init,y){
  
  th      <- max(par_init[1], 1e-10)
  # mu      <- par_init[2]
  mu      <- max(par_init[2], 1e-10)
  rr      <- par_init[3]
  
  zeroidx <- which(y==0)
  nzero   <- length(zeroidx)
  numSam  <- length(y)
  
  if(nzero>0){
    nz_y <- y[-zeroidx]
  }else{
    nz_y <- y
  }
  
  # pp <- (nzero/numSam-(th/(mu+th))^th)/(1-(th/(mu+th))^th)
  pp  <- exp(rr)/(1+exp(rr))
  
  ll1 <- nzero*log(pp + (1-pp)*(th/(mu+th))^th) 
  ll2 <- sum(log(1-pp)+(lgamma(th + nz_y) - lgamma(th) - lgamma(nz_y + 1) + th * log(ifelse(th>0, th,1e-10)) + nz_y * 
                          log(mu) - (th + nz_y) * log(th + mu)))
  ll  <- ll1 + ll2
  
  return(-ll)
}## end func


#' Fit zero-inflated negative binomial distribution using optimization
#'
#' This internal function fits a zero-inflated negative binomial (ZINB) distribution
#' to count data using maximum likelihood estimation via the \code{optim} function.
#' It handles zero-inflation and overdispersion commonly found in count data
#' from spatial transcriptomics and single-cell RNA sequencing.
#'
#' @param x A numeric vector of count values to be fitted
#' @param maxiter Maximum number of iterations for the optimization algorithm
#'               (default: 500)
#'
#' @return Returns a numeric vector with the following components:
#' \itemize{
#'   \item \code{p0}: Zero-inflation proportion (probability of structural zeros)
#'   \item \code{theta}: Dispersion parameter of the negative binomial component
#'   \item \code{mu}: Mean parameter of the negative binomial component
#'   \item \code{llk}: Log-likelihood value at the maximum
#'   \item \code{convergence}: Convergence code from \code{optim} (0 indicates success)
#'   \item \code{zinb}: Logical indicating if zero-inflated model was fitted (TRUE)
#' }
#'
#' @details
#' The zero-inflated negative binomial model accounts for two types of zeros:
#' \itemize{
#'   \item Structural zeros (due to zero-inflation)
#'   \item Sampling zeros (from the negative binomial distribution)
#' }
#'
#' The probability mass function is:
#' \deqn{P(X = k) = \pi \cdot I(k=0) + (1-\pi) \cdot f_{NB}(k; \mu, \theta)}
#' where:
#' \itemize{
#'   \item \eqn{\pi} is the zero-inflation probability (p0)
#'   \item \eqn{f_{NB}} is the negative binomial probability mass function
#'   \item \eqn{\mu} is the mean of the negative binomial component
#'   \item \eqn{\theta} is the dispersion parameter (inverse of overdispersion)
#' }
#'
#' The function uses the \code{optim} function with L-BFGS-B method for bounded
#' optimization, ensuring parameters stay within valid ranges (p0 in [0,1],
#' theta and mu > 0).
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats optim
#'
#' @examples
#' \dontrun{
#' # Generate example count data with zero-inflation
#' set.seed(123)
#' counts <- c(rep(0, 50), rnbinom(150, size = 2, mu = 5))
#' 
#' # Fit ZINB model
#' fit <- fit_zinb_optim(counts, maxiter = 1000)
#' 
#' # Extract parameters
#' p0 <- fit["p0"]        # Zero-inflation proportion
#' theta <- fit["theta"]  # Dispersion parameter
#' mu <- fit["mu"]        # Mean parameter
#' llk <- fit["llk"]      # Log-likelihood
#' 
#' # Check convergence
#' if (fit["convergence"] == 0) {
#'   print("Model converged successfully")
#' }
#' }

fit_zinb_optim <- function(x,
                           maxiter= 500){
  
  ## need to improve for the sparsematrix
  
  zeroidx <- which(x==0)
  n0    <- length(zeroidx)
  n     <- length(x)
  m     <- mean(x)
  v     <- stats::var(x)
  if(v>m){
    size = m^2/(v-m)
  }else{
    size = 1
  }
  
  if(n0!=0){
    param_init <- c(size,m,log(n0/(n-n0)))
    
    opt_out <- optim(param_init,zinb_loglik,y=x,control = list(maxit = maxiter), method="Nelder-Mead")
    pp_est <- exp(opt_out$par[3])/(1+exp(opt_out$par[3]))
    
    # if not converge or the proportion is problematic, then switch to NB
    if(pp_est<=0|pp_est==1|opt_out$convergence==1){
      opt_out   <- fit_nb_optim(x,maxiter= maxiter)
      fitres    <- c(0.0,opt_out,0)
    }else{
      fitres    <- c(pp_est,opt_out$par[1],opt_out$par[2],-opt_out$value,1-opt_out$convergence,1)
    }
  }else{
    opt_out   <- fit_nb_optim(x,maxiter= maxiter)
    fitres    <- c(0.0,opt_out,0)
  }
  
  names(fitres) <- c("p0","theta","mu","llk","convergence","zinb")
  return(fitres)
}## end func


#' Check if a matrix is empty
#'
#' This function checks whether a matrix is empty, defined as either having
#' zero dimensions (0x0) or being a 1x1 matrix containing only NA values.
#' Useful for validating matrix inputs in functions and pipelines.
#'
#' @param x A matrix or matrix-like object to check for emptiness
#'
#' @return Logical value indicating whether the matrix is empty (TRUE) or
#'         contains data (FALSE)
#'
#' @export
#'
#' @concept utils
#'
#' @examples
#' \dontrun{# Check various matrix types
#' IsMatrixEmpty(new("matrix"))  # 0x0 matrix
#' IsMatrixEmpty(matrix())       # 0x0 matrix  
#' IsMatrixEmpty(matrix(1:3))    # Non-empty matrix
#' IsMatrixEmpty(matrix(NA))     # 1x1 NA matrix
#' IsMatrixEmpty(matrix(NA, nrow = 0, ncol = 0))  # 0x0 matrix
#'
#' # Practical usage in functions
#' my_matrix <- matrix(1:4, nrow = 2)
#' if (!IsMatrixEmpty(my_matrix)) {
#'   print("Matrix contains data")
#' }
#'}
#'
IsMatrixEmpty <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
}## end func

#' Extract delimiter information from a string
#'
#' Parses a string (usually a cell name) and extracts specified fields based
#' on a given delimiter. This function is particularly useful for processing
#' cell barcodes or feature names that follow a structured naming convention.
#'
#' @param string Character string to parse. Typically a cell barcode or 
#'              feature identifier with delimited components
#' @param field Integer(s) indicating which field(s) to extract. Can be a
#'              single integer or a vector of multiple numbers. Fields are
#'              1-indexed (first field is 1)
#' @param delim Delimiter character to use for splitting the string. 
#'             Set to underscore by default
#'
#' @return A character string containing the extracted field(s). If multiple
#'         fields are specified, they are rejoined using the same delimiter
#'
#' @export
#'
#' @examples
#' # Extract single field from cell barcode
#' cell_name <- "sample1_batchA_AAACCTGAGCTATGCT"
#' ExtractField(cell_name, field = 1)  # Returns "sample1"
#' ExtractField(cell_name, field = 3)  # Returns "AAACCTGAGCTATGCT"
#'
#' # Extract multiple fields and rejoin
#' ExtractField(cell_name, field = c(1, 3))  # Returns "sample1_AAACCTGAGCTATGCT"
#'
#' # Use with different delimiters
#' cell_name2 <- "sample1-batchA-AAACCTGAGCTATGCT"
#' ExtractField(cell_name2, field = 2, delim = "-")  # Returns "batchA"
#'
#' # Practical usage in processing cell names
#' cell_names <- c("patient1_treatment_CTACGTG", "patient1_control_AGTCCTG")
#' samples <- sapply(cell_names, ExtractField, field = 2, USE.NAMES = FALSE)
#' print(samples)  # Returns "treatment" "control"
#'
#' # Extract multiple fields from multiple strings
#' fields <- lapply(cell_names, ExtractField, field = c(1, 3))
#' print(unlist(fields))  # Returns "patient1_CTACGTG" "patient1_AGTCCTG"
#'
ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field),
                                               split = "," )))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim ))
}## end func

