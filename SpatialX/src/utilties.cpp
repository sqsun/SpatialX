
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <unordered_map>
#include <fstream>
#include <string>
#include <Rinternals.h>
#include <iomanip>
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>
#include <ctime>
#include <Rcpp.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace RcppArmadillo;
//[[Rcpp::depends(RcppArmadillo)]]




// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;

//' Compute Euclidean distance matrix by rows
//' @param x A location matrix
// [[Rcpp::export]]
arma::mat ED_cpp(const arma::mat & x) {
  unsigned int outrows = x.n_rows, i = 0, j = 0;
  double d;
  arma::mat out = zeros<arma::mat>(outrows, outrows);
  
  for (i = 0; i < outrows - 1; i++) {
    arma::rowvec v1 = x.row(i);
    for (j = i + 1; j < outrows; j++) {
      d = sqrt(sum(pow(v1 - x.row(j), 2.0)));
      out(j, i) = d;
      out(i, j) = d;
    }
  }
  return out;
}// end func
//' Do inverse of sysmetric matrix 
//' @param Min A sysmetric matrix
//' 
//' @return A list
//' 
//' @export
// [[Rcpp::export]]
SEXP MatrixInverse(SEXP Min) {
  try {
    arma::mat M = as<mat>(Min);
    arma::vec eigval = zeros<vec>( M.n_rows );
    arma::mat eigvec = zeros<mat>( size(M) );
    eig_sym(eigval, eigvec, M, "dc");
    const uvec idx = find(eigval < 1e-8 );
    arma::vec tmp_value = ones<arma::vec>(idx.n_elem);
    eigval.elem( idx ) = tmp_value * 1e-8;
    arma::mat M1 = eigvec.each_row() % (1.0/eigval.t());
    M = M1 * eigvec.t();
    // return values
    return List::create(Named("eigval") = eigval, Named("eigvec") = eigvec, Named("inv_mat") = M);
    //return List::create(Named("eigval") = eigval, Named("eigvec") = eigvec);
  }
  catch (std::exception &ex)
  {
    forward_exception_to_r(ex);
  }
  catch (...)
  {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
}// end func


//' Do inverse of sysmetric matrix 
//' @param Min A sysmetric matrix
//' 
//' @return A list
//' 
//' @export
// [[Rcpp::export]]
SEXP SysMatEigen(SEXP Min) {
  try {
    arma::mat M = as<arma::mat>(Min);
    arma::vec eigval = zeros<arma::vec>( M.n_rows );
    arma::mat eigvec = zeros<arma::mat>( size(M) );
    eig_sym(eigval, eigvec, M, "dc");
    const uvec idx = find(eigval < 1e-8 );
    arma::vec tmp_value = ones<arma::vec>(idx.n_elem);
    eigval.elem( idx ) = tmp_value * 1e-8;
    arma::mat M1 = eigvec.each_row() % eigval.t();
    M = M1 * eigvec.t();
    // return values
    return List::create(Named("eigval") = eigval, Named("eigvec") = eigvec, Named("kernel_mat") = M);
    //return List::create(Named("eigval") = eigval, Named("eigvec") = eigvec);
  }
  catch (std::exception &ex)
  {
    forward_exception_to_r(ex);
  }
  catch (...)
  {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
}// end func


//' Do inverse of sparse sysmetric matrix 
//' @param Min A sysmetric matrix
//' @param num_topin The number of top eigen values
//' 
//' @return A list
//' 
//' @export
// [[Rcpp::export]]
SEXP SparseSysMatEigen(SEXP Min, SEXP num_topin) {
  try {
    arma::sp_mat M = as<sp_mat>(Min);
    int num_top = Rcpp::as<int>(num_topin);
    arma::vec eigval = zeros<arma::vec>( num_top );
    arma::mat eigvec = zeros<arma::mat>(M.n_rows, num_top );
    eigs_sym(eigval, eigvec, M, num_top);
    const uvec idx = find(eigval < 1e-8 );
    arma::vec tmp_value = ones<arma::vec>(idx.n_elem);
    eigval.elem( idx ) = tmp_value * 1e-8;
    arma::mat M1 = eigvec.each_row() % eigval.t();
    M = M1 * eigvec.t();
    // return values
    return List::create(Named("eigval") = eigval, Named("eigvec") = eigvec, Named("kernel_mat") = M);
    //return List::create(Named("eigval") = eigval, Named("eigvec") = eigvec);
  }
  catch (std::exception &ex)
  {
    forward_exception_to_r(ex);
  }
  catch (...)
  {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
}// end func

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> LogLSFactorCpp(Eigen::SparseMatrix<double> data, int scale_factor, bool display_progress = true){
  Progress p(data.outerSize(), display_progress);
  Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      it.valueRef() = log1p(double(it.value()) / colSums[k] * scale_factor);
      //it.valueRef() = log(( (double(it.value()) / colSums[k]) +pseudo_count)* scale_factor);
    }
  }
  return data;
}// end func



// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> LogSizeFactorCpp(Eigen::SparseMatrix<double> data, Eigen::VectorXd size_factor, int pseudo_count, bool display_progress = true){
  Progress p(data.outerSize(), display_progress);
  
  //Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
  
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      it.valueRef() = log(double(it.value()) / size_factor[k] + pseudo_count);
    }
  }// end for
  return data;
}// end func


// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd FastSparseRowScaleData(Eigen::SparseMatrix<double> mat, bool scale = true, bool center = true,
                                   double scale_max = 10, bool display_progress = true){
  mat = mat.transpose();
  Progress p(mat.outerSize(), display_progress);
  Eigen::MatrixXd scaled_mat(mat.rows(), mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double colMean = 0;
    double colSdev = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
    {
      colMean += it.value();
    }
    colMean = colMean / mat.rows();
    if (scale == true){
      int nnZero = 0;
      if(center == true){
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
          nnZero += 1;
          colSdev += pow((it.value() - colMean), 2);
        }
        colSdev += pow(colMean, 2) * (mat.rows() - nnZero);
      }
      else{
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
          colSdev += pow(it.value(), 2);
        }
      }
      colSdev = sqrt(colSdev / (mat.rows() - 1));
    }
    else{
      colSdev = 1;
    }
    if(center == false){
      colMean = 0;
    }
    Eigen::VectorXd col = Eigen::VectorXd(mat.col(k));
    scaled_mat.col(k) = (col.array() - colMean) / colSdev;
    for(int s=0; s<scaled_mat.col(k).size(); ++s){
      if(scaled_mat(s,k) > scale_max){
        scaled_mat(s,k) = scale_max;
      }
    }
  }
  return scaled_mat.transpose();
}


// [[Rcpp::export(rng = false)]]
NumericMatrix FastRowScaleData(Eigen::Map<Eigen::MatrixXd> mat, bool display_progress = true){
  Progress p(mat.cols(), display_progress);
  NumericMatrix std_mat(mat.rows(), mat.cols());
  for(int i=0; i < mat.cols(); ++i){
    p.increment();
    Eigen::ArrayXd r = mat.col(i).array();
    double colMean = r.mean();
    NumericMatrix::Column new_col = std_mat(_, i);
    for(int j=0; j < new_col.size(); j++) {
      new_col[j] = (r[j] - colMean) ;
    }// end for
  }// end for
  return std_mat;
}// end func


//' Compute Euclidean distance matrix by columns
//'
//' Used in sc3-funcs.R distance matrix calculation
//' and within the consensus clustering.
//'
//' @param x A numeric matrix.
// [[Rcpp::export]]
Rcpp::NumericMatrix ED2(const Rcpp::NumericMatrix & x) {
	unsigned int outcols = x.ncol(), i = 0, j = 0;
	double d;
	Rcpp::NumericMatrix out(outcols, outcols);

	for (j = 0; j < outcols - 1; j++) {
		Rcpp::NumericVector v1 = x.column(j);
		for (i = j + 1; i < outcols; i++) {
			d = sqrt(sum(pow(v1 - x.column(i), 2.0)));
			out(i, j) = d;
			out(j, i) = d;
		}
	}
	return out;
}// end func


//' compute D-alpha*W for CAR prior
//'
//'
//' @param A symmetric matrix
//' @param sigma Kernel parameter
//' @param alpha Tunning parameter
//' @param thr_car Threashold to cut
//' @export
// [[Rcpp::export]]
arma::mat ComputeDmW(arma::mat A, double sigma, double alpha, double thr_car) {
	// gaussian kernel 
	A = exp(-A/(sigma*A.max()) );
	//A = exp(-A/(2*sigma) );
	A.diag().zeros();
	A.elem( find(A > thr_car) ).ones();
	
	// quantity D^{1/2}
	arma::rowvec D_row = pow(sum(A), -0.5);
	A.each_row() %= D_row;
	arma::colvec D_col = conv_to< colvec >::from(D_row);
	A.each_col() %= D_col;
	arma::mat res = eye(A.n_cols, A.n_cols) - alpha*A;
	return(res);
}// end func


/* Performs column scaling and/or centering. Equivalent to using scale(mat, TRUE, apply(x,2,sd)) in R.
 Note: Doesn't handle NA/NaNs in the same way the R implementation does, */

// [[Rcpp::export(rng = false)]]
NumericMatrix Standardize(Eigen::Map<Eigen::MatrixXd> mat, bool display_progress = true){
  Progress p(mat.cols(), display_progress);
  NumericMatrix std_mat(mat.rows(), mat.cols());
  for(int i=0; i < mat.cols(); ++i){
    p.increment();
    Eigen::ArrayXd r = mat.col(i).array();
    double colMean = r.mean();
    double colSdev = sqrt((r - colMean).square().sum() / (mat.rows() - 1));
    NumericMatrix::Column new_col = std_mat(_, i);
    for(int j=0; j < new_col.size(); j++) {
      new_col[j] = (r[j] - colMean) / colSdev;
    }
  }
  return std_mat;
}// end func


/* Calculates the variance of rows of a matrix */
// [[Rcpp::export(rng = false)]]
NumericVector FastRowVar(Eigen::Map<Eigen::MatrixXd> x){
  NumericVector out(x.rows());
  for(int i=0; i < x.rows(); ++i){
    Eigen::ArrayXd r = x.row(i).array();
    double rowMean = r.mean();
    out[i] = (r - rowMean).square().sum() / (x.cols() - 1);
  }
  return out;
}// end func

// [[Rcpp::export(rng = false)]]
NumericVector FastRowMean(Eigen::Map<Eigen::MatrixXd> x){
  NumericVector out(x.rows());
  for(int i=0; i < x.rows(); ++i){
    Eigen::ArrayXd r = x.row(i).array();
    out[i] =  r.mean();
  }// end for
  return out;
}// end func

/* Calculate the variance in non-logspace (return answer in non-logspace) */
// [[Rcpp::export(rng = false)]]
Eigen::VectorXd FastSparseRowVar(Eigen::SparseMatrix<double> mat){
  int ncols = mat.cols();
  Eigen::VectorXd rowdisp(mat.rows());
  mat = mat.transpose();
  
  bool display_progress = true;
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating gene variances" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double rm = 0;
    double v = 0;
    int nnZero = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      rm += (it.value());
    }
    rm = rm / ncols;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      v += pow((it.value()) - rm, 2);
      nnZero += 1;
    }
    v = (v + (ncols - nnZero) * pow(rm, 2)) / (ncols - 1);
    rowdisp[k] = v;
  }
  return(rowdisp);
}// end func


//Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
// [[Rcpp::export(rng = false)]]
Eigen::VectorXd FastSparseRowMean(Eigen::SparseMatrix<double> mat){
  Eigen::VectorXd rowSums = mat * Eigen::VectorXd::Ones(mat.cols());
  return(rowSums/mat.cols());
}// end func

//------------------------------------
// rowMeans, colMeans of a sparse matrix 
//------------------------------------

// [[Rcpp::export]]
arma::rowvec sp_means_Rcpp(arma::sp_mat sp_data, bool rowMeans = false) {

  arma::sp_mat norm_col_sums;

  arma::mat tmp_mat;

  if (rowMeans) {

    norm_col_sums = arma::mean(sp_data, 1);

    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.col(0));}

  else {

    norm_col_sums = arma::mean(sp_data, 0);

    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.row(0));
  }

  arma::rowvec tmp_vec = arma::conv_to< arma::rowvec >::from(tmp_mat);

  return tmp_vec;
}



//------------------------------------
// rowSums, colSums of a sparse matrix
//------------------------------------


// [[Rcpp::export]]
arma::rowvec sp_sums_Rcpp(arma::sp_mat sp_data, bool rowSums = false) {

  arma::mat tmp_mat;

  arma::sp_mat norm_col_sums;

  if (rowSums) {

    norm_col_sums = arma::sum(sp_data, 1);

    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.col(0));}

  else {

    norm_col_sums = arma::sum(sp_data, 0);

    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.row(0));
  }

  arma::rowvec tmp_vec = arma::conv_to< arma::rowvec >::from(tmp_mat);

  return tmp_vec;
}

// [[Rcpp::export]]
arma::rowvec sp_vars_Rcpp(arma::sp_mat sp_data, bool rowVars = false) {

  arma::sp_mat norm_col_Vars;

  arma::mat tmp_mat;

  if (rowVars) {

    norm_col_Vars = arma::var(sp_data,0, 1);

    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_Vars.col(0));}

  else {

    norm_col_Vars = arma::var(sp_data, 0,0);

    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_Vars.row(0));
  }

  arma::rowvec tmp_vec = arma::conv_to< arma::rowvec >::from(tmp_mat);

  return tmp_vec;
}

