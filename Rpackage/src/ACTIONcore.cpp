#include <ACTIONcore.h>
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

#define ARMA_USE_CXX11_RNG

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat runsimplexRegression(mat A, mat B) {	

	mat X(A.n_cols, B.n_cols);
	ACTIONcore::simplexRegression(A, B, X.memptr());
				
	return X;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat runSPA(mat A, int k) {	

	SPA_results res = ACTIONcore::SPA(A, k);
	uvec selected_columns = res.selected_columns;
	
	vec cols(k);
	for(int i = 0; i < k; i++) {
		cols[i] = selected_columns[i] + 1;
	}
			
	return cols;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List runAA(mat A, mat W0) {	

	field<mat> decomposition = ACTIONcore::AA (A, W0);
	mat C = decomposition(0);
	mat H = decomposition(1);
	
	List out_list;		
	out_list["C"] = C;
	out_list["H"] = H;
			
	return out_list;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List runACTION(mat cell_signatures, int k_min, int k_max, int numThreads) {	

	ACTION_results trace = ACTIONcore::runACTION(cell_signatures, k_min, k_max, numThreads);

	List res;
	
	List C(k_max);
	for (int i = k_min; i <= k_max; i++) {
		C[i-1] = trace.C[i];
	}
	res["C"] = C;	

	List H(k_max);
	for (int i = k_min; i <= k_max; i++) {
		H[i-1] = trace.H[i];
	}
	res["H"] = H;
	
		
	return res;
}
