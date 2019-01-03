#ifndef ACTIONcore_H
#define ACTIONcore_H

#include <omp.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <unordered_map>
#include <ctime>

#define ARMA_DONT_USE_WRAPPER
#undef ARMA_BLAS_CAPITALS
#define ARMA_BLAS_UNDERSCORE
#define ARMA_64BIT_WORD
#define ARMA_BLAS_LONG_LONG

#include "armadillo"

using namespace std;
using namespace arma;

struct SPA_results {
	uvec selected_columns;
	vec column_norms;
};


struct ACTION_results {
	vector<uvec> selected_cols;
	vector<mat> H;
	vector<mat> C;
};

namespace ACTIONcore {
	SPA_results SPA(mat M, int k); // Solve convex-NMF using the Successive Projection Algorithm (SPA)
	void simplexRegression(mat &A, mat &B, double *X_ptr); // min_{X} (|| AX - B ||) s.t. simplex constraint
	field<mat> AA (mat X, mat Z0); // Robust archetypal analysis method
	ACTION_results runACTION(mat S_r, int k_min, int k_max, int numThreads); // Main ACTION function	
}
#endif
