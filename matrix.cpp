// crout reduction for tridiagonal matrix
//https://en.wikipedia.org/wiki/Crout_matrix_decomposition


// A is tridiagonal input matrix
// L is output for lower triangluar matrix
// U is output for upper triangular
// n x n matrix input
void crout(double** A, double** L, double** U, int n) {
	int i,j,k;
	double sum = 0;
	
	for (i = 0; i < n; i++) {
		U[i][i] = 1;
	}

	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}

		for (i = j; i < n; i++) {
			sum = 0;
			for(k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {
				// determinant 0
				exit(EXIT_FAILURE);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
	}
}
