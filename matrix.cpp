// tridiagonal matrix solver
//https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm

#include <iostream>
#include <cmath>
#include <cassert>
#include "matrix.h"
#include "cn_solver.h"

// a is lower / upper diagonal vector as matrix is symmetrix. real, only (n-1)
// b is diagonal. complex
// d is output vector
vector<complex<double>> tridiag_solve(vector<double> a, vector<complex<double>> b, vector<complex<double>> d) {
	assert(b.size() == d.size());
	assert(a.size() == d.size() - 1);
	int n = d.size();
	vector<complex<double>> x(n); //output vec
	
	for (int j=1; j < n; j++) {
		complex<double> w = a[j-1] / b[j-1];
		b[j] = b[j] - w * a[j-1];
		d[j] = d[j] - w * d[j-1];
	}
	
	x[n-1] = d[n-1] / b[n-1];
	for (int j=n-2; j >= 0; j--) {
		x[j] = (d[j] - a[j] * x[j+1]) / b[j];
	}
	
	return x;
}


void test_ex() {
	vector<double> a(1);
	vector<complex<double>> b(2);
	vector<complex<double>> d(2);

	a[0] = 1.0;
	b[0] = {1.0, 1.0};
	b[1] = {1.0, -1.0};
	d[0] = {0.0, 2.0};
	d[1] = {1.0, 1.0};
		
	vector<complex<double>> x = tridiag_solve(a, b, d);
	
	cout << "hello world" << endl;
	cout << x[0] << endl;
	cout << x[1];
}