#ifndef MATRIX_H 
#define MATRIX_H

#include <vector>
#include <complex>
using namespace std;


vector<complex<double>> tridiag_solve(vector<double> a, vector<complex<double>> b, vector<complex<double>> d);

#endif