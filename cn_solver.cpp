// reference https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method
#include <complex>
#include <cmath>
#include <iostream>
#include <cstring>
#include "matrix.h"
#include "cn_solver.h"

const double C1 = dt / (dx * dx * 2.0);


double totalsum(vector<complex<double>> &phi) {
	double sum = 0;
	for (int i = 0; i < NX; i++) {
		sum += norm(phi[i]);
	}
	return sum;
}


vector<double> normalise(vector<complex<double>> &phi) {
	vector<double> normalised(NX);
	double sum = 0;
	for (int i = 0; i < NX; i++) {
		sum += norm(phi[i]);
	}
	for (int i = 0; i < NX; i++) {
		normalised[i] = norm(phi[i]) / sum;
	}
	
	return normalised;
 }


void norm_wavefunc(vector<complex<double>> &phi) {
	double sum = 0;
	for (int i=0; i < NX; i++) {
		sum += norm(phi[i]);
	}
	
	for (int i=0; i < NX; i++) {
		phi[i] = phi[i] / sqrt(sum);
	}
}


void log_norm(vector<complex<double>> &phi, int f) {
	FILE * fp;
	//char fname[20] = "phi_";
	//fname[4] = 'a' + f;
	//strcpy(fname + 5, ".data");

	fp = fopen("phi_n.data","a");
	fprintf(fp, "#%d\n", f);
	vector<double> normd = normalise(phi);
	for (int i=0; i<NX; i++) {
		fprintf(fp, "%.3e\n", normd[i]);
	}
	fprintf(fp, "\n\n");
	fclose(fp);
}

void log_real(vector<complex<double>> &phi, int f) {
	FILE * fp;
	//char fname[20] = "phi_";
	//fname[4] = 'a' + f;
	//strcpy(fname + 5, ".data");

	fp = fopen("phi_r.data","a");
	fprintf(fp, "#%d\n", f);
	for (int i=0; i<NX; i++) {
		fprintf(fp, "%.3e\n", real(phi[i]));
	}
	fprintf(fp, "\n\n");
	fclose(fp);
}


void log_complex(vector<complex<double>> &phi, int f) {
	FILE * fp;

	fp = fopen("phi_c.data","a");
	fprintf(fp, "#%d\n", f);
	for (int i=0; i<NX; i++) {
		fprintf(fp, "%.3e %.3e\n", real(phi[i]), imag(phi[i]));
	}
	fprintf(fp, "\n\n");
	fclose(fp);
}


// takes phi_n and returns phi_n+1  (n is time index)
// V is potential function
vector<complex<double>> step(vector<complex<double>> &phi, double (*V)(int)) {
	int n = phi.size();
	complex<double> li = {0.0, 1.0};
	phi[0] = phi[n-1] = 0.0; // b.c.
	
	vector<double> a(n-1); // upper/lower diagonal of matrix
	for (int i=0; i < n-1; i++) {
		a[i] = C1;
	}
	
	vector<complex<double>> b(n);
	for (int i=0; i<n; i++) {
		b[i] = li - 2.0 * C1 - 0.5 * V(i);
	}

	vector<complex<double>> d(n);
	for (int i=1; i < n-1; i++) {
		d[i] = -C1 * phi[i-1] + phi[i] * (2.0 * C1 + 0.5 * V(i)) - C1 * phi[i+1] + li * phi[i];
	}
	
	d[0] = -C1 * phi[1];
	d[n-1] = -C1 * phi[n-2];
	
	vector<complex<double>> next_phi = tridiag_solve(a,b,d); 

	next_phi[0] = next_phi[n-1] = 0.0;


	return next_phi;
	
// solve TX = D where T_jj = b, T_jj+1 = T_j+1j = a 
}


void run(vector<complex<double>> &phi, double (*V)(int)) {
	FILE * fp;
	fp = fopen("phi_n.data","w");
	fp = fopen("phi_c.data","w");
	fclose(fp);

	cout << totalsum(phi) << endl;
	vector<complex<double>> next_phi;
	
	next_phi = step(phi, V);

	for (int t=1; t < NT; t++) {
		cout << totalsum(next_phi) << endl;
		next_phi = step(next_phi, V);
		log_norm(next_phi, t);
//		log_real(next_phi, t);
		log_complex(next_phi, t);
	}
}


