// time dependent schrodinger equation
#include <iostream>
#include <stdio.h>
#include <math.h>     
#include <complex>
#include <cmath>
#include <vector>
#include "matrix.h"
#include "cn_solver.h"


using namespace std;

double squared(int x);
double gaussian(int x);
double zero(int x);
void initialise(vector<complex<double> > &phi, double (*f)(int));
void initialise(vector<double> &phi, double (*f)(int));


double squared(int x) {
	if (x == 0 || x == NX-1) return 100000.0;
	double t = (double)x;
	double n = (double)NX;	
	return (t / n - 0.5) * (t / n - 0.5) * 5;
	}
double gaussian(int x) {
	return (double)(exp(-x*x / (double)(NX*NX/100)));
}
double zero(int x) {
	return 0;
}




void initialise(vector<complex<double> > &phi, double (*f)(int)) {
	for (int i = 0; i < NX; i++) {
		phi[i] = {f(i - (float)(NX/2)), 0.0}; //out of phase
		//R, I
	}
	phi[0] = {0, 0};
	phi[NX-1] = {0, 0};
	//b.c.
}


void initialise(vector<double> &phi, double (*f)(int)) {
	for (int i = 0; i < NX; i++) {
		phi[i] = f(i); 
	}
	phi[0] = 0;
	phi[NX-1] = 0;
	//b.c.
}



int main() {
	vector<complex<double> > phi(NX);
	vector<double> V(NX);
	initialise(phi, gaussian);
	initialise(V, squared);
	
	FILE * fp;
	FILE * fp2;
	FILE * fp3;
	fp = fopen ("complex.data","w");
	fp2 = fopen ("potential.data","w");
	fp3 = fopen ("norm.data","w");
	fprintf(fp, "#i x y\n");
	fprintf(fp2, "#V\n");
	fprintf(fp2, "#norm\n");
	
	for (int i=0; i<NX; i++) {
		fprintf(fp, "%d %.3e %.3e\n", i, real(phi[i]), imag(phi[i]));
	}
	
	for (int i=0; i<NX; i++) {
		fprintf(fp2, "%.3e\n", squared(i));
	}

	vector<double> normd = normalise(phi);
	for (int i=0; i<NX; i++) {
		fprintf(fp3, "%.3e\n", normd[i]);
	}

	
	fclose(fp);
	fclose(fp2);
	fclose(fp3);
	
	norm_wavefunc(phi);
	run(phi, squared);
	
}

///set term qt persist size 700,500
// plot 'complex.data' using 1:2, 'complex.data' using 1:3