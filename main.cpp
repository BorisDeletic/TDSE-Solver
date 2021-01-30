// time dependent schrodinger equation
#include <iostream>
#include <stdio.h>
#include <math.h>     
#include <complex>
#include <cmath>
#include <vector>
#include "matrix.h"
#include "cn_solver.h"

#define N 128


using namespace std;

double squared(int x);
double gaussian(int x);
void initialise(vector<complex<double> > &phi, double (*f)(int));
void initialise(vector<double> &phi, double (*f)(int));
vector<double> normalise(vector<complex<double> > &phi);


double squared(int x) {return (double)(x*x);}
double gaussian(int x) {
	return (double)(exp(-x*x / (float)(2*N)));
}


vector<double> normalise(vector<complex<double> > &phi) {
	vector<double> normalised(N);
	for (int i=0; i < N; i++) {
		normalised[i] = norm(phi[i]);
	}
	return normalised;
}


void initialise(vector<complex<double> > &phi, double (*f)(int)) {
	for (int i = 0; i < N; i++) {
		phi[i] = {f(i), -f(i - (float)(N/2))}; //out of phase
		//R, I
	}
	phi[0] = {0, 0};
	phi[N-1] = {0, 0};
	//b.c.
}


void initialise(vector<double> &phi, double (*f)(int)) {
	for (int i = 0; i < N; i++) {
		phi[i] = f(i); 
	}
	phi[0] = 0;
	phi[N-1] = 0;
	//b.c.
}



int main() {
	vector<complex<double> > phi(N);
	vector<double> V(N);
	initialise(phi, gaussian);
	initialise(V, squared);
	
	FILE * fp;
	FILE * fp2;
	fp = fopen ("complex.data","w");
	fp2 = fopen ("norm.data","w");
	fprintf(fp, "#i x y\n");
	fprintf(fp2, "#i phi^2\n");
	
	for (int i=0; i<N; i++) {
		fprintf(fp, "%d %.3e %.3e\n", i, real(phi[i]), imag(phi[i]));
	}
	
	vector<double> normd = normalise(phi);
	for (int i=0; i<N; i++) {
		fprintf(fp2, "%d %.3e\n", i, normd[i]);
	}
	
	fclose(fp);
	fclose(fp2);
	
}

///set term qt persist size 700,500