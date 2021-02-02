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
double quartic(int x);
double gaussian(int x);
complex<double> moving(int x);
complex<double> collision(int x);
double zero(int x);
void initialise(vector<complex<double>> &phi, complex<double> (*f)(int));
void initialise(vector<complex<double> > &phi, double (*f)(int));
void initialise(vector<double> &phi, double (*f)(int));


double quartic(int x) {
	if (x == 0 || x == NX-1) return 1000000.0;
	double t = (double)x;
	double n = (double)NX;	
	return (t / n - 0.5) * (t / n - 0.5) *(t / n - 0.5) *(t / n - 0.5) * 30;
}
double squared(int x) {
	if (x == 0 || x == NX-1) return 1000000.0;
	double t = (double)x;
	double n = (double)NX;	
	return (t / n - 0.5) * (t / n - 0.5) * 30;
}
double gaussian(int x) {
	return (double)(exp(-x*x / (double)(NX)));
}
double zero(int x) {
	if (x == 0 || x == NX-1) 
		return 100000.0;
	return 0.0;
}
double inverse(int x) {
	if (x == NX/2) return -100000.0;
	double p = (double) -1.0 / abs(x - NX/2);
	return p;
}
double wall(int x) {
	if (x == 0 || x == NX-1) return 100000;
	if (x > NX/2-500 && x < NX/2+1000) return 2;
	return 0;
}


complex<double> moving(int x) {
	double xpos = 500;
	double sd = NX * 20.0;
	double freq = 20.0;
	complex<double> li = {0.0, 1.0};
	complex<double> z = exp(li * (complex<double>)(x / freq)) * exp(-(x-xpos)*(x - xpos) / sd);
	return z;
}
complex<double> collision(int x) {
	double xpos = 3500.0;
	double sd = NX * 4.0;
	double freq = 100.0;
	complex<double> li = {0.0, 1.0};
	complex<double> z = exp(li * (complex<double>)(x / freq)) * exp(-(x-xpos)*(x - xpos) / sd);
	complex<double> z2 = exp(li * (complex<double>)(-x / freq)) * exp(-(x+xpos-NX)*(x + xpos-NX) / sd);
	return z + z2;
}


void initialise(vector<complex<double>> &phi, complex<double> (*f)(int)) {
	for (int i = 0; i < NX; i++) {
		phi[i] = f(i); //out of phase
		//R, I
	}
	phi[0] = {0, 0};
	phi[NX-1] = {0, 0};
	//b.c.
}


void initialise(vector<complex<double> > &phi, double (*f)(int)) {
	for (int i = 0; i < NX; i++) {
		phi[i] = {f(i - 100) + f(i - NX + 100), f(i-150)}; //out of phase
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
	vector<complex<double>> phi(NX);
	
	initialise(phi, moving);
	double (*V)(int) = wall;
	
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
		fprintf(fp2, "%.3e\n", V(i));
	}

	vector<double> normd = normalise(phi);
	for (int i=0; i<NX; i++) {
		fprintf(fp3, "%.3e\n", normd[i]);
	}

	
	fclose(fp);
	fclose(fp2);
	fclose(fp3);
	
	norm_wavefunc(phi);
	run(phi, V);
	
}

///set term qt persist size 700,500
// plot 'complex.data' using 1:2, 'complex.data' using 1:3