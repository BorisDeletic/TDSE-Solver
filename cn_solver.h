#ifndef CN_SOLVER_H
#define CN_SOLVER_H

const double dt = 10;
const double dx = 0.1;

#define NX 10000
#define NT 300

vector<complex<double>> step(vector<complex<double>> &phi, double (*V)(int));
void run(vector<complex<double>> &phi, double (*V)(int));

vector<double> normalise(vector<complex<double> > &phi);
void norm_wavefunc(vector<complex<double>> &phi);
double totalsum(vector<complex<double>> &phi);

void log_norm(vector<complex<double>> &phi, int f = 0);
void log_complex(vector<complex<double>> &phi, int f = 0);
void log_real(vector<complex<double>> &phi, int f = 0);


#endif