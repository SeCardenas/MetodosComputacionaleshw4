#include <iostream>
#include <cmath>
#include <fstream>

const double k = 1.62;
const double Cp = 820;
const double rho = 2.71;
const double nu = k/(Cp*rho);

const double d = 0.1; // en metros
const double L = 0.5; // en metros
const double T0_c = 10, T0_b = 50; //°C

const int N_x = 51;
const double eta = 0.5;

double get_r(int i, int j, int N, double dx);
double initial_cond(int i, int j, int N, double dx);
void step_fixed_borders(double **T_old, double **T_new, double eta, int N);
void step_periodic_borders(double **T_old, double **T_new, double eta, int N);
void step_open_borders(double **T_old, double **T_new, double eta, int N);

using namespace std;

int main() {
	int i, j, k;
	double dx = L/(N_x-1);
	double dt = eta*dx/nu*dx;
	//Creación de las matrices como un arreglo de apuntadores
	double **T_0 = new double *[N_x], **T_1 = new double *[N_x], **T_aux;
	for(i = 0; i<N_x; i++) {
		T_0[i] = new double[N_x];
		T_1[i] = new double[N_x];
	}
	//Condiciones iniciales
	for(i = 0; i<N_x; i++) {
		for(j = 0; j<N_x; j++) {
			T_0[i][j] = initial_cond(i, j, N_x, dx);
			T_1[i][j] = T_0[i][j];
		}
	}
}

double get_r(int i, int j, int N, double dx) {
	double center = 0.5*(N-1);
	double x = (i-center)*dx, y = (j-center)*dx;

	return sqrt(x*x + y*y);
}

double initial_cond(int i, int j, int N, double dx) {
	double r = get_r(i, j, N, dx);
	if(r>d) {
		return 10; //°C
	}
	return 100; //°C
}

void step_fixed_borders(double **T_old, double **T_new, double eta, int N) {
	for(int i = 0; i<N; i++) {
		for(int j = 0; j<N; j++) {
			if(i!=0 && i!=N-1 && j!=0 && j!=N-1)
				T_new[i][j] = T_old[i][j] + eta*(T_old[i+1][j] + T_old[i-1][j] + T_old[i][j+1] + T_old[i][j-1] - 4*T_old[i][j]);
		}
	}
}

void step_periodic_borders(double **T_old, double **T_new, double eta, int N) {
	for(int i = 0; i<N; i++) {
		for(int j = 1; j<N-1; j++) {
			T_new[i][j] = T_old[i][j] + eta*(T_old[(i+1)%N][j] + T_old[(i+N-1)%N][j] + T_old[i][(j+1)%N] + T_old[i][(j+N-1)%N] - 4*T_old[i][j]);
		}
	}
}

void step_open_borders(double **T_old, double **T_new, double eta, int N) {
	for(int i = 0; i<N; i++) {
		for(int j = 1; j<N-1; j++) {
			if(i!=0 && i!=N-1 && j!=0 && j!=N-1)
				T_new[i][j] = T_old[i][j] + eta*(T_old[i+1][j] + T_old[i-1][j] + T_old[i][j+1] + T_old[i][j-1] - 4*T_old[i][j]);
		}
	}
}
