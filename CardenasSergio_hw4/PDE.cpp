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
const double alpha = 0.5;

double get_r(int i, int j, int N, double dx);
double initial_cond(int i, int j, int N, double dx);

using namespace std;

int main() {
	int i, j, k;
	double dx = L/(N_x-1);
	double dt = alpha*dx/nu*dx;
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
			cout << T_0[i][j] << " ";
		}
		cout << endl;
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
