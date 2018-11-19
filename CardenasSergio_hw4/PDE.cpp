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
const double eta = 0.2;

double get_r(int i, int j, int N, double dx);
double initial_cond(int i, int j, int N, double dx);
double step_fixed_borders(double **T_old, double **T_new, double eta, double dx, int N);
double step_periodic_borders(double **T_old, double **T_new, double eta, double dx, int N);
double step_open_borders(double **T_old, double **T_new, double eta, double dx, int N);

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
	double mean_old = 10, mean_new;
	int t = 0;
	while(true) {
		t++;
		//avance
		mean_new = step_open_borders(T_0, T_1, eta, dx, N_x);
		cout << ""; //Para que no se trabe
		//cambio de apuntadores
		T_aux = T_1;
		T_1 = T_0;
		T_0 = T_aux;
		//condición de salida (equilibrio)
		if(abs(mean_old - mean_new) <= 0.00000002) break; //Para que no se quede infinitamente en el loop, la diferencia no es exactamente 0
		mean_old = mean_new;
	}
	for(i = 0; i<N_x; i++) {
		for(j = 0; j<N_x; j++) {
			cout << T_aux[i][j] << " ";
		}
		cout << endl;
	}
	cout << t << endl;
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

double step_fixed_borders(double **T_old, double **T_new, double eta, double dx, int N) {
	double r, sum = 0;
	int numPoints = 0;
	for(int i = 0; i<N; i++) {
		for(int j = 0; j<N; j++) {
			r = get_r(i, j, N, dx);
			if(i!=0 && i!=N-1 && j!=0 && j!=N-1 && r>d) {
				T_new[i][j] = T_old[i][j] + eta*(T_old[i+1][j] + T_old[i-1][j] + T_old[i][j+1] + T_old[i][j-1] - 4*T_old[i][j]);
			}
			//cálculo del promedio
			if(r>d) {
				sum += T_new[i][j];
				numPoints++;
			}
		}
	}
	return sum/numPoints;
}

double step_periodic_borders(double **T_old, double **T_new, double eta, double dx, int N) {
	double r, sum = 0;
	int numPoints = 0;
	for(int i = 0; i<N; i++) {
		for(int j = 0; j<N; j++) {
			r = get_r(i, j, N, dx);
			if(r>d) {
				T_new[i][j] = T_old[i][j] + eta*(T_old[(i+1)%N][j] + T_old[(i+N-1)%N][j] + T_old[i][(j+1)%N] + T_old[i][(j+N-1)%N] - 4*T_old[i][j]);
				//para el cálculo del promedio
				sum += T_new[i][j];
				numPoints++;
			}
		}
	}
	return sum/numPoints;
}

double step_open_borders(double **T_old, double **T_new, double eta, double dx, int N) {
	double r, sum, temp = 0;
	int numPoints = 0;
	double deriv_x, deriv_y;
	for(int i = 0; i<N; i++) {
		for(int j = 0; j<N; j++) {
			r = get_r(i, j, N, dx);

			if(r>d) {
				//Si es un borde, se usa segunda derivada de adelanto/atraso en x/y según sea el caso
				if(i == 0) {
					deriv_x = T_old[i+2][j] + T_old[i][j] - 2*T_old[i+1][j];
				}
				else if(i == N-1) {
					deriv_x = T_old[i][j] + T_old[i-2][j] - 2*T_old[i-1][j];
				}
				else {
					deriv_x = T_old[i+1][j] + T_old[i-1][j] - 2*T_old[i][j];
				}
				if(j == 0) {
					deriv_y = T_old[i][j+2] + T_old[i][j] - 2*T_old[i][j+1];
				}
				else if(j == N-1) {
					deriv_y = T_old[i][j] + T_old[i][j-2] - 2*T_old[i][j-1];
				}
				else {
					deriv_y = T_old[i][j+1] + T_old[i][j-1] - 2*T_old[i][j];
				}
				T_new[i][j] = T_old[i][j] + eta*(deriv_x + deriv_y);

				//se suma para el promedio
				sum += T_new[i][j];
				numPoints++;
			}
		}
	}
	return sum/numPoints;
}
