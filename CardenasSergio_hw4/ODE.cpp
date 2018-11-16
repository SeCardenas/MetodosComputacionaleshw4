#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

const int g = 10;
const double c = 0.2, m = 0.2;

double * x_prime(double * x, double * v);
double * v_prime(double * x, double * v);

void scalar_mul(double s, double * v, int N) {
	for(int i = 0; i<N; i++) {
		v[i] = s*v[i];
	}
}

double abs(double * v, int N) {
	double sum = 0;
	for(int i = 0; i<N; i++) {
		sum += v[i]*v[i];
	}
	return sqrt(sum);
}

int main() {
	double dt = 0.01;
}

/**
 * devuelve la derivada de la posición en términos de la posición actual y la velocidad actual.
 */
double * x_prime(double * x, double * v) {
	double * deriv = new double[2];
	deriv[0] = v[0];
	deriv[1] = v[1];
	return deriv;
}

/**
 * devuelve la derivada de la velocidad en términos de la posición actual y la velocidad actual.
 */
double * v_prime(double * x, double * v) {
	double * deriv = new double[2];
	deriv[0] = -c/m * abs(v, 2)*v[0]; // dv/dt en x
	deriv[1] = -g-c/m * abs(v, 2)*v[1]; // dv/dt en y
	return deriv;
}

double * RK4Step(double * x_old, double * v_old, double dt) {
	//k1
	double * k1_x = x_prime(x_old, v_old);
	double * k1_v = x_prime(x_old, v_old);

}
