#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

const int g = 10;
const double c = 0.2, m = 0.2;
# define PI           3.14159265358979323846

double * x_prime(double * x, double * v);
double * v_prime(double * x, double * v);
double * avg_k(double * k1, double * k2, double * k3, double * k4, double dt);

double * scalar_mul(double s, double * v, int N) {
	double * r = new double[N];
	for(int i = 0; i<N; i++) {
		r[i] = s*v[i];
	}
	return r;
}

double * add(double * v, double * w, int N) {
	double * r = new double[N];
	for(int i = 0; i<N; i++) {
		r[i] = v[i]+w[i];
	}
	return r;
}

double abs(double * v, int N) {
	double sum = 0;
	for(int i = 0; i<N; i++) {
		sum += v[i]*v[i];
	}
	return sqrt(sum);
}

double * append(double * v, double * w, int M, int N) {
	double * r = new double[M+N];
	int i;
	for(i = 0; i<M; i++) {
		r[i] = v[i];
	}
	for(i = 0; i<N; i++) {
		r[i+M] = w[i];
	}
	return r;
}

int main() {
	double dt = 0.001, v0 = 300, angle = 45;
	int N = 2000;
	double x_0[] = {0, 0};
	double v_0[] = {v0*cos(angle*PI/180), v0*sin(angle*PI/180)};
	double * x_now, * v_now, * x_aux, * v_aux, * x_x = new double[N+1], * x_y = new double[N+1];
	x_x[0] = x_0[0];
	x_y[0] = x_0[1];
	for(int i = 0; i<N; i++) {

	}
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

double * avg_k(double * k1, double * k2, double * k3, double * k4, double dt) {
	double * avg_k = new double[2];
	for(int i = 0; i<2; i++) {
		avg_k[i] = 1.0/6.0 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) * dt;
	}
	return avg_k;
}

double * RK4Step(double * x_old, double * v_old, double dt) {
	//k1
	double * k1_x = x_prime(x_old, v_old);
	double * k1_v = x_prime(x_old, v_old);

	//k2
	double * x_1 = add(x_old, scalar_mul(0.5*dt, k1_x, 2), 2);
	double * v_1 = add(v_old, scalar_mul(0.5*dt, k1_v, 2), 2);

	double * k2_x = x_prime(x_1, v_1);
	double * k2_v = v_prime(x_1, v_1);

	//k3
	double * x_2 = add(x_old, scalar_mul(0.5*dt, k2_x, 2), 2);
	double * v_2 = add(v_old, scalar_mul(0.5*dt, k2_v, 2), 2);

	double * k3_x = x_prime(x_2, v_2);
	double * k3_v = v_prime(x_2, v_2);

	//k4
	double * x_3 = add(x_old, scalar_mul(dt, k3_x, 2), 2);
	double * v_3 = add(v_old, scalar_mul(dt, k3_v, 2), 2);

	double * k4_x = x_prime(x_3, v_3);
	double * k4_v = v_prime(x_3, v_3);

	//1/6 * dt * (k1 + 2*k2 + 2*k3 + k4)
	double * avg_k_x = avg_k(k1_x, k2_x, k3_x, k4_x, dt);
	double * avg_k_v = avg_k(k1_v, k2_v, k3_v, k4_v, dt);

	double * x_new = add(x_old, avg_k_x, 2);
	double * v_new = add(v_old, avg_k_v, 2);
	return append(x_new, v_new, 2, 2);
}
