#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

const double c = 0.2, m = 0.2, g = 10;
# define PI           3.14159265358979323846

double * x_prime(double * x, double * v);
double * v_prime(double * x, double * v);
double * avg_k(double * k1, double * k2, double * k3, double * k4, double dt);
double * RK4_step(double * x_old, double * v_old, double dt);
double recorrido(double dt, double v0, double angle, ofstream * output);

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
	//parte 1
	double dt = 0.001, v0 = 300, angle = 45;
	ofstream output;
  	output.open ("proyectil.txt");

  	recorrido(dt, v0, angle, &output);
	output.close();

	//parte 2
	double angles[] = {10, 20, 30, 40, 50, 60, 70};
	int num_angles = sizeof(angles)/sizeof(*angles);
	double distance, max_distance = 0, angle_max_d = 0;

	output.open("proyectil2.txt");
	cout << "ángulos:" << endl;
	for(int i = 0; i<num_angles; i++) {
		cout << angles[i] << "°" << endl;
		distance = recorrido(dt, v0, angles[i], &output);
		if(distance > max_distance) {
			max_distance = distance;
			angle_max_d = angles[i];
		}
	}
	cout << "El ángulo con el que se consigue la mayor distancia recorrida es " << angle_max_d << "° con una distancia total de " << max_distance << endl;
	output.close();
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

double * RK4_step(double * x_old, double * v_old, double dt) {
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

	//se retorna x_new y v_new en un solo arreglo
	double * x_new = add(x_old, avg_k_x, 2);
	double * v_new = add(v_old, avg_k_v, 2);
	return append(x_new, v_new, 2, 2);
}

double recorrido(double dt, double v0, double angle, ofstream * output) {
	double x_0[] = {0, 0};
	double v_0[] = {v0*cos(angle*PI/180), v0*sin(angle*PI/180)};
	double * next, t = 0;

  	(*output) << t << " " << x_0[0] << " " << x_0[1] << " " << v_0[0] << " " << v_0[1] << endl;

	while(x_0[1] >= 0) {
		next = RK4_step(x_0, v_0, dt);
		x_0[0] = next[0];
		x_0[1] = next[1];
		v_0[0] = next[2];
		v_0[1] = next[3];
		t += dt;
		//guardar recorrido en un archivo
		(*output) << t << " " << x_0[0] << " " << x_0[1] << " " << v_0[0] << " " << v_0[1] << endl;

	}
	cout << "distancia recorrida con un ángulo inicial de " << angle << "°: " << x_0[0] << "m" << endl;
	return x_0[0];
}
