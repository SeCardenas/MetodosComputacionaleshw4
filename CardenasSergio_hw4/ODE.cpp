#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

const int g = 10;
const double c = 0.2, m = 0.2;

vector<double> derivativeX(vector<double> x, vector<double> v);
vector<double> derivativeV(vector<double> x, vector<double> v);

int main() {
	double dt = 0.01;
}

vector<double> derivativeX(vector<double> x, vector<double> v) {
	vector<double> deriv(2);
	deriv[0] = v[0];
	deriv[1] = v[1];
	return deriv;
}

vector<double> derivativeV(vector<double> x, vector<double> v) {
	vector<double> deriv(2);
	deriv[0] = -c/m * sqrt((v[0]*v[0] + v[1]*v[1]))*v[0]; //dv/dt en x
	deriv[1] = -g-c/m * sqrt((v[0]*v[0] + v[1]*v[1]))*v[1]; //dv/dt en y
	return deriv;
}
