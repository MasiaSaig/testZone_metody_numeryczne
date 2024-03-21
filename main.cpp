/*
sprawozdanie - sprzężonych gradientów i największego spadku
ćwiczenia - metoda największego spadku

wektor reszty: Ax_k - b = 0;

funkcja gradientu
Q = 1/2 * x^T * A*x - x^T * b
grad(Q) = Ax - b

x_k+1 = x_k + a*v_k
v_k = grad(Q) = -r_k - wektor reszt
r_k = b-Ax_k

x_k+1 = x_k + a

a = (r^T * r) / (r^T * Ar) - alpha mowi jak daleko idziemy w stronę gradientu

do{}
while((||r_k+1||-||r_k||) > tolerancja)

wykres reszt zalezny od iteracji




metoda bisekcji - pozwala okreslone własności własne
*/

#include "Matrix.h"
#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

using std::abs;
using std::sqrt;

int main() {
	const int n = 1000, m = 5;

	for (int i = 0; i < 131; ++i) { cout << i << ','; }

	cout << endl << endl;

	Matrix A(n, n);
	for (int i = 0; i < A.getRows(); ++i) {
		for (int j = 0; j < A.getColumns(); ++j) {
			if (abs(i - j) <= m) {
				A[i][j] = 1.0 / (1.0 + abs(i - j));
			}
			else {
				A[i][j] = 0;
			}
		}
	}

	Matrix b(n, 1);
	for (int i = 0; i < n; ++i) { b[i][0] = i; }
	//cout << "b:" << b << endl;

// Metoda najwiekszego spadku x = 0
	Matrix x(n, 1);
	Matrix r(n, 1);
	double a, war_normy_eukl_r, war_normy_eukl_x;
	for (int i = 0; i < n; ++i) { x[i][0] = 0; }
	//cout << "x:" << x << endl;


	//Matrix B({{1,2}, {3,4}});
	//cout << "B:" << B << endl;
	int iteracje = 0;
	do {
		war_normy_eukl_r = 0;
		r = b - (A * x);
		//cout << "r:" << r << endl;
		//Matrix licz = (r.transpose() * r);
		//Matrix mian = (r.transpose() * A * r);
		//cout << "trans(r):" << endl << trans(r);
		// wybieramy element [0][0] ponieważ wynik mnożenia tych macierzy 
		// jest macierzą o jednym elemencie, czyli tak naprawdę liczbą
		a = (r.transpose() * r)[0][0] / (r.transpose() * A * r)[0][0];
		//cout << "alpha = " << a << endl;
		x = x + (r * a);
		//cout << "x:" << x << endl;
		++iteracje;
		war_normy_eukl_r = sqrt(abs((r.transpose() * r)[0][0]));
		war_normy_eukl_x = sqrt(abs((x.transpose() * x)[0][0]));
		 cout << iteracje << ";" << war_normy_eukl_r << ";" << a << 
			 ";" << war_normy_eukl_x << endl;
	} while (abs(war_normy_eukl_r) > 0.000001);

	//cout << endl << "Iteracje: " << iteracje << " dla x = 0" << endl << endl;


// Metoda najwiekszego spadku x = 1
	for (int i = 0; i < n; ++i) { x[i][0] = 1; }

	iteracje = 0;
	do {
		r = b - (A * x);
		//Matrix licz = (r.transpose() * r);
		//Matrix mian = (r.transpose() * A * r);
		a = (r.transpose() * r)[0][0] / (r.transpose() * A * r)[0][0];
		x = x + (r * a);
		++iteracje;
		war_normy_eukl_r = sqrt(abs((r.transpose() * r)[0][0]));
		war_normy_eukl_x = sqrt(abs((x.transpose() * x)[0][0]));
		cout << iteracje << ";" << war_normy_eukl_r << ";" << a <<
			";" << war_normy_eukl_x << endl;
	} while (abs(war_normy_eukl_r) > 0.000001);

	cout << endl << "Iteracje: " << iteracje << " dla x = 1" << endl;


// Metoda sprzezonych gradientow dla macierzy wstegowej
	cout << "Metoda sprzezonych gradientow dla macierzy wstegowej" << endl;
	for (int i = 0; i < n; ++i) { b[i][0] = i+1; }
	for (int i = 0; i < n; ++i) { x[i][0] = 0; }
	Matrix v(n, 1);

	double alpha, beta;

	//v = r = b - (A * x);
	v = b;
	r = b;
	iteracje = 0;
	while ((r.transpose() * r)[0][0] > 0.000001) {
		alpha = (r.transpose() * r)[0][0] / (v.transpose() * A * v)[0][0];
		x = x + (v * alpha);
		Matrix temp_r = r;
		r = r - ((A * v) * alpha);
		beta = (r.transpose() * r)[0][0] / (temp_r.transpose() * temp_r)[0][0];
		v = r + (v * beta);
		++iteracje;
		cout << iteracje << ";" << sqrt(abs((r.transpose() * r)[0][0])) << ";"
			<< alpha << ";" << beta << ";" << sqrt(abs((x.transpose() * x)[0][0])) << endl;
	}
	cout << endl << "Iteracje: " << iteracje << " dla x = 0" << endl;
	
	return 0;
}
