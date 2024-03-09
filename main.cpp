#include "Matrix.h"
#include <iostream>
using namespace std;

#define N 4

int main(){
	Matrix A(N, N, 10);
	for (int i = 0; i < A.getRows(); ++i) {
		for (int j = 0; j < A.getColumns(); ++j) {
			A[i][j] = static_cast<double>(1) / (i + j + 2);
		}
	}
	cout << A;
	
	Matrix L, U;
	double wyznacznik_A = A.rozkladLU(L, U);
	cout << "Macierz A: \n" << A;
	cout << "Macierz L: \n" << L;
	cout << "Macierz U: \n" << U;
	cout << "Wyznacznik macierzy A, rowna sie " << wyznacznik_A << '.' << endl;

	// obliczenie macierzy A^(-1)
	// ustawianie wektorów x1, x2, x3, ...
	// X - to będzie macierz odwrotności
	Matrix X = A;
	X.setToUnitMatrix();

	for (int i = 0; i < A.getRows(); ++i) {
		std::vector<double>x_o; // wektor x_i
		for (int k = 0; k < A.getRows(); ++k) {
			x_o.push_back(0);
		}
		x_o[i] = 1;
		A.eliminacjaGaussaJordana(X[i], x_o);
	}
	cout << X;	// wypisuje macierz odwrotności

	// wskaźnik uwarunkowania macierzy
	// dla normy: ||A|| = max(a_ij) - oznacza to, największy wyraz w macierzy
	// wsk_uwarunk = ||A|| * ||A^(-1)||
	double wsk_uwarunkowania = A.getMaxElement() * X.getMaxElement();
	cout << "Wskaznik uwarunkowania macierzy A, rowna sie " << wsk_uwarunkowania << ".\n";

	return 0;
}
