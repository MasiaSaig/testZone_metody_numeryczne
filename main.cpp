/*

(h_11 h_12 ... )
(h_21 h_22 ... )
( ...   	   )

w_0=1
w_1=h_11 - lamb
w_2 = (h_22 - lamb)w_k-1 +(h_12)^2w_k-2

ile wartości wlasnych mniejszych ni lamb

życie metody bisekcji

liczba zmian znaków oznaza ilosc wartości mnijszych niż lamb
 jeśli wieksza jest liczba wartosci własnych niż liczba zmian znaków, to idziemy do lewego
jeśli mniejsza lub równa, to idzeimy do prawego



szukamy najwiekszej sumy || wiersza


metoda ?householdera? - transform macierzy do posatci trójdiagonalnej, zależnie czy jest symetralna
*/

#include "Matrix.h"
#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;


using std::abs;

int main() {
	const int N = 5, IT_MAX=50;
	double L=5;
	//double L = 10;

	double delta_x = ((2.0*L) / static_cast<double>(N));

	Matrix A(N,N);
	//A.zeroOut();
	for(int i=1; i<N; ++i){
		double val = -1.0 / (2.0*delta_x*delta_x); 
		A[i][i-1] = val;
		A[i-1][i] = val;
	}
	for(int i=0; i<N; ++i){
		A[i][i] = (1.0 / (delta_x * delta_x)) + ((i*delta_x - L)*(i*delta_x - L) / 2.0);
	}

	std::cout << A;

	const int K = 3;
	// szukanie maksymanej wartości, do stworzenia przedziału
	double S=A[0][0];
	for(int i=0; i<N; ++i){
		double suma = 0;
		for(int j=0; j<N; ++j){
			suma += abs(A[i][j]);
		}
		if(suma >= S) S=suma;
	}

	double l = -S;
	double p = S;
	double w[N]{};
	// i-ta wartosc wlasna
	double lambda = 0;	// (S+(-S)) / 2
	
	w[0] = 1; 
	w[1] = A[0][0] - lambda;
	
	double wartosci_wlasne[N]{};
	int zmiany_znaku;
	std::cout << std::endl << "Lambdy: ";
	for(int i=0; i<IT_MAX; ++i){
		zmiany_znaku = 0;
		for(int k=2; k<N; ++k){
			w[k] = (A[k-1][k-1]-lambda)*w[k-1] - (A[k-1][k])*(A[k-1][k])*w[k-2];
			if(w[k]*w[k-1] < 0) ++zmiany_znaku;	// dodaj 1, jeśli wystąpiła zmiana znaku
			if (w[k] == 0) wartosci_wlasne[i] = lambda;
		}
		if (zmiany_znaku <= i) l = lambda;
		else p = lambda;
		lambda = (p + l) / 2;
		std::cout << lambda << ", ";
	}
	std::cout << std::endl;
	std::cout << std::endl << "Wartosc wlasna(lambda): " << lambda << std::endl;
	// obliczanie (i-tego) wektora wlasnego z lambda
	double x[N]{};
	x[0] = 1;
	x[1] = (lambda - A[0][0]) / A[1][0];
	std::cout << "Wektor wlasny: " << x[0] << ',' << x[1] << ',';
	for(int i=2; i<N; ++i){
		x[i] = ((lambda - A[i][i])*x[i-1] - A[i-1][i-2]*x[i-2]) / A[i][i-1];
		std::cout << x[i] << ',';
	}
	std::cout << std::endl << std::endl;
	// unormowanie wektora
	std::cout << "Wektor wlasny unormowany: ";
	for (int i = 0; i < N; ++i) {
		x[i] /= static_cast<double>(N);
		std:cout << x[i] << ',';
	}

	return 0;
}
