#ifndef _MY_MATRIX_NUMERICAL_METHODS_H_
#define _MY_MATRIX_NUMERICAL_METHODS_H_
#pragma once

#include <vector>
#include <iostream>

class Matrix {
public:
	//~Matrix();
	Matrix(int rows=0, int cols=0,  int precision=5);
	Matrix(const Matrix& matrix);
	template <size_t r, size_t c>
	Matrix(const double (&array)[r][c]) : Matrix(r,c) {
		*this = array;
	}

	template <size_t r, size_t c>
	Matrix& operator= (const double (&array)[r][c]) {
		for (size_t i = 0; i < r; ++i) {
			for (size_t j = 0; j < c; ++j) {
				_matrix[i][j] = array[i][j];
			}
		}
		return *this;
	}
	Matrix& operator= (const Matrix& matrix_a);
	Matrix operator+ (const Matrix& matrix_a);
	Matrix operator- (const Matrix& matrix_a);
	Matrix operator* (const Matrix& matrix_a);

	std::vector<double>& operator[] (int idx);
	const std::vector<double>& operator[] (int idx) const;
	
	inline const double getRows() const { return _rows; }	//nie trzeba dodawaæ inline
	inline const double getColumns() const { return _columns; }
	const double getMaxElement() const;

	void zeroOut();
	void setToUnitMatrix();
	void clear();

	void eliminacjaGaussaJordana(std::vector<double>& c, std::vector<double> y) const;
	double rozkladLU(Matrix& L, Matrix& U) const; // zwraca wyznacznik macierzy

	static void showSimpleEquation(const Matrix& matrix_a, const Matrix& matrix_b, char sign, const Matrix& matrix_result);
	static void showSimpleEquation(const Matrix& matrix_a, const double y[], const double c[]=nullptr);
	friend std::ostream& operator<< (std::ostream& stream, const Matrix& matrix);
private:
	int _rows, _columns;
	std::vector<std::vector<double>> _matrix;

	int _output_precision;
};

std::ostream& operator<<(std::ostream& stream, const Matrix& matrix);

#endif
