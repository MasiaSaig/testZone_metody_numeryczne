#include "Matrix.h"
#include <iostream>
#include <iomanip>
#include <cassert>

using std::cout;
using std::endl;
using std::setw;

//Matrix::~Matrix(){}
Matrix::Matrix(int rows, int cols, int precision) :
	_rows{rows}, _columns{cols}, _output_precision{precision}
{
    // uzupe�nia macierz jedynkami
    for (int i = 0; i < _rows; ++i) {
        std::vector<double> row;
        for (int j = 0; j < _columns; ++j) {
            row.push_back(1);
        }
        _matrix.push_back(row);
    }
}
Matrix::Matrix(Matrix& matrix) : 
    _rows{matrix._rows}, 
    _columns{matrix._columns}, 
    _output_precision{matrix._output_precision}
{
    _matrix.insert(_matrix.begin(), matrix._matrix.begin(), matrix._matrix.end());
}

Matrix& Matrix::operator= (Matrix& obj) {
    if (this != &obj) {
        _rows = obj._rows;
        _columns = obj._columns;
        _matrix.insert(_matrix.begin(), obj._matrix.begin(), obj._matrix.end());
    }
	return *this;
}

Matrix Matrix::operator+(const Matrix& matrix_a)
{
    assert((this->_rows == matrix_a._rows && this->_columns == matrix_a._columns) && "Rozmiary macierzy si� niezgadzaj�!");
    Matrix temp(_rows, _columns);
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _columns; ++j) {
            _matrix[i][j] = this->_matrix[i][j] + matrix_a._matrix[i][j];
        }
    }
    return temp;
}
Matrix Matrix::operator-(const Matrix& matrix_a)
{
    assert((this->_rows == matrix_a._rows && this->_columns == matrix_a._columns) && "Rozmiary macierzy si� niezgadzaj�!");
    Matrix temp(_rows, _columns);
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _columns; ++j) {
            _matrix[i][j] = this->_matrix[i][j] - matrix_a._matrix[i][j];
        }
    }
    return temp;
}
Matrix Matrix::operator*(const Matrix& matrix_a)
{
    assert((this->_columns == matrix_a._rows) && "Rozmiary macierzy si� niezgadzaj�!");
    Matrix temp(_rows, _columns);
    for (int i = 0; i < _rows; ++i) {
        
        for (int j = 0; j < _columns; ++j) {
            double sum = 0;
            for (int k = 0; k < _columns; ++k) {
                sum += this->_matrix[i][k] * matrix_a._matrix[k][j];
            }
            temp._matrix[i][j] = sum;
        }
    }
    return temp;
}
std::vector<double>& Matrix::operator[] (int idx) {
    return _matrix.at(idx);
}
const std::vector<double>& Matrix::operator[] (int idx) const {
    return _matrix.at(idx);
}

double Matrix::getRows() const { return _rows; }
double Matrix::getColumns() const { return _columns; }

void Matrix::zeroOut() {
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _columns; ++j) {
            _matrix[i][j] = 0;
        }
    }
}
void Matrix::clear() {
    for (int i = 0; i < _rows; ++i) {
        _matrix[i].clear();
    }
    _matrix.clear();
}

void Matrix::eliminacjaGaussaJordana(double c[], double y[]) {
    double l, w;
    for (int k = 0; k < _rows; ++k) {
        l = _matrix[k][k];
        // podzielenie wiersza przez warto�� na diagonali
        for (int j = k; j < _columns; ++j) {
            _matrix[k][j] /= l;
        }
        y[k] /= l;
        // odejmowanie wierszy
        for (int i = 0; i < _rows; ++i) {
            if (i == k)
                continue;
            w = _matrix[i][k];
            for (int j = k; j < _columns; ++j) {
                _matrix[i][j] -= _matrix[k][j] * w;
            }
            y[i] = y[i] - y[k] * w;
        }
    }
    for (int i = 0; i < _rows; ++i) {
        c[i] = y[i];
    }
    showSimpleEquation(*this, y);
}
double Matrix::rozkladLU(Matrix& l, Matrix& u) const {
    Matrix L(_rows, _columns);
    Matrix U(_rows, _columns);
    L.zeroOut();
    U.zeroOut();
    for (int i = 0; i < _rows; ++i) { L[i][i] = 1; }
    int i, j, k;
    double sum;
    for (j = 0; j < _columns; j++)
    {
        for (i = 0; i <= j; i++)
        {
            sum = 0;
            for (k = 0; k < i; k++) sum += L[i][k] * U[k][j];
            U[i][j] = _matrix[i][j] - sum;
        }
        for (i = j + 1; i < _rows; i++)
        {
            sum = 0;
            for (k = 0; k < j; k++) sum += L[i][k] * U[k][j];
            if (U[j][j]) L[i][j] = (_matrix[i][j] - sum) / U[j][j];
        }
    }
    l = L;
    u = U;
    double determinant = 1;
    for (i = 0; i < _rows; ++i) {
        determinant *= u[i][i];
    }
    return determinant;
}

void Matrix::showSimpleEquation(const Matrix& matrix_a, const Matrix& matrix_b, char sign, const Matrix& matrix_result) {
    int rows = matrix_a._rows;
    int cols_a = matrix_a._columns, cols_b = matrix_b._columns, cols_res = matrix_result._columns;
    for (int i = 0; i < rows; ++i) {
        cout << "| ";
        for (int j = 0; j < cols_a; ++j) {
            cout << setw(matrix_a._output_precision) << matrix_a._matrix[i][j] << ' ';
        }
        if (i == (rows / 2))
            cout << '|' << sign;
        else
            cout << "| ";

        for (int j = 0; j < cols_b; ++j) {
            cout << setw(matrix_b._output_precision) << matrix_b._matrix[i][j] << ' ';
        }
        if (i == (rows / 2))
            cout << "|=";
        else
            cout << "| ";

        for (int j = 0; j < cols_res; ++j) {
            cout << setw(matrix_result._output_precision) << matrix_result._matrix[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
}

void Matrix::showSimpleEquation(const Matrix& matrix_a, const double y[], const double c[]){
    int rows = matrix_a._rows, cols = matrix_a._columns;
    for (int i = 0; i < rows; ++i) {
        cout << "| ";
        for (int j = 0; j < cols; ++j) {
            cout << setw(matrix_a._output_precision) << matrix_a._matrix[i][j] << ' ';
        }
        cout << "|";
        if (c) {
            cout << setw(matrix_a._output_precision) << c[i] << ' ';
        }else {
            cout << " c" << i << ' ';
        }
        cout << "|=";
        cout << setw(matrix_a._output_precision) << y[i] << ' ' << endl;
    }
    cout << endl;
}

std::ostream& operator<< (std::ostream& stream, const Matrix& matrix) {
    for (int i = 0; i < matrix._rows; ++i) {
        stream << "|";
        for (int j = 0; j < matrix._columns; ++j) {
            stream << setw(matrix._output_precision) << matrix._matrix[i][j] << ' ';
        }
        stream << "|" << endl;
    }
    std::cout << endl;
    return stream;
}
