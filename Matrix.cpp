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
Matrix::Matrix(const Matrix& obj) : 
    _rows{ obj._rows}, _columns{ obj._columns}, _output_precision{ obj._output_precision}
{
    _matrix.insert(_matrix.begin(), obj._matrix.begin(), obj._matrix.end());
}
Matrix::Matrix(Matrix&& obj) :
    _rows{ std::move(obj._rows) }, _columns{ std::move(obj._rows) }, _output_precision{ std::move(obj._output_precision) }
{
    _matrix = std::move(obj._matrix);
    // for (int i = 0; i < _rows; ++i) { _matrix[i] = std::move(obj._matrix[i]); }
}

Matrix& Matrix::operator= (const Matrix& obj) {
    if (this == &obj) 
        return *this;
    _rows = obj._rows;
    _columns = obj._columns;
    _output_precision = obj._output_precision;
    _matrix.insert(_matrix.begin(), obj._matrix.begin(), obj._matrix.end());
}
Matrix& Matrix::operator= (Matrix&& obj) {
    if (this == &obj)
        return *this;
    _rows = std::move(obj._rows);
    _columns = std::move(obj._columns);
    _output_precision = std::move(obj._output_precision);
    _matrix = std::move(obj._matrix);
}
Matrix Matrix::operator+(const Matrix& obj)
{
    assert((this->_rows == obj._rows && this->_columns == obj._columns) && "Rozmiary macierzy sie niezgadzaja!");
    Matrix temp(_rows, _columns);
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _columns; ++j) {
            _matrix[i][j] = this->_matrix[i][j] + obj._matrix[i][j];
        }
    }
    return temp;
}
Matrix Matrix::operator-(const Matrix& obj)
{
    assert((this->_rows == obj._rows && this->_columns == obj._columns) && "Rozmiary macierzy si� niezgadzaj�!");
    Matrix temp(_rows, _columns);
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _columns; ++j) {
            _matrix[i][j] = this->_matrix[i][j] - obj._matrix[i][j];
        }
    }
    return temp;
}
Matrix Matrix::operator*(const Matrix& obj)
{
    assert((this->_columns == obj._rows) && "Rozmiary macierzy si� niezgadzaj�!");
    Matrix temp(_rows, _columns, _output_precision+4);
    for (int i = 0; i < _rows; ++i) {
        
        for (int j = 0; j < _columns; ++j) {
            double sum = 0;
            for (int k = 0; k < _columns; ++k) {
                sum += this->_matrix[i][k] * obj._matrix[k][j];
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


const double Matrix::getMaxElement() const {
    double max = _matrix[0][0];
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _columns; ++j) {
            if (max < _matrix[i][j]) max = _matrix[i][j];
        }
    }
    return max;
}


void Matrix::zeroOut() {
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _columns; ++j) {
            _matrix[i][j] = 0;
        }
    }
}
void Matrix::setToUnitMatrix() {
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _columns; ++j) {
            if (i == j) _matrix[i][j] = 1;
            else _matrix[i][j] = 0;
        }
    }
}
void Matrix::clear() {
    for (int i = 0; i < _rows; ++i) {
        _matrix[i].clear();
    }
    _matrix.clear();
}

void Matrix::eliminacjaGaussaJordana(std::vector<double>& c, std::vector<double> y) const {
    Matrix A = *this;
    double l, w;
    for (int k = 0; k < _rows; ++k) {
        l = A[k][k];
        // podzielenie wiersza przez wartosc na diagonali
        for (int j = k; j < _columns; ++j) {
            A[k][j] /= l;
        }
        y[k] /= l;
        // odejmowanie wierszy
        for (int i = 0; i < _rows; ++i) {
            if (i == k)
                continue;
            w = A[i][k];
            for (int j = k; j < _columns; ++j) {
                A[i][j] -= A[k][j] * w;
            }
            y[i] = y[i] - y[k] * w;
        }
    }
    for (int i = 0; i < _rows; ++i) {
        c[i] = y[i];
    }
    //showSimpleEquation(*this, y);
}
double Matrix::rozkladLU(Matrix& L, Matrix& U) const {
    // skopiowanie macierzy
    U = *this;
    L = *this;
    L.zeroOut();
    for (int i = 0; i < _rows; ++i) { L[i][i] = 1; }

    for (int k = 0; k < _rows; ++k) {
        for (int i = k + 1; i < _rows; ++i) {
            double l = U[i][k] / U[k][k];
            L[i][k] = l;    // wyrazy macierz L, to obliczone wsółczynniki l_ik=a_ik/a_kk
            for (int j = 0; j < _columns; ++j) {
                // odejmujemy pierwszy wiersz pomnożony przez l od kolejnych wierszy
                U[i][j] = U[i][j] - U[k][j] * l;
            }
        }
    }
    double determinant = 1;
    for (int i = 0; i < _rows; ++i) {
        determinant *= U[i][i];
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
void Matrix::showInLatexForm() {
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _columns; ++j) {
            cout << _matrix[i][j] << " & ";
        }
        cout << "\\\\" << endl;
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
