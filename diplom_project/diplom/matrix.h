#pragma once
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class Matrix 
{
public:

Matrix();
~Matrix();

Matrix(int _row, int _column);

Matrix(int _row, int _column, vector<vector<double>>& arr);


Matrix(const Matrix &arr);

Matrix transpose();

Matrix minor( const Matrix& A, int row, int col);

double cofactor(const Matrix &A, int row, int col );

double det(const Matrix &A);

friend ostream& operator<<(ostream &out, const Matrix& x);

double det();
Matrix Inverse();



vector< vector<double> > matrix;
int row;
int column;

};
