#include "matrix.h"
Matrix::Matrix() : row(0), column(0), matrix(vector<vector<double>>(0)) {}

Matrix::~Matrix()
{

}

Matrix::Matrix(int _row, int _column)
{
    row = _row;
    column = _column;
    matrix = vector<vector<double>>(row, vector<double>(column));
}

Matrix::Matrix(int _row, int _column, vector<vector<double> > &arr)
{
    row = _row;
    column = _column;
    matrix.assign(arr.begin(), arr.end());
}

Matrix::Matrix(const Matrix &arr)
{
    row = arr.row;
    column = arr.column;
    matrix.assign(arr.matrix.begin(), arr.matrix.end());
}

Matrix Matrix::transpose()
{
    cout << "log: Matrix::Transpose/ start\n";
    Matrix tr(this->column, this->row);
    for(int i = 0; i < row; ++i)
    {
        for(int j = 0; j < column; ++j)
        {
            tr.matrix[j][i] = this->matrix[i][j];
        }
    }

    cout << "log: Matrix::Transpose/ end\n";
    return tr;
}

Matrix Matrix::minor(const Matrix &A, int row, int col)
{
    cout << "log: Matrix::Minor/ Start\n";
    int n = A.row;
    cout << "data log: n = " << n << endl;
    Matrix m(n-1, n-1);

    int _i = 0;
    for(int i = 0; i < n; ++i)
    {
        if( i == row )
            continue;
        int _j = 0;

        for(int j = 0; j < n; ++j)
        {
            if(j == col)
                continue;
            m.matrix[_i][_j++] = A.matrix[i][j];

        }
        _i++;
    }

    cout << "data log: minor = \n" << m;
    return m;
}

double Matrix::cofactor(const Matrix &A, int row, int col)
{
    cout << "log: Matrix::cofactor/ Start\n";
    cout << "Log data: row = " << row << " col = " << col << endl;
    cout << "log data: matrix: \n" << A << endl;
    const Matrix &_minor = minor(A, row, col);
    cout << "log data: minor: \n" << _minor << endl;
    return std::pow(-1, col+row) * det(_minor);
}

double Matrix::det(const Matrix &A)
{
    cout << "log: Matrix::det(Matrix)/ start\n";
    int n = A.row;
    cout << "log: Matrix::det(Matrix)/ n = " << n << "\n";
    if(n == 1)
    {
        return A.matrix[0][0];
    }

    double d = 0;

    cout <<"Matrix::det(Matrix) :: Data log: matrix = \n" << A;

    for(int i = 0; i < n; i++)
    {
        cout << "log: Matrix::det(Matrix)/ cofactor\n";
        cout << "data log`: i = " << i << endl;
        d += A.matrix[0][i] * cofactor(A, 0, i);
    }

    return d;
}

double Matrix::det()
{
    cout << "log: Matrix::det()/ start\n";

    if(this->row != this->column)
    {
        __throw_invalid_argument;
    }
    cout << "log: Matrix::det()/Matrix::det(Matrix)\n";
    return det(*this);
}

Matrix Matrix::Inverse()
{
    cout << "log: Matrix::Inverse/ start\n";
    cout << "log: Matrix::Inverse/Matrix::det()\n";
    double d = det();
    if(d == 0)
        __throw_invalid_argument;

    Matrix A(this->row, this->column);

    for (int i = 0 ; i < this->row; ++i)
    {
        for (int j = 0; j < this->column; ++j)
        {
            A.matrix[i][j] = cofactor(*this, i, j)/ d;
        }
    }
    cout << "log: Matrix::Inverse/Matrix::transpose() A = \n";
    cout << A;
    Matrix  trans = A.transpose();
    return trans;
}

ostream& operator<<(ostream &out,const Matrix& x){
    for(int i=0; i < x.row; i++){
        for(int j = 0; j  < x.column; j++) {
            out << x.matrix[i][j] << ' ' ;
        }
        out << endl;
    }
return out;
}
