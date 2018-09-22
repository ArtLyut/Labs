#include "simplex.h"


Simplex::Simplex(Matrix &A)
{
    size = A.row;

    if(A.row != A.column)
        __throw_invalid_argument;

    matrix = A;

    for(int i = 0; i < size; ++i) {
        if(matrix.matrix[i][size-1] != 1)
            __throw_invalid_argument;
    }
}

Simplex::Simplex(Simplex &simplex)
{
    size = simplex.size;
    matrix = simplex.matrix;
}

double Simplex::lamdaMax()
{
    cout << "log: Simplex::lamdaMax/Matrix::Inverse/\n";

    const Matrix& A = matrix.Inverse();


    cout <<"log: Simplex::lamdaMax/ A = \n" << A << endl;
    vector<double> lambda_i(size);

    for(int i = 0; i < size; i++)
    {
        lambda_i[i] = 0;
    }

    for(int j = 0; j  < size; ++j)
    {
        for(int i = 0; i < size - 1; ++i)
        {
            if(A.matrix[i][j] < 0)
                lambda_i[j] -= A.matrix[i][j];
        }

        lambda_i[j] -= A.matrix[size-1][j];
    }

    cout << "log: Simplex::lamdaMax/std::sort/\n";
    std::sort(lambda_i.begin(), lambda_i.end(), std::greater<double>());
    cout << "log: Simplex::lamdaMax/std::sort/ end\n";
    cout << "log data: lambda_i = ";
    for(auto i : lambda_i) cout<< i << ' ';
    cout << endl;
    return lambda_i[0];
}

double Simplex::xi()
{
    return size*lamdaMax()+1;
}
