#pragma once
#include "matrix.h"
#include <algorithm>
#include <functional>


class Simplex
{
public:

Simplex(Matrix& A);

Simplex(Simplex& simplex);

double lamdaMax();

double xi();

int size;
Matrix matrix;
};
