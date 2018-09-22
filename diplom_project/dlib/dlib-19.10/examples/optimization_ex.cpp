
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>
#include <iostream>

#include <vector>
#include <iomanip>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <random>

#include <fstream>
#include <algorithm>
#include <string>
#include <iomanip>



using namespace dlib;
using namespace std;
using namespace Eigen;
using namespace std;
using namespace boost;

#include <fstream>
#include <string>

void LogMSG(const std::string &s)
{
    return;
    fstream out;
    out.open("/Users/artem/diplomlogs/diplomlogs.txt", ios::out);
    cout << s << endl;
    out.close();
}


int SYMPLEX_DIM;
typedef matrix<double,0,1> column_vector;

//Генератор случайных чисел
double sample(double dummy)
{
    using namespace std::chrono;
    std::default_random_engine engine(
        system_clock::to_time_t(system_clock::now()));
    std::uniform_real_distribution<> distr(0, 1);
    return distr(engine);
}

//Получение случайной матрицы
void getRandomMatrix( MatrixXd& a, int n)
{
    using namespace std::chrono;
    std::default_random_engine engine(
        system_clock::to_time_t( system_clock::now()) );
    std::uniform_real_distribution<> distr(0, 1);
    auto gen_number = [&engine, &distr] () { return distr(engine); };
    for( int i = 0; i <= n; ++i )
    {
        for( int j = 0; j < n; ++j )
        {
            a(i, j) = gen_number();
        }
    }

    for( int i = 0; i <= n; ++i )
    {
        a(i, n) = 1;
    }

}


//Вычисление нормы проектора
double getNorm( MatrixXd& a)
{
    LogMSG("getNorm/ start");
    int dim = SYMPLEX_DIM + 1;
    LogMSG("getNorm/ MatrixXd::inverse");
    MatrixXd invMatrix = a.inverse();
    //cout << a << endl;
    //cout << "detA = " << a.determinant() << endl;

    int pow2 = 1 << (dim-1);
    double norm = 0;
    LogMSG( "getNorm/ start search norm" );
    LogMSG( "getNorn/ number of cube vertex = " + std::to_string( pow2 ) );
    for(int i = 0; i < pow2; ++i)
    {
        std::vector<double> x(dim, 0);
        x[dim-1]=1;
        int ind = 0;
        int d = i;
        while( d )
        {
            x[ind] = d % 2;
            ind++;
            d /= 2;
        }
        double sum = 0;
        for( int k = 0; k < dim; ++k )
        {
            double cur_sum = 0;
            for( int j = 0; j < dim; ++j )
            {
                cur_sum += invMatrix(j,k)*x[j];
            }
            sum += std::abs( cur_sum );
        }

        norm = std::max( sum, norm );
    }
    LogMSG( "getNorm/ return norm" );
    return norm;
}

//Составить матрицу по вектору начального приближения
MatrixXd getMatixByVect( const column_vector& m )
{
    LogMSG("getMatixByVect/ start");
    MatrixXd A = MatrixXd::Zero( SYMPLEX_DIM+1,SYMPLEX_DIM+1 );
    for(int i = 0; i <= SYMPLEX_DIM; ++i)
    {
        for( int j = 0; j < SYMPLEX_DIM; ++j )
        {
            A(i, j) = m( SYMPLEX_DIM*i + j );
        }
        A(i, SYMPLEX_DIM) = 1;
    }
    return A;
}

//Составить матрицу по вектору
column_vector getVectorByMatrix( MatrixXd& A )
{
    column_vector vect( SYMPLEX_DIM*( SYMPLEX_DIM+1 ) );
    for( int i = 0; i <= SYMPLEX_DIM; ++i )
    {
        for( int j = 0; j < SYMPLEX_DIM; ++j )
        {
            vect(SYMPLEX_DIM*i + j) = A(i, j);

        }
    }
    return vect;
}

//Целевая функция для оптимизации
double func( const column_vector& m )
{
    LogMSG( "func/ start" );
    LogMSG( "func/ getMatrixByVect" );
    MatrixXd A =getMatixByVect( m );
    LogMSG( "func/ getNorm" );
    return getNorm( A );
}

//Метод для посчета минимальной нормы на
//случайно сгенерированном наборе симплексов
double MonteCarlo( int dimension )
{
    double minimalNorm = 1e9+7;
    for( int i = 0; i < 1e10; ++i )
    {
        MatrixXd a = MatrixXd::Zero( dimension+1,dimension+1 );
        getRandomMatrix( a, dimension );
        minimalNorm = std::min(minimalNorm, getNorm( a ));
    }
    return minimalNorm;
}

void solve( int n, column_vector& starting_point )
{
    try
    {
        //n-мерный случай
        SYMPLEX_DIM = n;
        cout << starting_point << endl;
        LogMSG( "main/ find_min_box_constrained" );
        cout << endl <<  find_min_box_constrained(
                    bfgs_search_strategy(),
                    objective_delta_stop_strategy( 1e-15 ).be_verbose(),
                    func,
                    derivative( func,1e-15 ),
                    starting_point, 0, 1 );
        cout << "--------"<< endl;

        cout << "main/getMatixByVect" << endl;
        MatrixXd A1 = getMatixByVect( starting_point );
        cout << A1 << endl;
        cout << "main/getNorm" << endl;
        cout << "theta(" << n << ") = " << setprecision(11)
             << getNorm(A1) << endl << "---------------\n";
    }
    catch(...){
        cout << "Процесс вычесления прерван" << endl;
    }
}

void test_solve()
{
    return;
}

void ReadDataFromConsole( int& n, column_vector& r_data )
{
    std::cout << "ReadDataFromConsole" << endl;
    std::cin >> n;
    r_data = column_vector( n*(n+1) );
    for( int i = 0; i <= n; ++i )
    {
        for( int j = 0; j < n; ++j )
        {
            double coord;
            std::cin >> coord;
            r_data(n*i + j) = coord;
        }
    }

    return;
}

void ReadDataFromFile( const std::string& file_name,
                       int &n,
                       column_vector& r_data )
{
    std::cout << "ReadDataFromFile" << endl;
    std::fstream input_file( file_name, ios_base::in );
    input_file >> n;
    r_data = column_vector( n*(n+1) );
    for( int i = 0; i <= n; ++i )
    {
        for( int j = 0; j < n; ++j )
        {
            double coord;
            input_file >> coord;
            r_data( n*i + j ) = coord;
        }
    }
    return;
}

void solve_with_random_start_point( int n )
{
    MatrixXd B = MatrixXd::Zero( n+1, n+1 );
    getRandomMatrix( B, n );
    column_vector starting_point = getVectorByMatrix( B );
    solve( n, starting_point );
}

int main() try
{

    int TYPE_INPUT = 0;
    // type input. Default input from console type = 0, from file type = 1;
    std::string NAME_INPUT_FILE;
    std::cin >> TYPE_INPUT;
    int n;
    column_vector starting_point;

    switch( TYPE_INPUT )
    {
        case 0:
            ReadDataFromConsole( n, starting_point );
            break;
        case 1:
            std::cin >> NAME_INPUT_FILE;
            ReadDataFromFile( NAME_INPUT_FILE, n, starting_point );
            break;
        default:
            std::cout << "Error input, check type input!";
            return 0;

    }

    //Примеры решения и стартовые приближения для n = 1..6

    solve( n, starting_point );

    column_vector starting_point2 = {0.45, 0.45,
                                     0.8, 0.69,
                                     0.3, 0.9};
    solve(2, starting_point2);
    column_vector starting_point3 = {0.547079, 0.0126172, 0.0867465,
                                     0.490059, 0.478761,  0.635784,
                                     0.93395, 0.514259, 0.193468,
                                     0.217008, 0.73415, 0.467616};
    solve(3, starting_point3);
    column_vector starting_point4 = {1, 0.292919, 0, 0,
                                     1, 1, 0.707131, 1,
                                     0.499911, 0, 1, 0.500076,
                                     0, 0.292961, 0, 1,
                                     0, 1, 0.70711, 0};
    solve(4, starting_point4);

    column_vector starting_point5 = {1, 1, 1, 0.282777, 0.662421,
                                     0, 0.282777, 0.662421, 0, 1,
                                     0.337578, 0, 1, 1, 0.282777,
                                     0, 1, 0.282777, 0.662421, 0,
                                     0.717222, 0.662421, 0, 1, 1,
                                     1, 0, 0, 0, 0};
    solve(5, starting_point5);
    column_vector starting_point6 = {1, 0, 1, 0.091, 0, 0.4999,
                                     0, 1, 1, 1, 0, 0,
                                     0.5, 1, 1, 0.091, 1, 1,
                                     0.9106, 0.0896, 0.0895, 1, 0.91053, 0.9106,
                                     0, 0, 0.5, 0.09115, 1, 0,
                                     0, 0.49999, 0, 0.09104, 0, 1,
                                     1, 1, 0, 0.09105, 0.5001, 0};
    solve(6, starting_point6);


}
catch (std::exception& e)
{
    cout << e.what() << endl;
}
catch(...)
{
    cout << "unknown exception"<<endl;
}
