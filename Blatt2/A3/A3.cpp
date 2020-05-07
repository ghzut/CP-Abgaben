#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include "../profiler/profiler.cpp"
#include <fstream>
#include <thread>

using namespace std;


Eigen::VectorXd solve_normally(Eigen::MatrixXd matrix, Eigen::VectorXd vector, ofstream& file){
    Eigen::VectorXd x_vec;
    Profiler::start(0);
    x_vec = matrix.inverse() * vector;
    Profiler::stop(0);

    file << Profiler::getTimeInS( 0 ) << "  ";

    return x_vec;
}

Eigen::VectorXd partial_LU(Eigen::MatrixXd matrix, Eigen::VectorXd vector, ofstream& file){
    Eigen::VectorXd x_vec;
    Profiler::start(1);
    x_vec = matrix.partialPivLu().solve(vector);
    Profiler::stop(1);
    
    file << Profiler::getTimeInS( 1 ) << "  ";

    return x_vec;
}

Eigen::VectorXd full_LU(Eigen::MatrixXd matrix, Eigen::VectorXd vector, ofstream& file){
    Eigen::VectorXd x_vec;
    Profiler::start(2);
    x_vec = matrix.fullPivLu().solve(vector);
    Profiler::stop(2);
    file << Profiler::getTimeInS( 2 ) << "  ";

    return x_vec;
}

double calc_dev(Eigen::VectorXd& exact, Eigen::VectorXd& estimate){
    double abs_ex = exact.squaredNorm();
    double abs_est = estimate.squaredNorm();
    double deviation = abs(1 - abs_est/abs_ex)*100;
    return deviation;
}

int main(){
    
    srand(42);

    const auto nThreads = std::thread::hardware_concurrency();
    Eigen::setNbThreads(nThreads);
    Profiler::init(3);

    ofstream times("build/times.txt", ofstream::trunc);
    times << "#     Normal     Partial_LU     Full_LU \n\n";

    ofstream devs("build/deviations.txt", ofstream::trunc);
    devs << "#  Partial_LU     Full_LU \n\n";


    for(int N = 1; N <= 1000; N+=2){

    Eigen::MatrixXd M = Eigen::MatrixXd::Random(N,N);
    
    while(M.determinant() == 0){
            M = Eigen::MatrixXd::Random(M.rows(),M.cols());
          }

    Eigen::VectorXd b = Eigen::VectorXd::Random(N);

    times << N << "     ";
    Eigen::VectorXd x_vec_normal = solve_normally(M, b, times);
    Eigen::VectorXd x_vec_pLU = partial_LU(M, b, times);
    Eigen::VectorXd x_vec_fLU = full_LU(M, b, times);
    times << "\n";

    double partial_LU_err = calc_dev(x_vec_normal, x_vec_pLU);
    double full_LU_err = calc_dev(x_vec_normal, x_vec_fLU);

    devs << N << "  " << partial_LU_err << "    " << full_LU_err << "\n";

    Profiler::resetAll ();
    }
    times.flush();
    times.close();
    devs.flush();
    devs.close();

    return 0;
}