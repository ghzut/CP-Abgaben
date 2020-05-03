#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>
#include <stdexcept>
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


int main(){
 
    const auto nThreads = std::thread::hardware_concurrency();
    Eigen::setNbThreads(nThreads);
    Profiler::init(3);
    ofstream times("build/times.txt", ofstream::trunc);
    times << "#     Normal     Partial_LU     Full_LU \n\n";

    for(int N = 1; N <= pow(2, 13); N*=2){

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
    }
    times.flush();
    times.close();
    //results << x_vec_normal <<  x_vec_pLU  << x_vec_fLU << "\n" ;

    return 0;
}