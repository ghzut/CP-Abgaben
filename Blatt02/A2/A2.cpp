#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "profiler.cpp"
#include <fstream>
#include <thread>

using namespace std;
using namespace Eigen;

int main()
{
  //Multithreading zur Zeitersparnis.
  const auto nThreads = thread::hardware_concurrency();
  Eigen::setNbThreads(nThreads);

  Profiler::init(3);
  ofstream outfile("build/times.txt", ofstream::trunc);
  outfile << "#N,1,2,3\n";

  //logarithmisch steigendes N.
  for (int N = 1; N < 8193; N*=2)
  {
    //Initialisiere zufällige Problem der Größe N und stoppe die benötigte Zeit.
    Profiler::start(0);
    MatrixXd M = MatrixXd::Random(N,N);
    VectorXd b = VectorXd::Random(N);
    Profiler::stop(0);

    MatrixXd L(N,N);
    MatrixXd U(N,N);
    MatrixXd P(N,N);
    VectorXd x(N);

    //Führe die LU-Zerlegung durch und stoppe die Zeit.
    Profiler::start(1);
    L = M.partialPivLu().matrixLU().triangularView<Eigen::UpLoType::UnitLower>();
    U = M.partialPivLu().matrixLU().triangularView<Eigen::UpLoType::Upper>();
    P = M.partialPivLu().permutationP();
    Profiler::stop(1);

    //Löse das Problem Mx = b und stoppe die Zeit.
    Profiler::start(2);
    x = M.partialPivLu().solve(b);
    Profiler::stop(2);
    
    outfile << N << " " << Profiler::getTimeInS( 0 ) << " " << Profiler::getTimeInS( 1 ) << " " << Profiler::getTimeInS( 2 ) << "\n";
    Profiler::resetAll();
  }
  outfile.flush();
  outfile.close();

  ofstream outfile2("build/times_lin.txt", ofstream::trunc);
  outfile2 << "#N,1,2,3\n";

  //Hier nochmal dasselbe aus Interesse und für mehr Datenpunkte mit linear ansteigendem N.
  //Kann je nach CPU etwas lange dauern, bei mir einige Minuten.
  for (int N = 1; N < 2001; ++N)
  {
    Profiler::start(0);
    MatrixXd M = MatrixXd::Random(N,N);
    VectorXd b = VectorXd::Random(N);
    Profiler::stop(0);

    MatrixXd L(N,N);
    MatrixXd U(N,N);
    MatrixXd P(N,N);
    VectorXd x(N);

    Profiler::start(1);
    L = M.partialPivLu().matrixLU().triangularView<Eigen::UpLoType::UnitLower>();
    U = M.partialPivLu().matrixLU().triangularView<Eigen::UpLoType::Upper>();
    P = M.partialPivLu().permutationP();
    Profiler::stop(1);

    Profiler::start(2);
    x = M.partialPivLu().solve(b);
    Profiler::stop(2);

    outfile2 << N << " " << Profiler::getTimeInS( 0 ) << " " << Profiler::getTimeInS( 1 ) << " " << Profiler::getTimeInS( 2 ) << "\n";
    Profiler::resetAll();
  }
  outfile2.flush();
  outfile2.close();

  return 0;
}
