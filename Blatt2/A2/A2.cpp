#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "profiler.cpp"
#include <fstream>

using namespace std;
using namespace Eigen;

int main()
{
  Profiler::init(3);
  ofstream outfile("build/times.txt", ofstream::trunc);
  outfile << "#N,1,2,3\n";
  for (int N = 1; N < 8193; N*=2)
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
    outfile << N << " " << Profiler::getTimeInS( 0 ) << " " << Profiler::getTimeInS( 1 ) << " " << Profiler::getTimeInS( 2 ) << "\n";
  }
  outfile.flush();
  outfile.close();

  return 0;
}
