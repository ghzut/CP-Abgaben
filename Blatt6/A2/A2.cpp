#include <iostream>
#include <Eigen/Dense>
#include <fstream>

using namespace std;
using namespace Eigen;

double f(const VectorXd &x)
{
  return x.sum();
}

double bfgs(function<double(const VectorXd&)> f, const VectorXd &x0)//, MatrixXd C0, double epsilon)
{
  double result = f(x0);
  return result;
}

int main()
{
  VectorXd x0(5);
  x0 << 6, 8, 14, 4, 10;
  double test = bfgs(f, x0);
  cout << test << "= 42?" << endl;
  return 0;
}
