#include <iostream>
#include <math.h>
#include <fstream>
#include <complex>
#include "Eigen/Dense"
using namespace Eigen;
using namespace std;

typedef complex<double> cdouble;

cdouble omega_j_N(int j, double n_l)
{
  cdouble o = -1.;
  o = pow(o, 2*(j+1.)/n_l);
  return o;
}

VectorXcd init_v_f(double (*func)(double), int n)
{
  VectorXcd vec = VectorXcd::Zero(n);
  return vec;
}

void init_Mat(MatrixXcd &M, int n, const VectorXcd &v_f)
{
  M = MatrixXd::Zero(n,n);
  double n_l;
  for(int j = 0; j < n-1; ++j)
  {
    for(int l = 0; l < n-2; ++l)
    {
      n_l = double(l)/n;
      M(j, l) = omega_j_N(j, n_l);
    }
  }
}

VectorXcd v_dir_F(const VectorXcd &v_f, int dim)
{
  VectorXcd v_dir = VectorXcd::Zero(dim);
  for(int j = 0; j < dim; ++j)
  {
    for(int l = 0; l < dim; ++l)
    {
      cout << omega_j_N(j, double(l)/dim) << endl;
      v_dir(j) += omega_j_N(j, double(l)/dim)*v_f(l);
    }
    cout << endl;
  }
  return v_dir;
}

int main()
{
  //1. f_l = sqrt(1 + l), l aus 2^m mit m=3,4
  ofstream outfile1("A1/build/1_1.txt", ofstream::trunc);
  outfile1 << "# Direkt(Re,Im), FFT(Re,Im)";
  for(int m = 3; m < 5; ++m)
  {
    int dim = pow(2,m);
    VectorXcd v_f_m = VectorXcd::Zero(dim);
    for(int l = 0; l < dim; ++l)
    {
      v_f_m(l) = sqrt(1.+l);
    }
    VectorXcd v_dir, v_fft;
    v_dir = v_dir_F(v_f_m, dim);
    cout << endl;
    //cout << v_dir << endl << endl;
  }
  return 0;
}
