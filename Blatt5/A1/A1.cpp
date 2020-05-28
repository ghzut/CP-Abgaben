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
  o = pow(o, 2*(j+1.)*n_l);
  return o;
}

VectorXcd init_v_f(double (*func)(double), int n)
{
  VectorXcd vec = VectorXcd::Zero(n);
  return vec;
}

MatrixXcd init_Mat(int n)
{
  MatrixXcd M = MatrixXcd::Zero(n,n-1);
  double n_l;
  for(int j = 0; j < n-1; ++j)
  {
    for(int l = 0; l < n-2; ++l)
    {
      n_l = double(l)/n;
      M(j, l) = omega_j_N(j, n_l);
    }
  }
  cout << M << endl << endl;
  return M;
}

VectorXcd v_F_FFT(int n, const VectorXcd &v_f)
{
  MatrixXcd M = init_Mat(n);
  VectorXcd v_FFT = M * v_f;
  return v_FFT;
}

VectorXcd v_dir_F(const VectorXcd &v_f, int dim)
{
  VectorXcd v_dir = VectorXcd::Zero(dim);
  for(int j = 0; j < dim; ++j)
  {
    for(int l = 0; l < dim; ++l)
    {
      v_dir(j) += omega_j_N(j, double(l)/dim)*v_f(l);
    }
  }
  return v_dir;
}


int main()
{
  //1. f_l = sqrt(1 + l), l aus 2^m mit m=3,4
  ofstream outfile1("build/1_1.txt", ofstream::trunc);
  outfile1 << "# Direkt(Re,Im), FFT(Re,Im)\n";
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
    v_fft = v_F_FFT(dim, v_f_m);
    for(int i = 0; i < dim; ++i)
    {
      outfile1 << real(v_dir(i)) << " " << imag(v_dir(i)) << " " << real(v_fft(i)) << " " << imag(v_fft(i)) << "\n";
    }
    outfile1 << "\n";
  }
  outfile1.flush();
  outfile1.close();
  return 0;
}
