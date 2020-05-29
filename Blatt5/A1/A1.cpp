#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <complex>
#include "Eigen/Dense"
using namespace Eigen;
using namespace std;

typedef complex<double> cdouble;

double f1(int l)
{
  return sqrt(1.+l);
}


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

int lbar(int l, int m)
{
  VectorXd binNum = VectorXd::Zero(m);
  for (int i = 0; i < m; ++i)
  {
    binNum(i) = l % 2;
    l = l / 2;
  }
  VectorXd lqu(m);
  for (int i = 0; i < m; ++i)
  {
    lqu(i)=binNum(m-i-1);
  }
  int dec = 0;
  int base = 1;
  for (int i = 0; i < m; ++i)
  {
    dec += lqu(i)*base;
    base = base*2;
  }
  return dec;
}


MatrixXcd init_Mat(int n, const VectorXcd &v_f)
{
  MatrixXcd M = MatrixXcd::Zero(n,n);
  double n_l;
  for(int j = 0; j < n; ++j)
  {
    for(int l = 0; l < n; ++l)
    {
      n_l = double(l)/n;
      M(j, l) = omega_j_N(j, n_l) * v_f(lbar(l, log2(n)));
    }
  }
  return M;
}

VectorXcd v_F_FFT(int n, const VectorXcd &v_f)
{
  MatrixXcd M = init_Mat(n, v_f);
  VectorXcd v_FFT = VectorXcd::Zero(n);
  MatrixXcd M_temp = MatrixXcd::Zero(n,n);
  RowVectorXcd v_Fj;
  for(int j = 0; j < n; ++j)
  {
    for(int l = 0; l < n/2; ++l)
    {
      v_Fj(l) = (M(j,2*l) + M(j,2*l+1) * omega_j_N(j, 2));
    }
    M_temp.row(j) = v_Fj;
  }
  for(int i = 2; pow(2,i) < n; ++i)
  {
    for(int j = 0; j < pow(2,i); ++j)
    {
      for(int a = 0; a < n/pow(2,i); ++a)
      {
        v_Fj(a)= M_temp(j, 2 * a) + M_temp(j,2 * a + 1) * omega_j_N(j, pow(2,i)));
      }
      M_temp.row(j) = v_Fj;
      for(int b = 0; b < n/pow(2,i); ++b)
      M_temp.row(j+b*pow(2,i)) = v_Fj;
    }
  }
  for(int j = 0; j < n; ++j)
  {
    v_FFT(j) = M_temp(j,0);
  }
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
      v_f_m(l) = f1(l);
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
