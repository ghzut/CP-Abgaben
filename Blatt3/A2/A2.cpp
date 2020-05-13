#pragma once
#include<iostream>
#include<fstream>
#include <Eigen/Dense>
#include <Eigen/SVD>

using namespace std;
using namespace Eigen;

int kroeningerdelta(int a, int b)
{
  (a != b)? return 0 : return 1;
}

double k_j(int n, int j)
{
  return double(n-j);
}

double init_first_last[2](int n, bool first = true)
{
  if(first) return {-n+1., n-1};
  else return {1./n,-1./n};
}


MatrixXf initMatrix(int n)
{
  double m;
  MatrixXf A(n,n) = MatrixXf::Zero(n);
  A.block(0,0,1,2) = init_first_last(n, true);
  A.block(n,n-1,1,2) = init_first_last(n, false);
  for (int i = 1; i < n - 1; ++i)
  {
    m = i+1;
    for (int j = i-1; j < i + 2; ++j)
    {
      A(i,j) = -kroeningerdelta(i, j) * (k_j(n, j - 1) + k_j(n, j)) + kroeningerdelta(i - 1, j) * k_j(n, j) + kroeningerdelta(i+1,j) * k_j(n, i);
      A(i,j) /= m;
    }
  }
  return A;
}


int main()
{
  int n = 10;
  MatrixXf A(n,n) = initMatrix(n);
  JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);
  cout << "Die Eigenfrequenzen des Systems sind: ", sqrt(svd.SingularValues());

}
