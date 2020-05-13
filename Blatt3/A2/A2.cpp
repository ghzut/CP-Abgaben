#pragma once
#include<iostream>
#include<fstream>
#include <Eigen/Dense>
#include <Eigen/SVD>

using namespace std;
using namespace Eigen;

int kroeningerdelta(int a, int b)
{
  int c;
  c = (a != b)? 0:1;
  return c;
}

double k_j(int n, int j)
{
  return double(n-j);
}

void init_first_last(int n, double arr[2], bool first = true)
{
  if(first) arr = {-n+1., n-1};
  else arr = {1./n,-1./n};
}


MatrixXf initMatrix(int n)
{
  double m;
  double first[2], last[2]
  MatrixXf A(n,n) = MatrixXf::Zero(n);
  first = init_first_last(n, first true);
  last = init_first_last(n, last, false);
  A.block(0,0,1,2) = first;
  A.block(n,n-1,1,2) = last;
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
  MatrixXf A(n,n);
  A = initMatrix(n);
  JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);
  cout << "Die Eigenfrequenzen des Systems sind: ", sqrt(svd.SingularValues());

}
