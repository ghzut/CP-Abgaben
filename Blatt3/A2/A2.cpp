#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
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

float k_j(int n, int j)
{
  return float(n-j);
}

RowVector2f init_first_last(int n, bool first = true)
{
  RowVector2f rv;
  if(first) rv << -n+1., n-1.;
  else rv << 1./n, -1./n;
  return rv;
}



int main()
{
  int n = 3;
  MatrixXf A(n,n);
  float m;
  RowVector2f first;
  RowVector2f last;
  A = MatrixXf::Zero(n);
  first = init_first_last(n, true);
  last = init_first_last(n, false);
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
  JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);
  VectorXf ew(n);
  ew = svd.singularValues();
  for (int i = 0; i < n; ++i)
  {
    ew(i) = sqrt(ew(i));
  }
  cout << "Die Eigenfrequenzen des Systems sind: " << endl << ew;

}
