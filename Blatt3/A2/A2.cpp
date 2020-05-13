#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Core>

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

RowVectorXf init_first_last(int n, bool first = true)
{
  RowVectorXf rv(n);
  if(first)
  {
    rv(0) = -n+1.;
    rv(1) =  n-1.;
  }
  else
  {
    rv(n-2) = 1./n;
    rv(n-1) = -1./n;
  }
  return rv;
}


MatrixXf initMatrix(int n)
{
  float m;
  RowVector2f first;
  RowVector2f last;
  MatrixXf A(n,n);
  first = init_first_last(n, true);
  last = init_first_last(n, false);
  A.col(0) = first;
  A.col(n-1) = last;
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
  cout << "A:" << endl << A << endl;
  JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);
  VectorXf ew(n);
  ew = svd.singularValues();
  for (int i = 0; i < n; ++i)
  {
    ew(i) = sqrt(ew(i));
  }
  cout << "Die Eigenfrequenzen des Systems sind: " << endl << ew;

}
