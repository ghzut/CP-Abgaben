#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>

using namespace std;
using namespace Eigen;

int main()
{
  M = MatrixXd::Random();
  cout << M;

  return 0;
}
