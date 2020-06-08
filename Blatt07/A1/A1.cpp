#include <iostream>
#include <Eigen/Dense>
#include <fstream>

using namespace std;
using namespace Eigen;


VectorXd get_r(double t, const VectorXd &v0)
{
  return v0;
}

VectorXd get_v(function<VectorXd(double, const VectorXd& )> f, double t, const VectorXd &r0)
{
  return -r0;
}

VectorXd rk4(function<VectorXd(double, const VectorXd&)> func, const VectorXd &y0, double h)
{
  VectorXd sum_k;
  VectorXd k_1 = func(t0, r0);
  VectorXd k_2 = func(t0 + h/2., r0 + k_1/2.);
  VectorXd k_3 = func(t0 + h/2., r0 + k_2/2.);
  VectorXd k_4 = func(t0 + h/2., r0 + k_3);
  sum_k = k_1 + 2. * k_2 + 2. * k_3 + k_4;
  sum_k *= h/6.;
  return sum_k;
}

int main()
{
  VectorXd rk;
  VectorXd vk;
  VectorXd r0;
  VectorXd v0;
  double err = 10.;
  ofstream outfile("build/A1.txt", ofstream::trunc);
  outfile << "#h, i, err\n"
  for(double h = 1.; h > 1e-8; h/=10.)
  {
    r0 << 1.,1.;
    v0 << 0.,0.;
    vk = v0 + rk4(get_v, r0, h);
    rk = r0 + rk4(get_r, vk, h);
    err = (rk - r0).abs();
    int k = 0;
    outfile << h << " " << k+1 << " " << err << "\n";
    for(k = 1; k < 10; ++k)
    {
      vk += rk4(get_v, rk, h);
      rk += rk4(get_r, vk, h);
      err = (rk-r0).abs();
      outfile << h << " " << k+1 << " " << err << "\n";
    }
  }
  outfile.flush();
  outfile.close();
  return 0;
}
