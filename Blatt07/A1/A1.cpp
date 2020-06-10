#include <iostream>
#include <Eigen/Dense>
#include <fstream>

using namespace std;
using namespace Eigen;


Vector3d get_r(double t, const Vector3d &v)
{
  Vector3d new_v = -v*sin(t);
  return new_v;
}

Vector3d get_v(double t, const Vector3d &r)
{
  Vector3d new_v = -r*cos(t);
  return new_v;
}

Vector3d rk4(function<Vector3d(double, const Vector3d&)> func, const Vector3d &y, double h, double t0)
{
  Vector3d sum_k;
  Vector3d k_1 = func(t0, y);
  Vector3d k_2 = func(t0 + h/2., y + k_1/2.);
  Vector3d k_3 = func(t0 + h/2., y + k_2/2.);
  Vector3d k_4 = func(t0 + h, y + k_3);
  sum_k = k_1 + 2. * k_2 + 2. * k_3 + k_4;
  sum_k *= h/6.;
  return sum_k;
}

int main()
{
  int k;
  Vector3d rk;
  Vector3d vk;
  Vector3d r0;
  Vector3d v0;
  double err = 10.;
  ofstream outfile("build/A1.txt", ofstream::trunc);
  ofstream outfile2("build/A1_Mat.txt", ofstream::trunc);
  outfile2 << "#Matrix der Ergebnisvektoren\n";
  outfile << "#h, i, err\n";
  for(double h = 1.; h >= 1e-7; h/=10.)
  {
    MatrixXd M = MatrixXd::Zero(3,10);
    r0 << 1.,0.,0.;
    v0 << 0.,0.,0.;
    k = 0;
    vk = v0 + rk4(get_v, r0, h, double(k));
    rk = r0 + rk4(get_r, vk, h, double(k));
    err = (rk - r0).norm();
    outfile << h << " " << k+1 << " " << err << "\n";
    M.col(k) = rk;
    for(k = 1; k < 10; ++k)
    {
      vk += rk4(get_v, rk, h, double(k));
      rk += rk4(get_r, vk, h, double(k));
      err = (rk-r0).norm();
      outfile << h << " " << k+1 << " " << err << "\n";
      M.col(k) = rk;
    }
    outfile2 << M << "\n\n";
  }
  outfile.flush();
  outfile.close();
  outfile2.flush();
  outfile2.close();
  return 0;
}
