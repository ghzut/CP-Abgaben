#include <iostream>
#include <Eigen/Dense>
#include <fstream>

using namespace std;
using namespace Eigen;

double get_energy(const Vector3d &r)
{
  return pow(r.norm(),2)/2.;
}

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
  r0 << 42.,42.,42.;
  v0 << 0.,0.,0.;
  for(double h = 1.; h >= 1e-5; h/=10.)
  {
    MatrixXd M = MatrixXd::Zero(3,10);
    vk = v0 + rk4(get_v, r0, h, double(k));
    rk = r0 + rk4(get_r, vk, h, double(k));
    err = (rk - r0).norm();
    outfile << h << " " << 1 << " " << err << "\n";
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

  //Aufgabenteil c) Energieerhaltung, es wird das maximale h verwendet das die Toleranz aus b erfüllt
  //Für die Energie wird eigentlich eine Masse benötigt wegen E_i = m/2 * omega^2 *C_i^2. Diese wird hier m=1 gesetzt.
  double h = 1e-3;
  double energy = 0;// get_energy(r0);
  ofstream outfile3("build/A1_c.txt", ofstream::trunc);
  outfile3 << "#k, E\n";
  outfile3 << 0 << " " << energy << "\n";
  vk = v0 + rk4(get_v, r0, h, double(k));
  rk = r0 + rk4(get_r, vk, h, double(k));
  //energy = get_energy(rk);
  outfile3 << 1 << " " << energy << "\n";
  for(k = 1; k < 20; ++k)
  {
    vk += rk4(get_v, rk, h, double(k));
    rk += rk4(get_r, vk, h, double(k));
    //energy = get_energy(rk);
    outfile3 << k+1 << " " << energy << "\n";
  }
  outfile3.flush();
  outfile3.close();
  return 0;
}
