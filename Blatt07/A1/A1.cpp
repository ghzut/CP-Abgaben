#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include "math.h"

using namespace std;
using namespace Eigen;

//Bestimme die Gesamtenergie des Systems
double get_energy(const Vector3d &r, const Vector3d &v)
{
  return (pow(r.norm(),2.) + pow(v.norm(),2.))/2.;
}


//Für r beliebig und v0=0 fallen in r alle sin-Anteile weg die restlichen Konstanten sind dann bestimmt durch r0
Vector3d get_r_1(double t, const Vector3d &v)
{
  Vector3d new_r = -v*sin(t);
  return new_r;
}

Vector3d get_v_1(double t, const Vector3d &r)
{
  Vector3d new_v;
  new_v(0) = -r(0)*cos(t);
  return new_v;
}

//Für r beliebig, nicht parallel zu v0 und v0=/=0 fallen in r die sin-Anteile nicht weg, die restlichen Konstanten sind wieder bestimmt durch r0
Vector3d get_r_2(double t, const Vector3d &v)
{
  Vector3d new_r = -v*sin(t);
  new_r(1) = cos(t);
  return new_r;
}

Vector3d get_v_2(double t, const Vector3d &r)
{
  Vector3d new_v = -r*cos(t);
  new_v(1) += sin(t);
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
  for(int j = 0; j < 2; ++j)
  {
    double err = 10.;
    ofstream outfile("build/A1_"+to_string(j)+".txt", ofstream::trunc);
    ofstream outfile2("build/A1_Mat_"+to_string(j)+".txt", ofstream::trunc);
    outfile2 << "#Matrix der Ergebnisvektoren\n";
    outfile << "#h, i, err\n";
    r0 << 42.,42.,42.;
    if(j==0) v0 << 0.,0.,0.;
    else v0 << 0.,1.,0.;
    for(double h = 1.; h >= 1e-7; h/=10.)
    {
      MatrixXd M = MatrixXd::Zero(3,10);
      k = 0;
      if(j==0)
      {
        vk = v0 + rk4(get_v_1, r0, h, k*2*M_PI);
        rk = r0 + rk4(get_r_1, vk, h, k*2*M_PI);
      }
      else
      {
        vk = v0 + rk4(get_v_2, r0, h, k*2*M_PI);
        rk = r0 + rk4(get_r_2, vk, h, k*2*M_PI);
      }
      err = (rk - r0).norm();
      outfile << h << " " << k+1 << " " << err << "\n";
      M.col(k) = rk;
      for(k = 1; k < 10; ++k)
      {
        if(j==0)
        {
          vk += rk4(get_v_1, rk, h, k*2*M_PI);
          rk += rk4(get_r_1, vk, h, k*2*M_PI);
        }
        else
        {
          vk += rk4(get_v_2, rk, h, k*2*M_PI);
          rk += rk4(get_r_2, vk, h, k*2*M_PI);
        }
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
  }

  //Um zu zeigen, dass tatsächlich ein harmonischer Oszillator wird außerdem der Vektor bei t = (2n+1)*pi bestimmz
  v0 << 0.,0.,0.;
  double h = 1e-3;
  ofstream outfi("build/A1_Mat_komp.txt");
  outfi << "#Mat\n";
  MatrixXd M_komp = MatrixXd::Zero(3,10);
  k = 0;
  vk = v0 + rk4(get_v_1, r0, h, (k*2+1)*M_PI);
  rk = r0 + rk4(get_r_1, vk, h, (k*2+1)*M_PI);
  M_komp.col(k) = rk;
  for(k = 1; k < 10; ++k)
  {
    vk += rk4(get_v_1, rk, h, (k*2+1)*M_PI);
    rk += rk4(get_r_1, vk, h, (k*2+1)*M_PI);
    M_komp.col(k) = rk;
  }
  outfi <<  M_komp << endl;
  outfi.close();
  //Aufgabenteil c) Energieerhaltung, es wird das maximale h und das nächst kleinere verwendet, das die Toleranz aus b erfüllt
  //Für die Energie wird eigentlich eine Masse benötigt wegen E_i = m/2 * (v_i + omega^2 * x_i^2). Diese wird hier m=1 gesetzt.
  for(int i = 0; i <2; ++i, h/=10.)
  {
    double energy = get_energy(r0, v0);
    ofstream outfile3("build/A1_c_"+to_string(i)+".txt", ofstream::trunc);
    outfile3 << "#k, E\n";
    outfile3 << 0 << " " << energy << "\n";
    k = 0;
    vk = v0 + rk4(get_v_1, r0, h, k*2*M_PI);
    rk = r0 + rk4(get_r_1, vk, h, k*2*M_PI);
    energy = get_energy(rk, vk);
    outfile3 << 1 << " " << energy << "\n";
    for(k = 1; k < 20; ++k)
    {
      vk += rk4(get_v_1, rk, h, k*2*M_PI);
      rk += rk4(get_r_1, vk, h, k*2*M_PI);
      energy = get_energy(rk, vk);
      outfile3 << k+1 << " " << energy << "\n";
    }
    outfile3.flush();
    outfile3.close();
  }

  return 0;
}
