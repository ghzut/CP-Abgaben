#include<iostream>
#include<fstream>
#include<cmath>
using namespace std;

void stable()
{
  //Aufgabe 2
  ofstream out_sqrts ("0_1a.txt", ofstream::trunc);
  out_sqrts << "#x, df\n";
  ofstream out_frac_trig ("0_1b.txt", ofstream::trunc);
  out_frac_trig << "#x, df\n";
  ofstream out_sin ("0_1c.txt", ofstream::trunc);

  double result = 0.;
  double stable_result = 0.;
  for (double x = 0.01; x < 1000000.; x *= 10.)
  {
    result = 1./sqrt(x) - 1./sqrt(1 + x);
    stable_result = 1./(sqrt(x)*(x+1.)+sqrt(x+1.)*x);
    out_sqrts << abs(result-stable_result)/result << "\n";
    result = (1.-cos(1./x))/sin(1./x);
    stable_result = tan(1./(2.*x));
    out_frac_trig << abs(result-stable_result)/result << "\n";
    result = sin(x+1./x) - sin(x);
    stable_result = 2*sin(1./(2.*x))*cos(x+1./(2.*x));
    out_frac_trig << abs(result-stable_result)/result << "\n";
  }
  out_sqrts.flush();
  out_sqrts.close();
  out_frac_trig.flush();
  out_frac_trig.close();
  out_sin.flush();
  out_sin.close();
}
