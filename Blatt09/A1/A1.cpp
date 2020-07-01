#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include "math.h"
#include <random>

using namespace std;
using namespace Eigen;


double monte_carlo(int spin0, int schritte, double H)
{
  // random_device rd;
  mt19937 generator(2);


  uniform_real_distribution<double> distribution(0, 1);

  bool akzeptanz = false;
  int spin;
  double E=0., m=0.;


  for (int i = 0; i < schritte; i++)
  {
    spin = -spin0; //Vorschlag für den nächsten Schritt
    E = -(spin*H - spin0*H);

    //Überprüfen ob Vorschlag akzeptiert wird
    if (E < 0)
    {
      akzeptanz = true;
    }
    else
    {
      auto p = distribution(generator);

      if (p < exp(-E))
      {
        akzeptanz = true;
      }
      else
      {
        akzeptanz = false;
      }
    }

    //Berechnung der Magnetisierung
    if (akzeptanz == true)
    {
      if (spin < 0)
      {
        m += -1;
      }
      else
      {
        m += 1;
      }

      spin0 = spin;

    }
    else
    {
      if (spin0 < 0)
      {
        m += -1;
      }
      else
      {
        m += 1;
      }
    }

  }

  return m/schritte;
}


int main()
{
  int schritte, spin0;
  VectorXd H;
  double m;

  spin0 = 1; //Initialisierung mit beliebigem Spin
  schritte = 1e5;
  H = VectorXd::LinSpaced(1e4, -5, 5);

  ofstream file;
  file.open("build/mc.txt");
  file << "#H \t m \n\n";

  for (int i = 0; i < H.size(); i++)
  {
    m = monte_carlo(spin0, schritte, H(i));
    file << H(i) << "\t" << m << "\n";
  }

  file.flush();
  file.close();


  return 0;
}
