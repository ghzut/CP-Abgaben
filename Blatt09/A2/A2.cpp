#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include "math.h"
#include <random>

using namespace std;
using namespace Eigen;



MatrixXi spin(100, 100);
mt19937 generator(2); //Zufallszahlen-Generator mit konstantem Seed
uniform_real_distribution<double> distribution(0, 1);



void init(bool zufall)
{
  if (zufall == true) //zufällige Anfangsbedingungen
  {
    for (int i = 0; i < spin.rows(); i++)
    {
      for (int j = 0; j < spin.cols(); j++)
      {
        auto randomSpin = distribution(generator);

        if (randomSpin < 0.5)
        {
          spin(i,j) = -1;
        }
        else
        {
          spin(i,j) = 1;
        }
      }
    }
  }
  else //feste Anfangsbedingungen
  {
    for (int i = 0; i < spin.rows(); i++)
    {
      for (int j = 0; j < spin.cols(); j++)
      {
        spin(i,j) = -1;
      }
    }
  }
}


//Periodische Randbedingungen
int rb(int x, int y)
{

  if (x < 0)
  {
    x += spin.rows();
  }
  else if (x >= spin.rows())
  {
    x -= spin.rows();
  }

  if (y < 0)
  {
    y += spin.cols();
  }
  else if (y >= spin.cols())
  {
    y -= spin.cols();
  }

	return spin(x,y);
}


//Berechnung der Gesamtenergie über alle Spins
double energie()
{
  double E = 0;

  for (int x = 0; x < spin.rows(); x++)
  {
    for (int y = 0; y < spin.cols(); y++)
    {
      E += spin(x,y) * (rb(x, y+1) + rb(x-1, y) + rb(x, y-1) + rb(x+1, y));
    }
  }

  return -E;
}



void momentaufnahme(double beta, int sweeps, bool zufall, string filename)
{
  cout << "Start!\n";

  ofstream file;
  file.open("build/"+filename+"_anfang.txt");

  double dE;

  init(zufall);

  file << spin << "\n\n";
  file.close();

  for (int i = 0; i < sweeps; i++)
  {
    for(int j = 0; j < 10000; j++)
    {
      auto x = distribution(generator) * 100;
      auto y = distribution(generator) * 100;

      // Spin-Flip anbieten
      dE = 2 * spin(x, y) * (rb(x, y+1) + rb(x-1, y) + rb(x, y-1) + rb(x+1, y));

      if (dE <= 0 || distribution(generator) < exp(-beta * dE))
      {
        spin(x, y) *= -1;
      }
    }
  }

  file.open("build/"+filename+"_ende.txt");
  file << spin << "\n\n";

  file.close();

  cout << "Fertig!\n\n";

}

void ising(double beta, bool zufall, string filename)
{
  ofstream file;
  file.open("build/"+filename+".txt");
  file << "#schritt  energie\n\n";

  double dE, E, Eges;
  Eges = 0;

  init(zufall);

  E = energie();

  cout << "Start!\n\n";

  for (int i = 0; i < 500; i++)
  {
    for(int j = 0; j < 10000; j++)
    {

      auto x = distribution(generator) * 100;
      auto y = distribution(generator) * 100;

      // Spin-Flip anbieten
      dE = 2 * spin(x, y) * (rb(x, y+1) + rb(x-1, y) + rb(x, y-1) + rb(x+1, y));

      if (dE <= 0 || distribution(generator) < exp(-beta * dE))
      {
        spin(x, y) *= -1;
        E += dE;
      }
      // file << i+1 << "\t" << E/(spin.cols()*spin.rows()) << "\n";
    }

    Eges += E;
    file << i+1 << "\t" << Eges/((i+1)*spin.cols()*spin.rows()) << "\n";
    // cout << "Schritt  " << i+1 << "\n";

  }

  cout << "Fertig! \n\n";
  file.close();
}

int main()
{

  // a)

  momentaufnahme(1, 100000, false, "1kbt-a-fest");
  momentaufnahme(1, 100000, true, "1kbt-a-zufall");
  momentaufnahme(1/3, 100000, false, "3kbt-a-fest");
  momentaufnahme(1/3, 100000, true, "3kbt-a-zufall");

  // b)

  double Tc;
  int sweeps;

  Tc = 2 / log(1 + sqrt(2));
  sweeps = 30;

  ising(1/1.5, false, "15kbt-b-fest");
  ising(1/1.5, true, "15kbt-b-zufall");
  ising(1/Tc, false, "22kbt-b-fest");
  ising(1/Tc, true, "22kbt-b-zufall");
  ising(1/3, false, "3kbt-b-fest");
  ising(1/3, true, "3kbt-b-zufall");

  return 0;
}
