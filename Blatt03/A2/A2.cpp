//Anmerkung: Ich habe die Ruhelängen l_j nicht benutzt,
//da in den Differentialgleichungen nur die Auslenkungen
//aus der Ruhelage und nicht die absolute Position im Raum benötigt wird.

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;




int kroneckerdelta(int a, int b)
{
  int c;
  c = (a != b)? 0:1;
  return c;
}

//Funktion zur berechnung der benötigten Federkonstante.
double k_j(int n, int j)
{
  return double(n-j);
}

//Funktion zur Berechnung der 1. und letzten Reihe der Kopplungsmatrix.
//Diese müssen gesondert betrachtet werden, da die äußersten Massen nur einen Nachbarn haben.
void init_first_last(int n, MatrixXd &M)
{
  M(0,1) = -n+1.;
  M(0,0) =  n-1.;
  M(n-1,n-1) = 1./n;
  M(n-1,n-2) = -1./n;
}

//Funktion zur Initialisierung der gesamten Kopplungsmatrix mit variabler Dimension n x n.
void initMatrix(int n, MatrixXd &A)
{
  double m;
  init_first_last(n, A);
  for (int i = 1; i < n - 1; ++i)
  {
    m = i+1;
    for (int j = i-1; j < i + 2; ++j)
    {
      A(i,j) = -kroneckerdelta(i, j) * (k_j(n, j) + k_j(n, j + 1)) + kroneckerdelta(i - 1, j) * k_j(n, j+1) + kroneckerdelta(i+1,j) * k_j(n, j);
      A(i,j) /= -m;
    }
  }
}


int main()
{
  ofstream outfile("build/spektrum.txt", ofstream::trunc);
  outfile << "#n, w_i\n";
  int n = 10;

  //Um das Spektrum verschiedener Problemgrößen zu untersuchen, kann man verschiedene n untersuchen
  MatrixXd ew_Mat(n,n-2);
  for (int i = 2; i < n; ++i)
  {
    MatrixXd M(i,i);
    initMatrix(i, M);
    VectorXd ev = M.eigenvalues().real();
    for (int j = 0; j < i; ++j)
    {
      if(ev(j) > 0.00001) ew_Mat(j,i-2) = sqrt(ev(j)); //einige Egenwerte werden aufgrund von
      else ew_Mat(j,i-2) = 0.; //RUndungsfehlern als sehr kleine negative Zahlen zurückgegeben
    }
  }
  //Initialisierung der 10x10 Kopplungsmatrix und Bestimmung der Eigenwerte mithilfe von eigen.
  //Da die Matrix bereits tridiagonal ist kann sie mit n-1 Jacobi-Drehungen diagonalisiert werden.
  MatrixXd A(n,n);
  initMatrix(n, A);
  VectorXd ew = A.eigenvalues().real();
  for (int i = 0; i < n; ++i)
  {
    ew_Mat(i, n-3) = sqrt(ew(i));
  }
  outfile << ew_Mat << endl;
  outfile.close();
}
