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

// Lange Geschichte, ich sollte das so nennen.
// Nicht irritieren lassen, gemeint ist natürlich das Kronecker-Delta,
// um zu überprüfen in welchem Matrixeintrag wir sind.
int kroeningerdelta(int a, int b)
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
//Seltsamerweise funktionierte es nicht einen 1x2 Vektor zu initialisieren und per Matrix.block() zuzuweisen?
RowVectorXd init_first_last(int n, bool first = true)
{
  RowVectorXd rv(n);
  if(first)
  {
    rv(1) = -n+1.;
    rv(0) =  n-1.;
  }
  else
  {
    rv(n-1) = 1./n;
    rv(n-2) = -1./n;
  }
  return rv;
}

//Funktion zur Initialisierung der gesamten Kopplungsmatrix mit variabler Dimension n x n.
MatrixXd initMatrix(int n)
{
  double m;
  RowVectorXd first;
  RowVectorXd last;
  MatrixXd A = MatrixXd::Zero(n,n);
  first = init_first_last(n, true);
  last = init_first_last(n, false);
  A.row(0) = first;
  A.row(n-1) = last;
  for (int i = 1; i < n - 1; ++i)
  {
    m = i+1;
    for (int j = i-1; j < i + 2; ++j)
    {
      A(i,j) = -kroeningerdelta(i, j) * (k_j(n, j) + k_j(n, j + 1)) + kroeningerdelta(i - 1, j) * k_j(n, j+1) + kroeningerdelta(i+1,j) * k_j(n, j);
      A(i,j) /= -m;
    }
  }
  return A;
}


int main()
{
  ofstream outfile("build/spektrum.txt", ofstream::trunc);
  outfile << "#n, w_i\n";
  int n = 10;
  for (int i = 2; i < n; ++i)
  {
    MatrixXd M = initMatrix(i);
    VectorXd ev = M.eigenvalues().real();
    cout << M << "\n\n" << ev << "\n\n";
  }
  //Initialisierung der 10x10 Kopplungsmatrix und Bestimmung der Eigenwerte mithilfe von eigen.
  //Da die Matrix bereits tridiagonal ist kann sie mit n-1 Jacobi-Drehungen diagonalisiert werden.
  MatrixXd A = initMatrix(n);
  VectorXd ew = A.eigenvalues().real();
  for (int i = 0; i < n; ++i)
  {
    ew(i) = sqrt(ew(i));
  }
//  cout << "Die Eigenfrequenzen des 10x10 Systems sind: " << endl << ew;

}
