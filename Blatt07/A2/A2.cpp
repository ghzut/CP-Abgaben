#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <math.h>

using namespace std;
using namespace Eigen;

double pot(const VectorXd& r){
  double Pot = 0.5*r.dot(r);
  return Pot;
}

VectorXd dgl(const VectorXd& x, const VectorXd& x_punkt,const double& alpha){
  VectorXd DGL = -x-alpha*x_punkt;
  return DGL;
}

VectorXd next_step(const VectorXd& y, const double& alpha){
  unsigned int d = y.size()/2;
  VectorXd temp(2*d);

  temp.segment(0,d) = y.segment(d,d);
  temp.segment(d,d) = dgl(y.segment(0,d), y.segment(d,d),alpha);
  return temp; 
}

MatrixXd runge_kutta(VectorXd (*f)(const VectorXd& , const double&), double T, int N, double alpha, const VectorXd& r, const VectorXd& v){
  double h = T/N;
  unsigned int d = r.size();
  VectorXd tn(N+1), k1, k2, k3, k4, y(2*d), y_next(2*d);
  MatrixXd ergebnis(2*d, N+1);

  // Die Startwerte in die Matrix schreiben
  ergebnis.col(0).segment(0,d) = r;
  ergebnis.col(0).segment(d,d) = v;

  for (int i = 0; i<=N; i++){
    tn(i) = i*h;
  }

  y.segment(0,d) = r;
  y.segment(d,d) = v;
  for (int i = 1; i < N+1; i++){
    k1 = h*f(y, alpha);
    k2 = h*f(y+0.5*k1, alpha);
    k3 = h*f(y+0.5*k2, alpha);
    k4 = h*f(y+k3, alpha);
    y_next = y + 1./6.*(k1 + 2*k2 + 2*k3 + k4);
    y = y_next;
    ergebnis.col(i) = y;
  }
  return ergebnis;
}


void adams_bashfort(VectorXd (*f)(const VectorXd&, const double&), double T, int N, double alpha, VectorXd& x, VectorXd& x_punkt, ofstream &file, VectorXd &energie){
  double h = T/N;
  unsigned int d = x.size();
  VectorXd y(2*d), tn(N+1);
  MatrixXd ergebnis(2*d, N+1), runge(2*d, 4);

  y.segment(0,d) = x;
  y.segment(d,d) = x_punkt;
  runge = runge_kutta(next_step, T*3.0/N, 3, alpha, x, x_punkt);

  ergebnis.col(0) = runge.col(3) + h/24.0*(55*f(runge.col(3), alpha) - 59*f(runge.col(2), alpha) +37*f(runge.col(1), alpha) -9*f(runge.col(0), alpha));

  ergebnis.col(1) = ergebnis.col(0) + h/24.0*(55*f(ergebnis.col(0), alpha) - 59*f(runge.col(3), alpha) +37*f(runge.col(2), alpha) -9*f(runge.col(1), alpha));

  ergebnis.col(2) = ergebnis.col(1) + h/24.0*(55*f(ergebnis.col(1), alpha) - 59*f(ergebnis.col(0), alpha) +37*f(runge.col(3), alpha) -9*f(runge.col(2), alpha));

  ergebnis.col(3) = ergebnis.col(2) + h/24.0*(55*f(ergebnis.col(2), alpha) - 59*f(ergebnis.col(1), alpha) +37*f(ergebnis.col(0), alpha) -9*f(runge.col(3), alpha));

  // Energie für die ersten vier Schritte bestimmen
  energie(0) = 0.5*ergebnis.col(0).segment(d,d).squaredNorm() + pot(ergebnis.col(0).segment(0,d));
  energie(1) = 0.5*ergebnis.col(1).segment(d,d).squaredNorm() + pot(ergebnis.col(1).segment(0,d));
  energie(2) = 0.5*ergebnis.col(2).segment(d,d).squaredNorm() + pot(ergebnis.col(2).segment(0,d));
  energie(3) = 0.5*ergebnis.col(3).segment(d,d).squaredNorm() + pot(ergebnis.col(3).segment(0,d));

  for(int i = 4; i < N+1; i++){
    y = ergebnis.col(i-1) + h/24.0*(55*f(ergebnis.col(i-1), alpha) - 59*f(ergebnis.col(i-2), alpha) + 37*f(ergebnis.col(i-3), alpha) - 9*f(ergebnis.col(i-4), alpha));
    ergebnis.col(i) = y;
    energie(i) = 0.5*y.segment(d,d).squaredNorm()+ pot(y.segment(0,d));
  }

  for (int i = 0; i<=N; i++){
    tn(i) = i*h;
    file << tn(i) << " ";
  }
  file << endl;

  for(int i = 0; i<ergebnis.rows()/2; i++){
    for(int j = 0; j<ergebnis.cols(); j++){
      file << ergebnis(i, j) << " ";
    }
    file << endl;
  }
}


int main() {
  unsigned int d = 3;
  double alpha = 0., T = 20.;
  int N = 300;
  VectorXd x(d), x_punkt(d), y(2*d), energie(N+1);
  ofstream file;

  //Startvektoren
  x << 1, 1, 1;
  x_punkt << 1, 2, 1;
  y.segment(0,d) = x;
  y.segment(d,d) = x_punkt;

  // Aufgabenteil a)

  //Harmonischer Oszilator
  file.open("build/aufg1_a1.txt", ios::trunc);
  adams_bashfort(next_step, T, N, alpha, x, x_punkt, file, energie);
  file.close();

  //Aperiodischer Grenzfall
  alpha = 2;
  file.open("build/aufg1_a2.txt", ios::trunc);
  adams_bashfort(next_step, T, N, alpha, x, x_punkt, file, energie);
  file.close();

  //Kriechfall
  alpha = 4;
  file.open("build/aufg1_a3.txt", ios::trunc);
  adams_bashfort(next_step, T, N, alpha, x, x_punkt, file, energie);
  file.close();

  // Gedämpfte Schwingung
  alpha = 0.1;
  file.open("build/aufg1_a4.txt", ios::trunc);
  adams_bashfort(next_step, T, N, alpha, x, x_punkt, file, energie);
  file.close();

  // Aufgabenteil b)
  file.open("build/aufg1_b.txt", ios::trunc);
  file << "Zeit Energie" << endl;
  for(int i=0; i<=N; i++)
  {
      file << setprecision(10) << i*T/N << " " << energie(i) << endl;
  }
  file.close();
  return 0;
}
