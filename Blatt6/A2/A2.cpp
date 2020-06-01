#include <iostream>
#include <Eigen/Dense>
#include <fstream>

using namespace std;
using namespace Eigen;

double f1(const VectorXd& x)
{
  if (x.size() != 2)
  {
    cerr << "Dieser Vektor ist nicht geeignet für diese Funktion." << endl;
    return -1.;
  }
  else return pow(1. - x(0), 2) + 100.*pow(x(1) - pow(x(0), 2.), 2.);
}

VectorXd g1(const VectorXd& x)
{
  VectorXd grad = VectorXd::Zero(2);
  if (x.size() != 2)
  {
    cerr << "Dieser Vektor ist nicht geeignet für diesen Gradienten." << endl;
    return grad;
  }
  else
  {
    grad(0) = 2.* (-1. + x(0) + 200 * pow(x(0),3.) - 200 * x(0) * x(1));
    grad(1) = 200 * (x(1) - pow(x(0), 2.));
    return -grad;
  }
}

MatrixXd hesse1(const VectorXd& x)
{
  MatrixXd H = MatrixXd::Zero(2);
  if (x.size() != 2)
  {
    cerr << "Dieser Vektor ist nicht geeignet für diese Hesse-Matrix." << endl;
    return H;
  }
  else
  {
    H(0,0) = 1200. * pow(x(0), 2.) - 400. * x(1) + 2.;
    H(0,1) = -400. * x(0);
    H(1,0) = H(0,1);
    H(1,1) = 200.;
    return H;
  }
}


double f1_lambda(const VectorXd &x0, const VectorXd &b0, double l0)
{
  if (x0.size() != 2)
  {
    cerr << "Dieser Vektor ist nicht geeignet für diese Funktion." << endl;
    return -1.;
  }
  return pow(1.-x0(0)-l0*b0(0), 2.) + 100*pow(-pow(b0(0)*l0,2.)-(2*x0(0)*b0(0)-b0(1))*l0 + x0(1) - pow(x0(0),2.), 2.);
}


double erste_ableitung(function<double(const VectorXd&, const VectorXd&, double)> f, double x, const VectorXd &x0, const VectorXd &b0)
{
    double h = 0.0001;

    return (f(x0, b0, x+h) - f(x0, b0, x-h))/(2*h);
}

double zweite_ableitung(function<double(const VectorXd&, const VectorXd&, double)> f, double x, const VectorXd &x0, const VectorXd &b0)
{
    double h = 0.0001;

    return (f(x0, b0, x+h) - 2*f(x0, b0, x) + f(x0, b0, x-h))/(h*h);
}

VectorXd newton(function<double(const VectorXd&, const VectorXd&, double)> f, const VectorXd &x0, const VectorXd &b0)
{
  double dx = 1e4;
  double l_0 = 1.;
  VectorXd x_new(x0.size());
  while (dx > 0.1)
  {
    dx = erste_ableitung(f, l_0, x0, b0)/zweite_ableitung(f, l_0, x0, b0);
    l_0-=dx;
  }
  x_new = x0 + l_0 * b0;
  return x_new;
}

void bfgs(function<double(const VectorXd&)> f, function<VectorXd(const VectorXd&)> g, const VectorXd &x0, const MatrixXd &C0, const double epsilon, string init)
{
  if(C0.rows() != C0.cols())
  {
    cerr << "Diese Startmatrix kann nicht symmetrish sein." << endl;
    return;
  }
  if(C0.rows() != x0.size())
  {
    cerr << "Matrix und Vektor passen nicht zusammen." << endl;
    return;
  }
  ofstream outfile("build/A2_"+init+".txt", ofstream::trunc) ;
  outfile << "# iter, err\n";
  int n = C0.rows();
  double err;
  MatrixXd I = MatrixXd::Zero(n,n);
  for(int i = 0; i < n; ++i)
  {
    I(i,i) = 1.;
  }
  // Ersten Liniensuchschritt anwenden.
  VectorXd pk;
  VectorXd yk;
  VectorXd bk1;
  double rho;
  MatrixXd Ck = C0;
  VectorXd bk = g(x0);
  VectorXd xk = newton(f1_lambda, x0, bk);
  pk = xk - x0;
  bk1 = g(xk);
  yk = bk1 - bk;
  bk = bk1;
  rho = 1./(pk.transpose()*yk);
  int iter = 0;
  err = bk.norm();

  cout << xk << endl << endl;
  cout << bk << endl << endl;
  cout << Ck << endl << endl;

  outfile << iter << " " << err << "\n";
  while (err > epsilon)
  {
    ++iter;
    pk = -Ck * bk;
    xk += pk;
    cout << xk << endl << endl;
    bk1 = g(xk);
    cout << bk1 << endl << endl;
    yk = bk1 - bk;
    bk = bk1;
    rho = 1./(pk.transpose()*yk);
    Ck = Ck - rho * pk * (yk.transpose() * Ck) - rho * (Ck * yk) * pk.transpose() + pow(rho, 2.) * pk * (yk.transpose() * (Ck * yk)) * pk.transpose() + rho * pk * pk.transpose();
    cout << Ck << endl << endl;
    err = bk.norm();
    outfile << iter << " " << err << "\n";
  }
  outfile.flush();
  outfile.close();
  cout << "Der minimierte Vektor ist\n\n" << xk << endl;
}




int main()
{
  VectorXd x0(2);
  x0 << -1.,1.;
  MatrixXd I = MatrixXd::Zero(2,2);
  I << 1., 0., 0., 1.;
  double init_3 = f1(x0);
  MatrixXd C0_1 = 4*I;//hesse1(x0);
  C0_1 = C0_1.inverse();
  MatrixXd C0_3 = I / init_3;
  bfgs(f1, g1, x0, C0_1, 1e-5, "1");
  bfgs(f1, g1, x0, C0_3, 1e-5, "3");
  return 0;
}
