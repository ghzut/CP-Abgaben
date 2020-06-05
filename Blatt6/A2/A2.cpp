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
  MatrixXd H = MatrixXd::Zero(x.size(),x.size());
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

double f2(const VectorXd& x)
{
  if (x.size() != 2)
  {
    cerr << "Dieser Vektor ist nicht geeignet für diese Funktion." << endl;
    return -1.;
  }
  else return 1./(1. + exp(-10. * pow(x(0) * x(1) - 3., 2.))/(pow(x(0),2.) + pow(x(1),2.)));
}

VectorXd g2(const VectorXd x)
{
  if (x.size() != 2)
  {
    cerr << "Dieser Vektor ist nicht geeignet für diesen Gradienten." << endl;
    return -1.;
  }
  VectorXd grad(2);
  grad(0) = (2.*exp(10.*x(0)*x(1) + 30.)*(5.*pow(x(0),2.)*x(1) + x(0) + 5.*pow(x(1),3))/pow(exp(10.*x(0)*x(1))*(pow(x(0),2.) + pow(x(1),2.)) + exp(30.)),2.);
  grad(1) = (2.*exp(10.*x(0)*x(1) + 30.)*(5.*(pow(x(0),2.)+pow(x(1),2.))*x(0) + x(1) )/pow(exp(10.*x(0)*x(1))*(pow(x(0),2.) + pow(x(1),2.)) + exp(30.)),2.);
  return grad;
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

VectorXd g2_to_h2(const VectorXd& x, bool x_or_y)
{
  VectorXd h_vec;
  double h = 0.0001;
  if(x_or_y)
  {
    h_vec << h,0.;
  }
  else h_vec << 0.,h;
  return (g2(x + h_vec)-g2(x + h_vec))/(2*h);
}

MatrixXd hesse2(const VectorXd& x, bool diag)
{
  MatrixXd hesse2(x.size(),x.size());
  VectorXd hes_col = g2_to_h2(x, true);
  hesse2.col(0) = hes_col;
  hes_col = g2_to_h2(x, false);
  hesse2.col(1) = hes_col;
  if(diag)
  {
    hesse2(0,1) = 0.;
    hesse2(1,0) = 0.;
  }
  return hesse2;
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

void bfgs(function<double(const VectorXd&)> f, function<VectorXd(const VectorXd&)> g, const VectorXd &x0, const MatrixXd &C0, const double epsilon, string init, bool linie=false)
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
  outfile << "# iter, bk, r\n";
  int n = C0.rows();
  double err;
  MatrixXd I = MatrixXd::Zero(n,n);
  for(int i = 0; i < n; ++i)
  {
    I(i,i) = 1.;
  }
  Vector2d min; min << 1,1;
  // Ersten Liniensuchschritt anwenden.
  VectorXd pk;
  VectorXd yk;
  VectorXd bk1;
  double rho;
  MatrixXd Ck = C0;
  VectorXd bk = g(x0);
  VectorXd xk;
  if(linie)//Zur Überprüfung, ob der Algorithmus mit zusätzlichen Liniensuchschritten besser konvergiert
  {
  xk = newton(f1_lambda, x0, -C0*bk);
  pk = xk - x0;
  bk1 = g(xk);
  yk = bk1 - bk;
  rho = 1./(pk.transpose()*yk);
  bk = bk1;
  }
  else xk = x0;
  int iter = 0;
  err = bk.norm();
  double r = (min-xk).norm();
  outfile << iter << " " << err << " " << r << "\n";
  while (err > epsilon)
  {
    ++iter;
    if(linie)//Zur Überprüfung, ob der Algorithmus mit zusätzlichen Liniensuchschritten besser konvergiert
    {
      VectorXd temp;
      temp = newton(f1_lambda, xk, - Ck * bk);
      pk = temp - xk;
      xk = temp;
    }
    else
    {
      pk = -Ck * bk;
      xk = xk + pk;
    }
    r = (min-xk).norm();
    bk1 = g(xk);
    yk = bk1 - bk;
    bk = bk1;
    rho = 1./(pk.transpose()*yk);
    Ck = Ck - rho * pk * (yk.transpose() * Ck) - rho * (Ck * yk) * pk.transpose() + pow(rho, 2.) * pk * (yk.transpose() * (Ck * yk)) * pk.transpose() + rho * pk * pk.transpose();
    err = bk.norm();
    outfile << iter << " " << err << " " << r << "\n";
  }
  outfile.flush();
  outfile.close();
  cout << "Der minimierte Vektor ist\n\n" << xk << endl;
}



int main()
{
  //Funktion 1
  Vector2d x0(2);
  x0 << -1.,-1.;
  Matrix2d I;
  I << 1., 0., 0., 1.;
  double init_3 = f1(x0);
  Matrix2d C0_1 = hesse1(x0).inverse();
  Matrix2d C0_2 = hesse1(x0);
  C0_2(0,1) = 0.;
  C0_2(1,0) = 0.;
  C0_2(0,0) = 1./C0_2(0,0);
  C0_2(1,1) = 1./C0_2(1,1);

  //Ich habe hier durch f(x0) geteilt weil sonst die Konvrgenz ohne Liniensuchschritt nicht gegeben war.
  Matrix2d C0_3 = I / init_3;

  bfgs(f1, g1, x0, C0_1, 1e-5, "1");
  bfgs(f1, g1, x0, C0_2, 1e-5, "2");
  bfgs(f1, g1, x0, C0_3, 1e-5, "3");


  bfgs(f1, g1, x0, C0_1, 1e-5, "1_l", true);
  bfgs(f1, g1, x0, C0_2, 1e-5, "2_l", true);
  bfgs(f1, g1, x0, C0_3, 1e-5, "3_l", true);

  //Aufgabenteil d)

  x0 << 1.5,2.3;

  MatrixXd C0_2_1 = hesse2(x0, false);
  MatrixXd C0_2_2 = hesse2(x0, true);
  MatrixXd C0_2_3 = I * f2(x0);

  bfgs(f2, g2, x0, C0_2_1, 1e-5, "d_1_0");
  bfgs(f2, g2, x0, C0_2_2, 1e-5, "d_2_0");
  bfgs(f2, g2, x0, C0_2_3, 1e-5, "d_1_0");


  bfgs(f2, g2, x0, C0_2_1, 1e-5, "d_1_0_l", true);
  bfgs(f2, g2, x0, C0_2_2, 1e-5, "d_2_0_l", true);
  bfgs(f2, g2, x0, C0_2_3, 1e-5, "d_3_0_l", true);

  x0 << -1.7,-1.9;

  C0_2_1 = hesse2(x0, false);
  C0_2_2 = hesse2(x0, true);
  C0_2_3 = I * f2(x0);

  bfgs(f2, g2, x0, C0_2_1, 1e-5, "d_1_1");
  bfgs(f2, g2, x0, C0_2_2, 1e-5, "d_2_1");
  bfgs(f2, g2, x0, C0_2_3, 1e-5, "d_1_1");


  bfgs(f2, g2, x0, C0_2_1, 1e-5, "d_1_1_l", true);
  bfgs(f2, g2, x0, C0_2_2, 1e-5, "d_2_1_l", true);
  bfgs(f2, g2, x0, C0_2_3, 1e-5, "d_3_1_l", true);

  x0 << 0.5,0.6;

  C0_2_1 = hesse2(x0, false);
  C0_2_2 = hesse2(x0, true);
  C0_2_3 = I * f2(x0);

  bfgs(f2, g2, x0, C0_2_1, 1e-5, "d_1_2");
  bfgs(f2, g2, x0, C0_2_2, 1e-5, "d_2_2");
  bfgs(f2, g2, x0, C0_2_3, 1e-5, "d_1_2");


  bfgs(f2, g2, x0, C0_2_1, 1e-5, "d_1_2_l", true);
  bfgs(f2, g2, x0, C0_2_2, 1e-5, "d_2_2_l", true);
  bfgs(f2, g2, x0, C0_2_3, 1e-5, "d_3_2_l", true);

  return 0;
}
