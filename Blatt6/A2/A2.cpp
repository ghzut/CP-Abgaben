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
    grad(0) = 2. * x(0) * (1. + 200.* (x(1)- pow(x(0),2.))) - 1.;
    grad(1) = 200 * (x(1) - pow(x(0), 2.));
    return grad;
  }
}



double f1_lambda(const VectorXd &x0, const VectorXd &b0, double l0)
{
  if (x0.size() != 2)
  {
    cerr << "Dieser Vektor ist nicht geeignet für diese Funktion." << endl;
    return -1.;
  }
  return pow(1.-x0(0)-l0*b0(0), 2.) + 100*pow(-pow(b0(0)*l0,2.)-(2*x0(0)*b0(0)-b0(1))*l0 + x0(1) - pow(b0(0),2.), 2.);
}


double erste_ableitung(function<double(const VectorXd&, const VectorXd&, double)> f, double x, const VectorXd &x0, const VectorXd &b0)
{
    double h = 0.001;

    return (f(x0, b0, x+h) - f(x0, b0, x-h))/(2*h);
}

double zweite_ableitung(function<double(const VectorXd&, const VectorXd&, double)> f, double x, const VectorXd &x0, const VectorXd &b0)
{
    double h = 0.001;

    return (f(x0, b0, x+h) - 2*f(x0, b0, x) + f(x0, b0, x-h))/(h*h);
}

VectorXd newton(function<double(const VectorXd&, const VectorXd&, double)> f, const VectorXd &x0, const VectorXd &b0)
{
  double dx;
  double l_0 = 4.;
  VectorXd x_new(x0.size());
  dx = erste_ableitung(f, l_0, x0, b0)/zweite_ableitung(f, l_0, x0, b0);
  x_new = x0 + (l_0+dx) *b0;
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
  double rho;
  VectorXd bk = g(x0);
  VectorXd xk = newton(f1_lambda, x0, bk);
  VectorXd bk1 = g(xk);
  MatrixXd Ck = C0;
  int iter = 1;
  bk = bk1;
  err = bk.norm();
  outfile << iter << " " << err << "\n";
  while (err > epsilon && iter < 50)
  {
    ++iter;
    pk = Ck * bk;
    xk += pk;
    bk1 = g(xk);
    yk = bk1 - bk;
    bk = bk1;
    rho = 1./(pk.transpose()*yk);
    Ck = Ck - rho*(Ck*yk)*pk.transpose() + pk*(yk.transpose()*Ck) + pk*pk.transpose() * pow(rho,2.)* (yk.transpose()*(Ck*yk)) + rho*pk*pk.transpose();
    /*if(err < 100*bk.norm())
    {
      cout << "Nope, so nicht. Fehler wird deutlich größer." << endl;
      break;
    }*/
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
  MatrixXd C0_3 = init_3 * I;
  bfgs(f1, g1, x0, C0_3, 1e-5, "3");
  return 0;
}
