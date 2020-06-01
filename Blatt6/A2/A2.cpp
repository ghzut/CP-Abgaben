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
  VectorXd bk1(n);
  VectorXd bk = g(x0);
  VectorXd pk = C0 * bk;
  VectorXd xk = x0 + pk;
  VectorXd yk = bk;
  double rho = 1./(pk.transpose()*yk);
  MatrixXd Ck = (C0*yk)*pk.transpose() + pk*(yk.transpose()*C0) - (yk.transpose()*(C0*yk)*rho)*pk*pk.transpose() - pk*pk.transpose();
  Ck *= rho;
  Ck = C0 - Ck;
  err = bk.norm();
  int iter = 0;
  while (err > epsilon)
  {
    ++iter;
    bk1 = g(xk);
    pk = Ck * bk1;
    xk += pk;
    yk = bk1 - bk;
    bk = bk1;
    rho = 1./(pk.transpose()*yk);
    Ck = Ck/rho - (Ck*yk)*pk.transpose() + pk*(yk.transpose()*Ck) + (yk.transpose()*(Ck*yk)*rho)*pk*pk.transpose() + pk*pk.transpose();
    Ck *= rho;
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
