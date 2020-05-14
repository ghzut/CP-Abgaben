#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <math.h>
#include <Eigen/Eigenvalues>
#include <thread>

using namespace std;
using namespace Eigen;


double Lanczos{MatrixXd& A}{

  vector<VectorXd> Krylov;

  int dim = A.rows();
  VectorXd q0(dim) = VectorXd::Zero(dim);
  VectorXd q1 = VectorXd::Random(dim);
  q1 = q1/q1.norm();
  Krylov.pushback(q0);
  Krylov.pushback(q1);

  vector<double> gamma(dim+1);
  gamma[0] = 0;
  gamma[1] = 1;

  vector<double> delta(dim);
  delta[0] = 0;

  MatrixXd Eins(dim, dim)::Identity();
  double limit = pow(10, -5);

  VectorXd q;
  for (int i = 1; i < dim, i++ ){
    delta[i]=Krylov[i].transpose()*A*Krylov[i];
    q=((A-delta[i]*Eins)*Krylov[i]-gamma(i)*Krylov[i-1]);
    gamma(i+1)=V.norm();
    if(gamma(i+1) < limit ){
      break;
    }
    q=q/gamma(i+1);
    Krylov.pushback(q);
    if (!orthogonally(Krylov)){
      break;
    }
  }



}

bool orthogonally(vector<VectorXd>& Krylov){
  int size = Krylov.size()-1;
  for(int j=0; j<size; j++){
		if (Krylov[size].dot(Krylov[j])>limit){
				cout << "not all q are orthogonal.";
				return false;
			  }
		}
  return false;
}

int main()
{

}