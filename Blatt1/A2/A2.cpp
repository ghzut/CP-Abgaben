#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>
#include<Eigen/Dense>
using namespace std;

typedef std::vector<vector<double>> double_vec;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;
typedef Eigen::Matrix<double, 2, 2> Matrix2d;


void inputdata(const char *path, vector<string> &v_line)
{
  ifstream datafile(path);
  string line;
  while (getline(datafile,line))
  {
    v_line.push_back(line);
  }
}

int main()
{
  double_vec v_v_var;
  vector<string> v_var_str;
  double var;
  string data_path = "data/data_x_y.txt";
  inputdata(data_path.c_str(), v_var_str);

  for (int i = 0; i < v_var_str.size(); ++i)
  {
    stringstream ss(v_var_str.at(i));
    vector<double> v_var;
    while (ss >> var)
    {
      v_var.push_back(var);
    }
    v_v_var.push_back(v_var);
  }
  int size = v_v_var.at(0).size();
  VectorXd vec_x(size);
  VectorXd vec_y(size);
  for (int i = 0; i < size; ++i)
  {
    vec_x(i,0) = v_v_var.at(0).at(i);
    vec_y(i,0) = v_v_var.at(1).at(i);
  }
  Eigen::Matrix<double, 10, 2> M_x_y;
  M_x_y.col(0) = vec_x;
  M_x_y.col(1) = vec_y;

  Eigen::Matrix<double, 10, 2> A;

  A.col(0) = vec_x;
  for (int i = 0; i < size; ++i)
  {
    A(i,1)=1.;
  }
  Matrix2d n = A.transpose()*A;
  Matrix2d L;
  Matrix2d U;
  Matrix2d P;
  VectorXd b(2); b = A.transpose() * vec_y;
  VectorXd x(2);

  cout << "The matrix:\n" << A << endl;
  cout << "The matrix transpose:\n" << A.transpose() << endl;
  cout << "The quadratic matrix:\n" << n << endl;
  cout << "The new b vector is:\n" << b <<endl;

  x = n.partialPivLu().solve(b);
  L = n.partialPivLu().matrixLU().triangularView<Eigen::UpLoType::UnitLower>();
  U = n.partialPivLu().matrixLU().triangularView<Eigen::UpLoType::Upper>();
  P = n.partialPivLu().permutationP();

  cout << "x =\n" << x << "\n\n";
  cout << "L =\n" << L << "\n\n";
  cout << "U =\n" << U << "\n\n";
  cout << "P =\n" << P << "\n\n";

  ofstream outfile("data/py_input_m_n.txt", ofstream::trunc);
  outfile << "#m, n\n";
  outfile << x;
  outfile.flush();
  outfile.close();

  ofstream outfile2("data/py_input_x_y.txt", ofstream::trunc);
  outfile2 << "#x, y\n";
  outfile2 << M_x_y;
  outfile2.flush();
  outfile2.close();
  return 0;
}
