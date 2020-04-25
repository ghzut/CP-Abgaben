#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>
#include<Eigen/Dense>
using namespace std;

typedef std::vector<vector<double>> double_vec;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;

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

  VectorXd vec_x(v_v_var.at(0).size());
  VectorXd vec_y(v_v_var.at(1).size());
  for (int i = 0; i < v_v_var.at(0).size();++i)
  {

    vec_x(i,0) = v_v_var.at(0).at(i);
    vec_y(i,0) = v_v_var.at(1).at(i);
  }

  cout << "The vector is:\n" << vec_x << endl;

  return 0;
}
