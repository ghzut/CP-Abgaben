#include<iostream>
#include<fstream>
using namespace std;

typedef vector<vector<double>> double_vec;

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
  inputdata(data_path.c_str(), v_vars)

  for (int i = 0; i < v_var_str.size(), ++i)
  {
    stringstream ss(v_var_str.at(i))
    vector<double> v_var;
    while (ss >> var)
    {
      v_var.push_back(var);
    }
    v_v_var.push_back(v_var);
    v_var.clear();
  }
  

  return 0;
}
