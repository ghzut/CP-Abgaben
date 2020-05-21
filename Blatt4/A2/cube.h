#include <iostream>
#include <Eigen/Dense>
#include <fstream>

using namespace std;
using namespace Eigen;

class Cube {
    private:
    	double a = 1;
        double b = -1;
        int n = 100;
        double h = (b-a)/n;

    public:

    double arg (int k) {
        return a-h/2 + k*h;
    }

    double integrate1D(const double& x, double (*f)(const double&, const double&, const double&, const double&), double y_s,    double z_s){
        double result = 0;

        for (int k = 1; k <= n; ++k){
            result += f(x, arg(k), y_s, z_s);
        }
        result = result * h;
        return result;
    }

    double integrate2D(const double x, double (*f)(const double&, const double&, const double&, const double&), double z_s){
        double result = 0;

        for (int k = 1; k <= n; ++k){
            result += integrate1D(x, f, arg(k), z_s);
        }
        result = result * h;
        return result;
    }

    double integrate3D(const double& x, double (*f)(const double&, const double&, const double&, const double&)){
        double result = 0;

        for (int k = 1; k <= n; ++k){
            result += integrate2D(x, f, arg(k));
        }
        result = result * h;
        return result;
    }

};

