#include <iostream>
#include <Eigen/Dense>
#include <fstream>

using namespace std;
using namespace Eigen;


class Cube {
    public:

    double a;
    double b;
    int n;
    
    double h(){
        return (b-a)/n;
    } 

    double arg (int k) {
        return a-h()/2 + k*h();
    }

    double integrate1D(const double& x, double (*f)(const double&, const double&, const double&, const double&), double y_s,    double z_s){
        double result = 0;

        for (int k = 1; k <= n; ++k){
            result += f(x, arg(k), y_s, z_s);
        }
        result = result * h();
        return result;
    }

    double integrate2D(const double x, double (*f)(const double&, const double&, const double&, const double&), double z_s){
        double result = 0;

        for (int k = 1; k <= n; ++k){
            result += integrate1D(x, f, arg(k), z_s);
        }
        result = result * h();
        return result;
    }

    double integrate3D(const double& x, double (*f)(const double&, const double&, const double&, const double&)){
        double result = 0;

        for (int k = 1; k <= n; ++k){
            result += integrate2D(x, f, arg(k));
        }
        result = result * h();
        return result;
    }

};

//Eigentlich sollten die beiden Funktionen Teil der Klasse sein...wenn ich aber später versuche via cube.f1 bzw cube.f2 darauf zuzugreifen (z.B. Zeile 81)) , kriege ich einen "invalid-use-of-non-static-member-function"-error. Ich habe es jetzt nicht anders behoben bekommen, kennt ihr eine gute Lösung?

double f1(const double& x, const double& x_s, const double& y_s, const double& z_s){
        return 1/sqrt((x-x_s)*(x-x_s)+y_s*y_s+z_s*z_s);
}

double f2(const double& x, const double& x_s, const double& y_s, const double& z_s){
            return x/sqrt((x-x_s)*(x-x_s)+y_s*y_s+z_s*z_s);
}


int main() 
{
    VectorXd x = 0.1*VectorXd::LinSpaced(70, 11, 80);
    double phi=0, phi_b=0; 

    Cube cube; 
    cube.a = -1;
    cube.b = 1;
    cube.n = 100;

    ofstream file_1, file_2;
    file_1.open("build/ausserhalb_a.txt");
    file_2.open("build/ausserhalb_b.txt");

    file_1 << "# x      phi \n\n";
    file_2 << "# x      phi \n\n";

    for (int i=0; i<x.size(); i++)
    {
        phi = cube.integrate3D(x(i), f1);
        phi_b = cube.integrate3D(x(i), f2);
        file_1 << x(i) << "       " << phi << "\n";
        file_2 << x(i) << "       " << phi_b << "\n";
    }

    file_1.close();
    file_2.close();

    file_1.open("build/innerhalb_a.txt");
    file_2.open("build/innerhalb_b.txt");

    file_1 << "# x      phi \n\n";
    file_2 << "# x      phi \n\n";

    x = 0.1*VectorXd::LinSpaced(11, 0, 10);

    for (int i=0; i<x.size(); i++)
    {
        phi = cube.integrate3D(x(i), f1);
        phi_b = cube.integrate3D(x(i), f2);
        file_1 << x(i) << "       " << phi << "\n";
        file_2 << x(i) << "       " << phi_b << "\n";
    }

    file_1.close();
    file_2.close();

    return 0;
}