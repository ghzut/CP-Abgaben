#include "cube.h"

//Wenn ich die Funktionen im Header inkludiere kriege ich immer "invalid use of non-static member function" als Fehlermeldung wenn ich sie später benutzen will, daher stehen sie erstmal hier. Ich habe die Lösung nicht rausgefunden, sonst ständen sie im header :(
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