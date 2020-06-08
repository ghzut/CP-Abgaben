#include <iostream>
#include <Eigen/Dense>
#include <fstream>

using namespace std;
using namespace Eigen;

double f1(double x, double x_s, double y_s, double z_s)
{
    return 1/sqrt((x-x_s)*(x-x_s)+y_s*y_s+z_s*z_s);
}

double f2(double x, double x_s, double y_s, double z_s)
{
    return x/sqrt((x-x_s)*(x-x_s)+y_s*y_s+z_s*z_s);
}

double integrate1D(double x, double (*f)(double, double, double, double), double a, double b, int n, double y_s, double z_s)
{
    double h, result; 
    h = (b-a)/n;
    result = 0;

    for (int k = 1; k<n+1; k++)
    {
        result += f(x, a-h/2 + k*h, y_s, z_s);
    }

    result = result * h;

    return result;
}

double integrate2D(double x, double (*f)(double, double, double, double), double a, double b, int n, VectorXd y_s, double z_s)
{
    double h, result;
    h = (b-a)/n;
    result = 0;

    for (int k = 1; k<n+1; k++)
    {
        result += integrate1D(x, f, a, b, n, y_s(k-1), z_s);
    }

    result = result * h;

    return result;
}
 
double integrate3D(double x, double (*f)(double, double, double, double), double a, double b, int n)
{
    double h, result;
    h = (b-a)/n;
    result = 0;

    VectorXd y_s = VectorXd::LinSpaced(n, a, b);
    VectorXd z_s = y_s;

    for (int k = 1; k<n+1; k++)
    {
        result += integrate2D(x, f, a, b, n, y_s, z_s(k-1));
    }

    result = result * h;

    return result;
}

int main()
{

    VectorXd x = 0.1*VectorXd::LinSpaced(70, 11, 80);
    double phi=0, phi_b=0;



    ofstream file_1, file_2;
    file_1.open("build/ausserhalb_a.txt");
    file_2.open("build/ausserhalb_b.txt");

    file_1 << "# x      phi \n\n";
    file_2 << "# x      phi \n\n";

    for (int i=0; i<x.size(); i++)
    {
        phi = integrate3D(x(i), &f1, -1, 1, 100);
        phi_b = integrate3D(x(i), &f2, -1, 1, 100);
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
        phi = integrate3D(x(i), &f1, -1, 1, 100);
        phi_b = integrate3D(x(i), &f2, -1, 1, 100);
        file_1 << x(i) << "       " << phi << "\n";
        file_2 << x(i) << "       " << phi_b << "\n";
    }

    file_1.close();
    file_2.close();

    return 0;
}