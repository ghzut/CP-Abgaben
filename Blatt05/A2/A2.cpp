#include <iostream>
#include <Eigen/Dense>
#include <fstream>

using namespace std;
using namespace Eigen;


//Funktion, die minimiert werden soll
double f(double x)
{
    return x*x-2;
}


//Funktionen zur numerischen Bestimmung der ersten und zweiten Ableitung einer Funktion
double erste_ableitung(double (*f)(double), double x)
{
    double h = 0.001;

    return (f(x+h) - f(x-h))/(2*h);
}

double zweite_ableitung(double (*f)(double), double x)
{
    double h = 0.001;

    return (f(x+h) - 2*f(x) + f(x-h))/(h*h);
}


//Minimierungsroutinen
void intervallhalbierung(double (*f)(double), double a, double b, double c, double x_c)
{
    double d;
    int n = 0;

    ofstream file;
    file.open("build/intervallhalbierung.txt");
    file << "# n      a       b       c\n\n";
    file << n << "  " << a << "  " << b << "  " << c << "\n";

    
    while (c-a > x_c)
    {
        if (b-a < c-b)
        {
            d = (c+b)/2;

            if (f(d) < f(b))
            {
                a = b;
                b = d;
            }
            else
            {
                c = d;
            }
        }
        else
        {
            d = (a+b)/2;

            if (f(d) < f(b))
            {
                c = b;
                b = d;
            }
            else
            {
                a = d;
            }
            
        }

        n += 1;

        file << n << "  " << a << "  " << b << "  " << c << "\n";
    }

    file.close();
}

void newton(double (*f)(double), double x_0, double x_c)
{
    double dx;
    int n = 0;

    dx = erste_ableitung(f, x_0)/zweite_ableitung(f, x_0);

    ofstream file;
    file.open("build/newton.txt");
    file << "# n      x_0\n\n";
    file << n << "  " << x_0 << "\n";


    while (dx > x_c)
    {
        n += 1;
        x_0 = x_0 - dx;
        
        file << n << "  " << x_0 << "\n";

        dx = erste_ableitung(f, x_0)/zweite_ableitung(f, x_0);
    }
}



int main()
{   
    double a, b, c, x_0, x_c;

    x_c = 1e-9;
    a = -0.5;
    b = -0.1;
    c = 2;
    x_0 = 1;

    intervallhalbierung(&f, a, b, c, x_c);
    newton(&f, x_0, x_c);

    return 0;
}