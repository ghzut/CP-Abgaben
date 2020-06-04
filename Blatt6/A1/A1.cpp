#include <iostream>
#include <Eigen/Dense>
#include <fstream>

using namespace std;
using namespace Eigen;

//Definition der Rosenbrock-Funktion
double rosenbrock(double x1, double x2)
{
    return (1-x1)*(1-x1) + 100*(x2-x1*x1)*(x2-x1*x1);
}

//Funktion aus Aufgabe b)
double f(double x1, double x2)
{
    return 1/(1 + exp(-10*(x1*x2-3)*(x1*x2-3)/(x1*x1+x2*x2)));
}

double erste_ableitung(double (*f)(double, double), double x1, double x2, int k)
{
    double h = 0.001;

    if (k==0)
    {
        return (f(x1+h, x2) - f(x1-h, x2))/(2*h);
    }
    else
    {
        return (f(x1, x2+h) - f(x1, x2-h))/(2*h);
    }   
}

double intervallhalbierung(double (*min)(double, double, double, double, double, double (*f)(double, double)), 
                           double (*f)(double, double), double x1, double x2, double g1, double g2, 
                           double a, double b, double c, double x_c)
{
    double d;
    int n = 0;

    //ofstream file;
    //file.open("build/intervallhalbierung.txt");
    //file << "# n      a       b       c\n\n";
    //file << n << "  " << a << "  " << b << "  " << c << "\n";

    
    while (c-a > x_c)
    {
        if (b-a < c-b)
        {
            d = (c+b)/2;

            if (min(d, x1, x2, g1, g2, f) < min(b, x1, x2, g1, g2, f))
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

            if (min(d, x1, x2, g1, g2, f) < min(b, x1, x2, g1, g2, f))
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

        //file << n << "  " << a << "  " << b << "  " << c << "\n";
    }

    //file.close();

    return (c-a)/2;
}

double minimize(double lambda, double x1, double x2, double g1, double g2, double (*f)(double, double))
{
  return f(x1 + lambda*g1, x2 + lambda*g2);
}

void gradientenverfahren(double (*f)(double, double), Vector2d x_0)
{
    VectorXd g(x_0.size());
    double lambda;
    int n=0;

    g(0) = erste_ableitung(f, x_0(0), x_0(1), 0)*(-1);
    g(1) = erste_ableitung(f, x_0(0), x_0(1), 1)*(-1);

    ofstream file;
    file.open("build/gradientenverfahren.txt");
    file << "# n      x1       x2       g1      g2\n\n";
    file << n << "  " << x_0(0) << "  " << x_0(1) << "  " << g(0) << "  " << g(1) << "\n";


    while (g.norm() > 0.005)
    {

    //for(int i = 0; i<1500000; i++){
        lambda = intervallhalbierung(minimize, rosenbrock, x_0(0), x_0(1), g(0), g(1), -50, 0, 50, 1e-6);
        x_0 = x_0 + lambda * g;

        n += 1;

        g(0) = erste_ableitung(f, x_0(0), x_0(1), 0)*(-1);
        g(1) = erste_ableitung(f, x_0(0), x_0(1), 1)*(-1);

        if (n%10000==0)
        {
            file << n << "  " << x_0(0) << "  " << x_0(1) << "  " << g(0) << "  " << g(1) << "\n";
        }
    } 
  
    file << n << "  " << x_0(0) << "  " << x_0(1) << "  " << g(0) << "  " << g(1) << "\n";

    file.close();
  
}

void konjugiert(double (*f)(double, double), Vector2d x_0)
{
    Vector2d g, g_neu, p;
    double lambda, mu;
    int n=0;

    g(0) = erste_ableitung(f, x_0(0), x_0(1), 0)*(-1);
    g(1) = erste_ableitung(f, x_0(0), x_0(1), 1)*(-1);

    p = g;

    ofstream file;
    file.open("build/konjugiert.txt");
    file << "# n      x1       x2       g1      g2\n\n";
    file << n << "  " << x_0(0) << "  " << x_0(1) << "  " << g(0) << "  " << g(1) << "\n";


    while (g.norm() > 0.005)
    {

    //for(int i = 0; i<1000000; i++){
        lambda = intervallhalbierung(minimize, rosenbrock, x_0(0), x_0(1), p(0), p(1), -50, 0, 50, 1e-6);
        x_0 = x_0 + lambda * p;

        n += 1;

        g_neu(0) = erste_ableitung(f, x_0(0), x_0(1), 0)*(-1);
        g_neu(1) = erste_ableitung(f, x_0(0), x_0(1), 1)*(-1);

        mu = (g_neu.dot(g_neu))/(g.dot(g));
        g = g_neu;
        p = g + mu * p;

        //if (n%10000==0)
        //{
            file << n << "  " << x_0(0) << "  " << x_0(1) << "  " << g(0) << "  " << g(1) << "\n";
        //}
        
    } 

    // file << n << "  " << x_0(0) << "  " << x_0(1) << "  " << g(0) << "  " << g(1) << "\n";
  

    file.close();
}

int main()
{   

    Vector2d x_0;
    x_0 << -1,-1;

    gradientenverfahren(&rosenbrock, x_0);
    konjugiert(&rosenbrock, x_0);

    return 0;
}