#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;


double f1(double x)
{
    return exp(-x)/x;
}

double f2(double x)
{
    return x * sin(1/x);
}

double mittelpunkt(double (*f)(double), double a, double b, int n)
{
    double h, result;
    h = (b-a)/n;
    result = 0;

    for (int k = 1; k<n+1; k++)
    {
        result += f(a-h/2 + k*h);
    }

    result = result * h;

    return result;
}

double trapez(double (*f)(double), double a, double b, int n)
{
    double h, result;
    h = (b-a)/n;
    result = 0;

    //cout << result << "\n";

    for (int k = 1; k<n; k++)
    {
        result += f(a+k*h);
    }

    result = h * (1/2 * f(a) + 1/2 * f(b) + result);

    return result;
}


double simpson(double (*f)(double), double a, double b, int n)
{
	double h, x[n+1], result = 0;
	int j;
	h = (b-a)/n;
	
	x[0] = a;
	
	for(j=1; j<=n; j++)
	{
		x[j] = a + h*j;
	}
	
	for(j=1; j<=n/2; j++)
	{
		result += f(x[2*j - 2]) + 4*f(x[2*j - 1]) + f(x[2*j]);
	}
	
	return result*h/3;
}

void test(double (*integrate)(double (*f)(double), double, double, int), double (*y)(double), double a, double b, string filename)
{
    double err, result, new_result;
    int n = 2;

    ofstream file;
    file.open("build/"+filename+".txt");
    file << "# n            err \n\n"; 

    err = 1e6;
    result = integrate(y, a, b, n);
    new_result = 0;

    while (err >= 1e-4)
    {
        n = 2*n;
        new_result = integrate(y, a, b, n);
        err = abs(result-new_result)/result;
        result = new_result;

        file << n << "          " << err << "\n";        
    }

    file.flush();
    file.close();
}

int main()
{

    test(&mittelpunkt, &f1, 1, 100, "mittelpunktsregel_f1");
    test(&trapez, &f1, 1, 100, "trapezregel_f1");
    test(&simpson, &f1, 1, 100, "simpsonregel_f1");
    
    test(&mittelpunkt, &f2, 1, 100, "mittelpunktsregel_f2");
    test(&trapez, &f2, 1, 100, "trapezregel_f2");
    test(&simpson, &f2, 1, 100, "simpsonregel_f2");

    //cout << mittelpunkt(&f1, 1, 100, 10000000) << "\n";
    //cout << trapez(&f2, 0, 1, 10000000) << "\n";
    //cout << simpson(&f1, 1, 100, 10000000) << "\n";

    return 0;
}
