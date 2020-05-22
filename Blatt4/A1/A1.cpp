#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>

using namespace std;
//Zu integrierende Funktionen a)
double f1(double x)
{
    return (exp(x)-exp(-x))/x;
}


//Zu integrierende Funktionen b)
double f2(double x)
{
    return 2*exp(-pow(x,2));
}


//Zu integrierende Funktionen c)
double f3(double x)
{
    return 2*(sin(x))/x;
}

//Um Polstellen zu umgehen wird ein offenes Verfahren benötigt
double mittelpunkt(double (*f)(double), double a, double b, int n)
{
    double h, result;
    h = (b-a)/n;
    result = 0;

    for (int k = 1; k<n+1; k++)
    {
        result += f(a-h/2. + k*h);
    }

    result = result * h;

    return result;
}

//Wir verwenden die Simpsonregel für Aufgabe b, auf Grund der schnellen Konvergenz
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
  result = f(x[0]);
	for(j=1; j<n/2; j++)
	{
		result += 4*f(x[2*j - 1]) + 2*f(x[2*j]);
	}
  result += 4*f(x[n-1]) + f(x[n]);
	return result*h/3;
}

double get_Int(double (*integrate)(double (*f)(double), double, double, int), double (*func)(double), double a, double max_err, double b)
{
  double temp, new_res;
  double err = 10000.;
  int n=2;
  temp = integrate(func, a, b, n);
  while (err > max_err)
  {
    n = 2*n;
    new_res = integrate(func, a, b, n);
    err = abs(temp-new_res);
    temp = new_res;
  }
  return temp;
}



void int_a(double (*func)(double),double a, double max_err, double limit, string part)
{
  ofstream outfile("build/A1"+part+".txt", ofstream::trunc);
  outfile << "#result\n";
  outfile.precision(8);
  double result;
  result = get_Int(&mittelpunkt, func, a, max_err, limit);
  outfile << result << endl;
  outfile.close();
}

void int_b_c(double (*integrate)(double (*f)(double), double, double, int),double (*func)(double), double a, double max_err, double limit, double increment, string part)
{
  ofstream outfile("build/A1"+part+".txt", ofstream::trunc);
  outfile << "#i, int, int 2, err\n";
  outfile.precision(9);
  double result, result2;

  for (double i = 1.; i < limit; i*=increment)
  {
    result = get_Int(integrate, func, a, max_err, i);
    result2 = get_Int(integrate, func, a, max_err, i+1.);//increment/2.);
    outfile << i << " " << result << " " << result2 << " " << abs(result2-result) << "\n";
  }
  outfile.flush();
  outfile.close();
}
int main()
{
  cout << setprecision(10);
  int_a(&f1, 0., 1e-7, 1., "a");
  int_b_c(&simpson, &f2, 0., 1e-10, pow(2.,5.), 2., "b");
  int_b_c(&mittelpunkt, &f3, 0., 1e-7, pow(5.,10.), 5., "c");
  return 0;
}
