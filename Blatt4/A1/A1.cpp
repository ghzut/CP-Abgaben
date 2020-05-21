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
    return 2*(sin(x)+sin(1./x))/x;
}

//Um Polstellen zu umgehen wird ein offenes Verfahren benötigt
double mittelpunkt(double (*f)(double), double a, double b, int n)
{
    double h, result;
    h = (b-a)/n;
    result = 0;

    for (int k = 1; k<n+1; k++)
    {
        result += 2*f(a-h/4. + k*h);
        result -= f(a-2.*h/4. + k*h);
        result += 2*f(a-3.*h/4. + k*h);
    }

    result = result * 1./3.* h;

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
    if(abs(temp-new_res) < err)
    {
    err = abs(temp-new_res);
    temp = new_res;
    cout << err << endl;
    }
    else continue;
  }
  return temp;
}



void int_a_c(double (*func)(double),double a, double max_err, double limit, string part)
{
  ofstream outfile("build/A2"+part+".txt", ofstream::trunc);
  outfile << "#result\n";
  outfile.precision(8);
  double result;
  result = get_Int(&mittelpunkt, func, a, max_err, limit);
  outfile << result << endl;
  outfile.close();
}

void int_b(double a, double max_err, double limit)
{
  ofstream outfile("build/A2b.txt", ofstream::trunc);
  outfile << "#i, int, err\n";
  outfile.precision(9);
  double result, result2;

  for (double i = 1.; i < limit; i*=2.)
  {
    result = get_Int(&simpson, &f2, a, max_err, i);
    result2 = get_Int(&simpson, &f2, a, max_err, i+1.);
    outfile << i << " " << result << " " << abs(result2-result) << "\n";
  }
  outfile.flush();
  outfile.close();
}
int main()
{
  int_a_c(&f1, 0., 1e-7, 1., "a");
  cout << endl << endl;
  int_b(0., 1e-10, pow(2.,7.));
  //cout << endl << endl;
  //int_a_c(&f3, 0., 1e-7, 1., "c");
  return 0;
}
