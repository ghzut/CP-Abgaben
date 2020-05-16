#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

//Zu integrierende Funktionen a)
double f1(double x)
{
    return exp(x)/x;
}


//Zu integrierende Funktionen b)
double f2(double x)
{
    return exp(-x)/sqrt(x);
}


//Zu integrierende Funktionen c)
double f3(double x)
{
    return sin(x)/x;
}

//Wir verwenden die Simpsonregel, da diese die wenigsten Iterationen benötigt
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

void integrate_exe(double (*f)(double), double a, double b, double max_err, string aufgabe)
{
  string outfilename = "build/A2"+aufgabe+".txt";
  ofstream outfile(outfilename, ofstream::trunc);
  outfile << "#n, int, err\n";
  double err, temp_result, new_result;

  err = 10000;
  temp_result = simpson(y, a, b, n);
  new_result = 0;

  while (err >= max_err)
  {
    n = 2*n;
    new_result = integrate(y, a, b, n);
    err = abs(temp_result-new_result)/temp_result;
    temp_result = new_result;
    outfile << n << " " << temp_result << " " << err << "\n";
  }
  outfile.flush();
  outfile.close();
}

int main()
{
  integrate_exe(&f1, -1., 1., 1e-7, "a");

  integrate_exe(&f1, 0., 1e, 1e-8, "b");

  integrate_exe(&f1, -1., 1., 1e-7, "c");

  return 0;
}
