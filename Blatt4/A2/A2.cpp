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
    return exp(-pow(x,2));
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

long double get_Int(double (*f)(double), double a, double max_err, double b)
{
  long double temp, new_res;
  double err = 10000.;
  double n=2.;
  temp = simpson(f, a, b, n);
  /*while (err >= max_err)
  {
    n = 2*n;
    new_res = simpson(f, a, b, n);
    err = abs(temp-new_res);
    temp = new_res;
  }*/
  return temp;
}

void integrate_b(double a, double max_err, double limit)
{
  string outfilename = "build/A2b.txt";
  ofstream outfile(outfilename, ofstream::trunc);
  outfile << "#i, int, err\n";
  long double result, result2;

  for (double i = 10.; i < limit; i*=10)
  {
    result = get_Int(&f2, a, max_err, limit);
    result2 = get_Int(&f2, a, max_err, 2*limit);
    outfile << i << " " << result2 << " " << abs(result2-result) << "\n";
  }
  outfile.flush();
  outfile.close();
}
int main()
{
  integrate_b(0., 1e-4, 1e7);

  return 0;
}
