#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>

using namespace std;
//Zu integrierende Funktionen a)
long double f1(long double x)
{
    return exp(x)/x;
}


//Zu integrierende Funktionen b)
long double f2(long double x)
{
    return 2*exp(-pow(x,2));
}


//Zu integrierende Funktionen c)
long double f3(long double x)
{
    return sin(x)/x;
}

//Wir verwenden die Simpsonregel, da diese die wenigsten Iterationen benötigt
long double simpson(long double (*f)(long double), long double a, long double b, int n)
{
	long double h, x[n+1], result = 0;
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

long double get_Int(long double (*f)(long double), long double a, long double max_err, long double b)
{
  long double temp, new_res;
  long double err = 10000.;
  long double n=2.;
  temp = simpson(f, a, b, n);
  while (err > max_err)
  {
    n = 2*n;
    new_res = simpson(f, a, b, n);
    err = abs(temp-new_res);
    temp = new_res;
    cout << temp << " " << err << endl;
  }
  cout << endl;
  return temp;
}

void integrate_b(long double a, long double max_err, long double limit)
{
  string outfilename = "build/A2b.txt";
  ofstream outfile(outfilename, ofstream::trunc);
  outfile << "#i, int, err\n";
  long double result, result2;

  for (long double i = 1.; i < limit; ++i)
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
  cout << fixed;
  cout << setprecision(8);
  integrate_b(0., 0.1, 20.);

  return 0;
}
