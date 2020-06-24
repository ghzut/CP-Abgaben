#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <stdlib.h>

using namespace std;
using namespace Eigen;

double PotentialLJ(const Vector2d &r, double rc)
{
    double r2 = r.squaredNorm();
    double r6 = 1/pow(r2,3);
    double pot;
    if (r2 > rc*rc){
        pot = 0;
    }
    else {
        pot = 4 * (r6*r6 - r6);
    }
    return pot;
}

Vector2d KraftLJ(const Vector2d &r, double rc)
{
    double r2 = r.squaredNorm();
    double r6 = 1./pow(r2,3);
    Vector2d kraft;

    if (r2 > rc*rc)
    {
        kraft = Vector2d::Zero();
    }
    else 
    {
            kraft = 48*r*((r6*r6)/r2 - 0.5*r6/r2);
    }

    return kraft
}

void periodische_RB(vector<Vector2d> &r, double L)
{
    for(uint i = 0; i < r.size(); ++i)
    {
      for(uint j = 0; j < r.at(i).size(); ++j)
      {
        if(r.at(i)(j) > L || r.at(i)(j) < 0) 
        {
            r.at(i)(j) -= floor(r.at(i)(j)/L) * L;
        }
      }
    }
}

Vector2d v_center(const vector<Vector2d> &v, int N)
{
    Vector2d vs = Vector2d::Zero();
    for(int i = 0; i < v.size(); ++i)
    {
      vs += v.at(i);
    }
    return vs/N;
}

double ekin(const vector<Vector2d> &v)
{
    double m = 1.;
    double ekin = 0.;
    for(uint i = 0; i < v.size(); ++i)
    {
      ekin += m*v.at(i).squaredNorm();
    }

    return ekin/2.;
}

double T (const vector<Vector2d> &v, uint N)
{
    double Nf = 2N-2;
    double T = 2*ekin(v)/Nf;
    return T;
}



auto init(double L, int N, double particlesPerRow, double T)
{
    vector<Vector2d> r, v;
    
    Vector2d r_vec;
    for(int n = 0; n < particlesPerRow; ++n)
    {
      for(int m = 0; m < particlesPerRow; ++m)
      {
        r_vec << 1 + 2 * n, 1 + 2 * m;
        r.push_back(r_vec);
      }
    }

    srand(42);
    Vector2d v_vec;
    for(int i = 0; i < N; ++i)
    {
      v_vec << double(rand()%10), double(rand()%10);
      v.push_back(v_vec);
    }

    Vector2d v_s = v_center(v, N);
    for (Vector2d& n : v)
    {
        n -= v_s;
    }

    double Nf = 2*N-2;                              // Anzahl Freiheitsgrade
    double skal = T*Nf/(2*ekin(v));
    for (int n=0; n<N; n++)
    {
        v.at(n) = skal*v.at(n);
    }
    struct vecs{
        vector<Vector2d> Ort, Geschwindigkeit;
    };
    return vecs{r, v};
}

vector<Vector2d> verlet_r(Vector2d (*f)(Vector2d &, double &), double (*pot)(const Vector2d &, double), const vector<Vector2d> &r_n, vector<Vector2d> r_nminus1, double h, int N, double L, double &Epot)
{
    vector<Vector2d> r_nplus1;
    Vector2d r_nplus1_vec, a;
    a = beschleunigung(f, pot, r_n, N, L, Epot );

    for (int i = 0; i < N; i++)
    {
        r_nplus1_vec = 2*r_n.at(i) - r_nminus1.at(i) + a.at(i)*h*h;
        r_nplus1.push_back(r_nplus1_vec);
    }

    return r_nplus1;
}

vector<Vector2d> verlet_v(vector<Vector2d> r_nminus1, vector<Vector2d> r_nplus1, double h, int N)
{
    vector<Vector2d> v_n;
    Vector2d v_n_vec;

    for (int i = 0; i < N; i++)
    {
        v_n_vec = (r_nplus1.at(i) - r_nminus1.at(i))/(2*h);
        v_n.push_back(v_n_vec);
    }

    return v_n;
}

vector<Vector2d> Beschleunigung(Vector2d (*f)(Vector2d &, double &), double (*pot)(const Vector2d &, double), const vector<Vector2d> &r_n, int N, double L, double &Epot)
{   
    double rc = L/2;
    Vector2d xx;
    xx << 0,0;

    vector<Vector2d> a(N, xx);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i != j)
            {
                Vector2d shift = Vector2d::Zero();

                for (int k = -1; k <= 1; k++)
                {
                    shift(0) = k*L;

                    for (int l = -1; l <= 1; l++)
                    {
                        shift(1) = l*L;

                        Vector2d temp=r.at(i)-(r.at(j)+shift);
						a.at(i)+=f(rel, rc);

                        if(i < j){
                            Epot += pot(rel, rc);
                        }
                    }
                }
            }
        }
    }
    return a;
}

void Aequilibrierung(Vector2d (*f)(const Vector2d &, double), double (*pot)(const Vector2d &, double), int N, double L, double particlesPerRow, double T, Vector2d, uint steps, double h ){
    auto rn, vn = init (L, N, particlesPerRow, T);
    vector<Vector2d> r_minus1, an;
    Vector2d r_temp;
    double Epot = 0;
    an = Beschleunigung(f, pot, rn, N, L, Epot );

    for(int i = 0; i < N; ++i){
        r_temp = rn.at(i)-vn.at(i)*h + 0.5 * an.at(i)*h*h;
        r_minus1.push_back(r_temp);

    }

}




int main()
{


    return 0;
}