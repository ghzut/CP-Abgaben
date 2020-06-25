#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <stdlib.h>
#include <fstream>

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

    return kraft;
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

VectorXd naechsterNachbar(VectorXd r1, VectorXd r2, double L)
{
    double rc = L/2;
    VectorXd shift = VectorXd::Zero(2, 1);     // Verschiebungsvektor
    VectorXd diff = VectorXd::Zero(2, 1);   // Abstandsvektor zw. Teilchen und Bildteilchen

    // Pruefe zunaechst die Bildraeume

    // Bildraum rechts
    shift << L, 0;
    diff = r1-(r2+shift);
    if (diff.norm()<rc) {
        return r2+shift;
    }
    // Bildraum links
    diff = r1-(r2-shift);
    if (diff.norm()<rc) {
        return r2-shift;
    }

    // Bildraum oben
    shift << 0, L;
    diff = r1-(r2+shift);
    if (diff.norm()<rc) {
        return r2+shift;
    }
    // Bildraum unten
    diff = r1-(r2-shift);
    if (diff.norm()<rc) {
        return r2-shift;
    }

    // Bildraum oben rechts
    shift << L, L;
    diff = r1-(r2+shift);
    if (diff.norm()<rc) {
        return r2+shift;
    }
    // Bildraum unten links
    diff = r1-(r2-shift);
    if (diff.norm()<rc) {
        return r2-shift;
    }

    // Bildraum unten rechts
    shift << L, -L;
    diff = r1-(r2+shift);
    if (diff.norm()<rc) {
        return r2+shift;
    }
    // Bildraum oben links
    diff = r1-(r2-shift);
    if (diff.norm()<rc) {
        return r2-shift;
    }

    // Kein Bildteilchen nah genug, also muss das reale Teilchen am naechsten sein.
    return r2;
}

vector<Vector2d> update_kraft(const vector<Vector2d> &r, int N, double L)
{
    double rc = L/2;
    Vector2d xx;
    xx << 0,0;
    vector<Vector2d> F(N, xx);
    // Folgende Groessen als Zwischenspeicher benoetigt
    VectorXd delta_r = VectorXd::Zero(2, 1);     // Abstandsvektor
    VectorXd r_next = VectorXd::Zero(2, 1);     // nächster Nachbar
    VectorXd F_temp = VectorXd::Zero(2, 1);     // Kraft zwischen j und i

    for (int i=0; i<(N-1); i++)
    {
        for (int j=i+1; j<N; j++)   // Keine Selbstwecheelwirkung oder Doppelzaehlung
        {
            r_next = naechsterNachbar(r.at(i), r.at(j), L);
            delta_r = r.at(i)-r_next;
            if(delta_r.norm()<rc)
            {
                F_temp = KraftLJ(delta_r, L);
                F.at(i) = F.at(i) + F_temp;
                F.at(j) = F.at(j) - F_temp;
            } 
        }
    }

    return F;
}

Vector2d v_center(const vector<Vector2d> &v, int N)
{
    Vector2d vs = Vector2d::Zero();
    for(uint i = 0; i < v.size(); ++i)
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
    double Nf = 2*N-2;
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
        cout << r_vec << "\n\n";
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

vector<Vector2d> Beschleunigung(Vector2d (*f)(const Vector2d &, double), double (*pot)(const Vector2d &, double), const vector<Vector2d> &r_n, int N, double L, double &Epot)
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

                        Vector2d temp=r_n.at(i)-(r_n.at(j)+shift);
						a.at(i)+=f(temp, rc);

                        if(i < j){
                            Epot += pot(temp, rc);
                        }
                    }
                }
            }
        }
    }
    return a;
}

vector<Vector2d> verlet_r(Vector2d (*f)(const Vector2d &, double), double (*pot)(const Vector2d &, double), const vector<Vector2d> &r_n, vector<Vector2d> r_nminus1, double h, int N, double L, double &Epot)
{
    vector<Vector2d> r_nplus1, a;
    Vector2d r_nplus1_vec;
    a = update_kraft(r_n, N, L);
    for (int i = 1; i<N; i++)
    {
        cout << a.at(i) << "\n\n";
    }

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


void Aequilibrierung(Vector2d (*f)(const Vector2d &, double), double (*pot)(const Vector2d &, double), int N, double L, double particlesPerRow, double T0, uint steps, double h )
{
    ofstream file;
    file.open("build/aequi.txt");
    file << "t  Ekin    Epot    T   vs_x    vs_y \n\n";

    auto [r_n, v_n] = init(L, N, particlesPerRow, T0);
    vector<Vector2d> r_minus1, a_n, r_nplus1;
    Vector2d r_temp;
    double Epot = 0;
    a_n = update_kraft(r_n, N, L);

    file << 0 << "\t" << ekin(v_n) << "\t" << Epot << "\t" << T(v_n, N) << "\t" << 0 << "\t" << 0 << "\n";

    for(int i = 0; i < N; ++i){
        r_temp = r_n.at(i)-v_n.at(i)*h + 0.5 * a_n.at(i)*h*h;
        r_minus1.push_back(r_temp);
        //cout << r_temp << "\n\n";
    }

    for (uint i = 1; i < steps; i++)
    {
        cout << "i = " << i << "\n\n"; 

        r_nplus1 = verlet_r(f, pot, r_n, r_minus1, h, N, L, Epot);
        //for (int i = 1; i<N; i++)
        //{
        //    cout << r_nplus1.at(i) << "\n\n";
//
        //}
        periodische_RB(r_nplus1, L);
        v_n = verlet_v(r_minus1, r_nplus1, h, N);

        Vector2d vs = v_center(v_n, N);
        double vs_x, vs_y;
        vs_x = vs(0);
        vs_y = vs(1);


        file << i*steps << "\t" << ekin(v_n) << "\t" << Epot << "\t" << T(v_n, N) << "\t" << vs_x << "\t" << vs_y << "\t" << "\n";

    }

    file.close();

}

 


int main()
{
    int N = 16;
    double L, particlesPerRow, T, h;
    L = 8;
    particlesPerRow = 4;
    T = 1;
    h = 0.01;
    uint steps = 10;


    Aequilibrierung(KraftLJ, PotentialLJ, N, L, particlesPerRow, T, steps, h);

    return 0;
}