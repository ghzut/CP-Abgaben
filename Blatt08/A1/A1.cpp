#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <math.h>
#include <random>
#include <string>


using namespace std;
using namespace Eigen;

VectorXd KraftLJ(VectorXd r, double rc)
{
    VectorXd f = VectorXd::Zero(2,1);
    if (r.norm()>rc)
    {
        f << 0, 0;
    }
    else
    {
        f = r*24*(2*pow(r.norm(),-14)-pow(r.norm(),-8));
    }
    return f;
}

double PotentialLJ(Vector2d delta_r)
{
    return 4*(pow(delta_r.norm(), -12)-pow(delta_r.norm(), -6));
}

VectorXd naechsterNachbar(VectorXd r1, VectorXd r2, double L)
{
    double rc = L/2;
    VectorXd shift = VectorXd::Zero(2, 1);
    VectorXd diff = VectorXd::Zero(2, 1);

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

    // Kein Bildteilchen nah genug => reale Teilchen
    return r2;
}


void update_kraft(MatrixXd &r, MatrixXd &F, int N, double L)
{
    double rc = L/2;
    F.setZero();

    Vector2d delta_r = Vector2d::Zero();
    Vector2d r_next = Vector2d::Zero();
    Vector2d F_temp = Vector2d::Zero();
    for (int i=0; i<(N-1); i++)
    {
        for (int j=i+1; j<N; j++)   // Keine Selbstwecheelwirkung oder Doppelzaehlung
        {
            r_next = naechsterNachbar(r.col(i), r.col(j), L);
            delta_r = r.col(i)-r_next;
            if(delta_r.norm()<rc)
            {
                F_temp = KraftLJ(delta_r, L);
                F.col(i) = F.col(i) + F_temp;
                F.col(j) = F.col(j) - F_temp;
            }
        }
    }
}

void periodische_RB(MatrixXd &r, int N, double L)
{


    /* Das funktioniert bei T=100 leider nicht mehr */

    //for(int i = 0; i < 2; ++i)
    //{
    //  for(int j = 0; j < N; ++j)
    //  {
    //    if(r(i,j) > L || r(i,j) < 0)
    //    {
    //        r(i,j) -= floor(r(i,j)/L) * L;
    //    }
    //  }
    //}

    for (int n=0; n<N; n++)
    {
        // Teilchen nach links abgedriftet
        if(r(0,n) < 0)
        {
            r(0,n) += L;
        }
        // Teilchen nach rechts abgedriftet
        if(r(0,n) > L)
        {
            r(0,n) -= L;
        }
        // Teilchen nach unten abgedriftet
        if(r(1,n) < 0)
        {
            r(1,n) += L;
        }
        // Teilchen nach oben abgedriftet
        if(r(1,n) > L)
        {
            r(1,n) -= L;
        }
    }
}


void init(int N, double L, MatrixXd &r, MatrixXd &v, double Tinit)
{
    srand(42);

    for (int n=0; n<4; n++)
    {
        for (int m=0; m<4; m++)
        {
            // Gitteranordnung
            r(0, 4*n+m) = 1 + 2 * n;
            r(1, 4*n+m) = 1 + 2 * m;
            // Zufaellige Geschwindigkeiten
            v(0, 4*n+m) = double(rand()%11);
            v(1, 4*n+m) = double(rand()%11);
        }
    }

    // Nun Schwerpunktsbewegung auf Null setzen
    Vector2d vmean = 1./N * v.rowwise().sum();
    for (int n=0; n<N; n++)
    {
        v.col(n) = (v.col(n) - vmean);
    }

    double Nf = 2*N-2;
    double skal = Tinit*Nf/v.colwise().squaredNorm().sum();
    for (int n=0; n<N; n++)
    {
        v.col(n) = sqrt(skal)*v.col(n);
    }
}


double pot(MatrixXd &r, int N, double L)
{
    double rc = L/2;
    VectorXd delta_r = VectorXd::Zero(2, 1);
    VectorXd r_next = VectorXd::Zero(2, 1);
    double Epot = 0.;

    for (int i=0; i<(N-1); i++)
    {
        for (int j=i+1; j<N; j++)   // Keine Selbstwecheelwirkung oder Doppelzaehlung
        {
            r_next = naechsterNachbar(r.col(i), r.col(j), L);
            delta_r = r.col(i)-r_next;
            if(delta_r.norm()<rc)
            {
                Epot += PotentialLJ(delta_r);
            }
        }
    }
    return Epot;
}

auto aequilibrierung(int N, double L, double T0, double t_aequi, double h, bool thermo, string filename)
{
    ofstream file;

    setprecision(2);
    if (thermo)
    {
        file.open("build/aequi_isokinetisch"+filename+".txt", ios::trunc);
    }
    else
    {
        file.open("build/aequi"+filename+".txt", ios::trunc);
    }

    file.precision(10);
    file << "#t vsx vsy Ekin Epot T \n\n";

    int Nf = 2*N-2;            // Anzahl Freiheitsgrade (in 2D)
    MatrixXd r = MatrixXd::Zero(2, N);
    MatrixXd v = MatrixXd::Zero(2, N);

    init(N, L, r, v, T0);

    double t, Ekin, Epot, T;
    Vector2d vS;
    t = 0;

    MatrixXd F_alt, F;

    F = MatrixXd::Zero(2,N);
    update_kraft(r, F, N, L);

    while (t < t_aequi)
    {
        vS = 1./N * v.rowwise().sum();
        Ekin = 0.5*v.colwise().squaredNorm().sum();
        Epot = pot(r, N, L);
        T = 2*Ekin/Nf;

        file << t << " " << vS(0) << " " << vS(1) << " " << Ekin << " " << Epot << " " << T << "\n";
        F_alt = F;

        for (int i=0; i<N; i++)
        {
            r.col(i) = r.col(i)+v.col(i)*h+0.5*h*h*F.col(i);
        }

        periodische_RB(r, N, L);
        update_kraft(r, F, N, L);

        for (int i=0; i<N; i++)
        {
            v.col(i) = v.col(i)+0.5*h*(F.col(i)+F_alt.col(i));
        }

        if (thermo)
        {
            double skal = T0*Nf/v.colwise().squaredNorm().sum();
            for (int i=0; i<N; i++)
            {
                v.col(i) = sqrt(skal)*v.col(i);
            }
        }

        t += h;
    }

    file.close();

    struct vecs{
        MatrixXd Ort, Geschwindigkeit, Kraft;
    };

    return vecs{r, v, F};

}

VectorXd Paarkorrelation(const MatrixXd &r, const MatrixXd &v, double L, int N, int N_h) {
    VectorXd g = VectorXd::Zero(N_h);
	double dr=L/(2*N_h);

	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < j; i++)
		{
			Vector2d shift=Vector2d::Zero();
			for (int i = -1; i <= 1; ++i)
			{
				shift(0)=i*L;
				for (int j = -1; j <= 1; ++j)
				{
					shift(1)=j*L;
					VectorXd rel=r.col(i)-(r.col(j)+shift);
					double r =rel.norm();
					if (r<L/2)
					{
						int l = (int) (r/dr);
						g(l)+=L*L/(N*N*M_PI*(-(l*dr)*(l*dr)+((l+1)*dr)*((l+1)*dr)));
					}
				}
			}
		}
	}
}

void simulation(int N, double L, double T0, double t_aequi, double t_max, double h, bool thermo, string filename)
{
    ofstream file;

    double t = t_aequi;

    auto [r, v, F] = aequilibrierung(N, L, T0, t_aequi, h, thermo, filename);
    int Nf = 2*N-2;
    MatrixXd F_alt;

    if (thermo)
    {
        file.open("build/messung_isokinetisch"+filename+".txt", ios::trunc);
    }
    else
    {
        file.open("build/messung"+filename+".txt", ios::trunc);
    }

    file.precision(10);
    file << "#t  T \n";
    
    while (t < t_max+t_aequi)
    {
        file << t-t_aequi;
        double T = v.colwise().squaredNorm().sum()/Nf;
        file << " " << T << "\n";
        // TODO: Paarkorrelationsfunktion

        F_alt = F;
        for (int i=0; i<N; i++)
        {
            r.col(i) = r.col(i)+v.col(i)*h+0.5*h*h*F.col(i);
        }
        periodische_RB(r, N, L);
        update_kraft(r, F, N, L);
        for (int i=0; i<N; i++)
        {
            v.col(i) = v.col(i)+0.5*h*(F.col(i)+F_alt.col(i));
        }

        // Isokinetisches Thermostat
        if (thermo)
        {
            double skal = T0*Nf/v.colwise().squaredNorm().sum();
            for (int i=0; i<N; i++)
            {
                v.col(i) = sqrt(skal)*v.col(i);
            }
        }
        t += h;
    }

    file.close();

}


int main()
{
    unsigned int N = 16;
    double L = 8;
    double h = 0.01;
    double t_aequi = 500;
    double t_max = 5000;

    double T0 = 1;
    double T1 = 0.01;
    double T2 = 100;

    omp_set_num_threads(6);

    #pragma omp parallel sections
    {
        #pragma omp section
            simulation(N, L, T0, t_aequi, t_max, h, false, "1");

        #pragma omp section
            simulation(N, L, T0, t_aequi, t_max, h, true, "1");

        #pragma omp section
            simulation(N, L, T1, t_aequi, t_max, h, false, "001");

        #pragma omp section
            simulation(N, L, T1, t_aequi, t_max, h, true, "001");
        
        #pragma omp section
            simulation(N, L, T2, t_aequi, t_max, h, false, "100");

        #pragma omp section
            simulation(N, L, T2, t_aequi, t_max, h, true, "100");
    }

    return 0;
}
