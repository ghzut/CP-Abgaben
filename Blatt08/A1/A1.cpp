#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <math.h>  // sqrt()
#include <random>
#include <string>

 
using namespace std;
using namespace Eigen;


/*  Berechnet die Kraft des Lennard-Jones-Potentials
 *  Keine Kraft bei Abstaenden groesser des Cutoff
 *  INPUT       r       2D-Ortsvektor zwischen zwei Teilchen
 *              rc      Cutoff
 *  OUTPUT      ljkraft Kraftvektor nach dem LJ-Potential
 */
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

/*  Findet Naechsten Nachbar von Teilchen r2 zu Teilchen r1
 *  Bei einem Cutoff von L/2 wechselwirkt Teilchen 1 bei r1 nur mit dem nächstgelegenen
 *  Teilchen 2, da alle anderen Bildteilchen von j weiter als L/2 enfternt liegen.
 *  INPUT       r1      Ortsvektor
 *              r2      Ortsvektor
 *              L       Kantenlaenge des Simulationsvolumens
 */
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


/*  Berechne wirkende Kraft auf jedes Teilchen
 *  INPUT       r       Ortsvektoren der Teilchen
 *              N       Anzahl der Teilchen
 *              F       (alte) Kraftvektoren der Teilchen
 *              L       Kantenlaenge des Gebiets [0,L]x[0,L]
 *  OUTPUT      F       wird ueberschrieben
 */
void update_kraft(MatrixXd &r, MatrixXd &F, unsigned int N, double L)
{
    double rc = L/2;
    F.setZero();
    // Folgende Groessen als Zwischenspeicher benoetigt
    Vector2d delta_r = Vector2d::Zero();     // Abstandsvektor
    Vector2d r_next = Vector2d::Zero();     // nächster Nachbar
    Vector2d F_temp = Vector2d::Zero();     // Kraft zwischen j und i

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


/*  Periodische Randbedingungen: Verhindere, dass Teilchen aus dem Simulationsgebiet
 *  heraus laufen
 *  EIN SIMULATIONSGEBIET, SIE ZU KNECHTEN, SIE ALLE ZU FINDEN,
 *  INS DUNKEL ZU TREIBEN UND EWIG ZU BINDEN.
 *  INPUT       r       2xN Matrix mit Ortsvektoren der Teilchen
 *              N       Anzahl der Teilchen
 *              L       Kantenlaenge des Gebiets [0,L]x[0,L]
 */
void periodische_RB(MatrixXd &r, unsigned int N, double L)
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


/*  Initialisiere Molekulardynamik Simulation
 *  INPUT       N       Anzahl Teilchen
 *              L       Kantenlaenge der Box
 *              r       2xN-Matrix fuer Orte
 *              v       2xN-Matrix fuer Geschwindigkeiten
 *              Tinit   Starttemperatur
 *  OUTPUT      r und v werden ueberschrieben
 */
void init(unsigned int N, double L, MatrixXd &r, MatrixXd &v, double Tinit)
{
    //double sqrtN = sqrt(N);
    srand(42);

    for (int n=0; n<4; n++)
    {
        for (int m=0; m<4; m++)
        {
            // Ordne Punkte auf einem Gitter an
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

    // Geschwindigkeiten werden mit @skal umskaliert, um Temperatur auf Tinit zu setzen.
    double Nf = 2*N-2;                              // Anzahl Freiheitsgrade
    double skal = Tinit*Nf/v.colwise().squaredNorm().sum();
    for (int n=0; n<N; n++)
    {
        v.col(n) = sqrt(skal)*v.col(n);
    }
}


/*  Berechne potentielle Energie
 *  INPUT       r       2xN Matrix der Ortsvektoren
 *              N       Anzahl der Teilchen
 *              L       Kantenlaenge des Simulationsvolumens
 *  OUTPUT      Epot    Potentielle Energie der Konfiguration
 */
double pot(MatrixXd &r, unsigned int N, double L)
{
    double rc = L/2;
    // Folgende Groessen als Zwischenspeicher benoetigt
    VectorXd delta_r = VectorXd::Zero(2, 1);     // Abstandsvektor
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
                // Lennard-Jones-Potential
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
    file << "t vsx vsy Ekin Epot T \n\n";

    int Nf = 2*N-2;            // Anzahl Freiheitsgrade (in 2D)
    MatrixXd r = MatrixXd::Zero(2, N);
    MatrixXd v = MatrixXd::Zero(2, N);

    init(N, L, r, v, T0);

    double t, Ekin, Epot, T; 
    Vector2d vS;
    t = 0;

    MatrixXd F_alt, F;
    //r_alt = MatrixXd::Zero(2,N);

    F = MatrixXd::Zero(2,N);
    update_kraft(r, F, N, L);

    //for (int i=0; i<N; i++)
    //{
    //    r_alt.col(i) = r.col(i) - v.col(i) * h + 0.5 * h * h * F.col(i);
    //}

    while (t < t_aequi)
    {
        vS = 1./N * v.rowwise().sum();
        Ekin = 0.5*v.colwise().squaredNorm().sum();
        Epot = pot(r, N, L);
        T = 2*Ekin/Nf;

        file << t << " " << vS(0) << " " << vS(1) << " " << Ekin << " " << Epot << " " << T << "\n";

        F_alt = F;
        // Neue Ortsvektoren
        for (int i=0; i<N; i++)
        {
            r.col(i) = r.col(i)+v.col(i)*h+0.5*h*h*F.col(i);
        }
        // periodische RB: verhindere Herauslaufen aus Box
        periodische_RB(r, N, L);
        // Neue Kraefte
        update_kraft(r, F, N, L);
        // neue Geschwindigkeiten
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

    struct vecs{
        MatrixXd Ort, Geschwindigkeit, Kraft;
    };

    return vecs{r, v, F};

}

VectorXd Paarkorrelation(const MatrixXd &r, const MatrixXd &v, double L, int N, int N_h) {
	//int row= Yn.rows();
	//int col= Yn.cols();
	//int N=g.size();
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
						//cout << r/dr << "\t" << l << "\n" << "\n";
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
    //MatrixXd r = MatrixXd::Zero(2, N);
    //MatrixXd v = MatrixXd::Zero(2, N);

    auto [r, v, F] = aequilibrierung(N, L, T0, t_aequi, h, thermo, filename);

    int Nf = 2*N-2; 

    MatrixXd F_alt;

    setprecision(2);
    if (thermo)
    {
        file.open("build/messung_isokinetisch"+filename+".txt", ios::trunc);
    }
    else
    {
        file.open("build/messung"+filename+".txt", ios::trunc);
    }

    file.precision(10);
    
    //file.precision(10);
    file << "t  T \n";
    while (t < t_max+t_aequi)
    {
        // Messe Observablen des Systems
        file << t-t_aequi;
        // berechne Temperatur
        double T = v.colwise().squaredNorm().sum()/Nf;
        file << " " << T << "\n";
        // TODO: Paarkorrelationsfunktion implementieren und hier einfuegen

        F_alt = F;
        // Neue Ortsvektoren
        for (int i=0; i<N; i++)
        {
            r.col(i) = r.col(i)+v.col(i)*h+0.5*h*h*F.col(i);
        }
        // periodische RB: verhindere Herauslaufen aus Box
        periodische_RB(r, N, L);
        // Neue Kraefte
        update_kraft(r, F, N, L);
        // neue Geschwindigkeiten
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

    // Einstellbare Parameter
    unsigned int N = 16;    // Anzahl Teilchen
    double L = 8;           // Kantenlaenge der Box
    double h = 0.01;        // Zeitschritt des Verlet-Algorithmus
    double T0 = 1;       // Temperatur bei Initialisierung t=0
    double t_aequi = 50;     // Äquilibrierungszeit
    double t_max = 500;      // Maximale Zeit der Simulation


    simulation(N, L, T0, t_aequi, t_max, h, false, "1");
    simulation(N, L, T0, t_aequi, t_max, h, true, "1");

    T0 = 0.01;

    simulation(N, L, T0, t_aequi, t_max, h, false, "0.01");
    simulation(N, L, T0, t_aequi, t_max, h, true, "0.01");

    T0 = 100;

    simulation(N, L, T0, t_aequi, t_max, h, false, "100");
    simulation(N, L, T0, t_aequi, t_max, h, true, "100");


    // Aufgabenteil b) Betrachte Äquilibrierungsphase
    //cout << "MD-Simulation 1 - T=1" << endl;
    //md_simulation(r, v, N, L, h, Tinit,
    //        tequi, "build/equilibration.txt", true,
    //        tmax, "build/messung_T1.txt",
    //        thermostat);
    //cout << endl;

    // Aufgabenteil c) Betrachte T und Paarkorrelationsfunktion
    //Tinit = 100;
    //// passe bei so hohen Temperaturen die Schrittweite an
    //h = 0.001;
    //tequi = 5;
    //tmax = 50;
    //cout << "MD-Simulation 2 - T=100" << endl;
    //md_simulation(r, v, N, L, h, Tinit,
    //        tequi, "build/dummy.txt", false,
    //        tmax, "build/messung_T1e2.txt",
    //        thermostat);
    //// cout << r << endl << v << endl;
    //cout << endl;
    //h = 0.01;
    //tequi = 50;
    //tmax = 500;
//
    //Tinit = 0.01;
    //cout << "MD-Simulation 3 - T=0.01" << endl;
    //md_simulation(r, v, N, L, h, Tinit,
    //        tequi, "build/dummy.txt", false,
    //        tmax, "build/messung_T1e-2.txt",
    //        thermostat);
    //cout << endl;
//
//
    //// Aufgabenteil d) Thermostat
    //thermostat = true;
    //cout << "MD-Simulation 4 - Thermostat eingebaut" << endl;
    //md_simulation(r, v, N, L, h, Tinit,
    //        tequi, "build/thermostat_equi.txt", true,
    //        tmax, "build/thermostat_mess.txt",
    //        thermostat);
    //cout << endl;
    //// Konfiguration am letzten Zeitschritt
    //file.open("build/endkonfig.txt", ios::trunc);
    //save_data(r, v, file);
    //file.close();
//
    //cout << "Ende des Programms!" << endl;
    return 0;
}
