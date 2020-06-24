#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <stdlib.h>

using namespace std;
using namespace Eigen;

Vector2d verlet(function<Vector2d(Vector2d)> F, Vector4d y, Vector4d y_m_1, double &c)
{
  double m = 1.;
  double h = 0.01;
  Vector4d ret;
  Vector2d r_n_m_1 << y_m_1(0),y_m_1(1);
  Vector2d r_n << y(0), y(1);
  Vector2d v_n_m_1 << y_m_1(2),y_m_1(3);
  Vector2d v_n << y(2), y(3);
  Vector2d a_n = F(r_n)/m;
  Vector2d r_n_p_1 = 2 * r_n - r_n_m_1 + a_n * pow(h,2.);
  Vector2d v_n_p_1 = r_n_p_1 - r_n_m_1;
  v_n_p_1 /= 2*h;
  ret << r_n_p_1(0), r_n_p_1(1), v_n_p_1(0), v_n_p_1(1);
  ++c;
  return ret;
}

// =================================================================================================
//                      PROGRAMMSTRUKTUR
//
//                          ===========
//                          | Dataset |
//                          ===========
//                              |
//                           ========
//                           | Data |
//                           ========
//                              |
// ==============             ======                =============
// | Thermostat | ----------- | MD | -------------- | Potential |
// ==============             ======                =============
//                              |
//                           ========
//                           | main |
//                           ========
//
// - Die Klasse MD beinhaltet die primäre Logik des Algorithmus
// - Die Klasse Thermostat und Potential sind davon getrennt, damit unterschiedliche
//   Thermostate und Potentiale durch Vererbung implementiert und flexibel verwendet
//   werden können
// - Die Klasse Data speichert die in der MD-Simulation gespeicherten Daten und
//   kümmert sich um das Abspeichern
//      - Dataset ist eine Datensatz aus Zeit, Temperatur, ...
//      - Data hält mehrere Datensätze und noch einige Daten, die nicht zeitaufgelöst
//        gespeichert werden
//      - Statt umständlich getter und setter zu verwenden, sind die Member von Data
//        und Dataset public, da es sich um simple Datencontainer handelt
// - main() ruft MD mit den für die Aufgabenteile notwendigen Parametern auf
//
// Hinweise zu den verwendeten Vectoren:
//      - Wegen der Performance verwenden wir Vector2d statt VectorXd
//        Das macht jedoch den möglichen Übergang zu 3d-Zuständen umständlicher als mit VectorXd
//      - Für Listen von Daten wird std::vector werden
// =================================================================================================


// ================================ Potential-Klasse ===============================================

// Virtuelle Klasse, aus der konkrete Potentiale vererbt werden können
// (Hier nur Lennard-Jones nötig, aber so kann man schnell andere Potentiale implementieren)
class Potential
{
    public:
        virtual double      V               ( double r2 ) const = 0;  // Virtuelle Funktion
        virtual Vector2d    F               ( Vector2d r ) const = 0; // Virtuelle Funktion
};

class PotentialLJ: public Potential
{
    public:
        double              V               ( double r2 ) const;  // Überschreibt virtuelle Funktion
        Vector2d            F               ( Vector2d r ) const; // Überschreibt virtuelle Funktion
};

// Für Potential reicht das Quadrat der Vektorlänge, womit eine Wurzelberechnung gespart wird
double PotentialLJ::V ( double r2 ) const
{
    double r6 = 1./pow(r2,3.);
    if(r2 <= L/2.) return 4*(pow(r6,2.)-r6)
    else return 0.;
}

Vector2d PotentialLJ::F ( Vector2d r ) const
{
    double r2 = r.squaredNorm();
    double r6 = 1./pow(r2,3.);
    return 24*r*(r6/r2-pow(r6,2.)/r2);
}

// ------------------------------ Ende Potential-Klasse ------------------------------------------
// ================================ Thermostat-Klasse ===============================================

// Virtuelle Klasse, aus der konkrete Thermostate vererbt werden können
class Thermostat
{
    public:
        virtual void        rescale         ( vector<Vector2d>& v, double T ) const = 0;
};

// Kein Thermostat
class NoThermostat: public Thermostat
{
    public:
        void                rescale         ( vector<Vector2d>& v, double T ) const {} // Macht nichts
};

// Isokinetisches Thermostat für Aufgabe d)
class IsokinThermostat: public Thermostat
{
    public:
        void                rescale         ( vector<Vector2d>& v, double T ) const;
};

void IsokinThermostat::rescale( vector<Vector2d>& v, double T ) const
{
    double v_sum = 0.;
    for (Vector2d& n : v)
    {
        v_sum += n.squaredNorm();
    }
    double scale = 2*(N-1)*T/v_sum;
    scale = sqrt(scale);
    for (uint i = 0; i < v.size(); i++)
    {
        v.at() *= scale;
    }
}

// ------------------------------ Ende Thermostat-Klasse ------------------------------------------
// ================================ Data-Structs ===============================================

// Datensatz für zeitaufgelöste Daten
// (structs im Grunde gleich zu class, aber alle Member sind standardmäßig public)
struct Dataset
{
    double                  t, T, Ekin, Epot;
    Vector2d                vS;
};

// Rückgabedaten der MD-Simulation
// Konstruktor Data data(n); reserviert Speicher und füllt Paarkorrelationsfunktion mit 0en
struct Data
{
    vector<Dataset>         datasets; // Zeitaufgelöste Datensätze
    vector<double>          rBin, g;  // Gemittelte Paarkorrelationsfunktion
    vector<Vector2d>        r;        // Momentaufnahme der finalen Position
                                      //    Für Aufgabe e) kann es sinnvoll sein, r stattdessen
                                      //    in den zeitaufgelösten Datasets abzuspeichern

                            Data            ( uint n, uint numBins, double binSize );
    void                    save            ( const string& filenameSets,
                                              const string& filenameG,
                                              const string& filenameR ) const;
};

Data::Data( uint n, uint numBins, double binSize ):
    datasets( n ),  // Initializer list, weil sie Konstruktoren der Member auftruft
    rBin( numBins ),
    g( numBins, 0. ),
    r( 0 )
{
  /*TODO*/
}

void Data::save ( const string& filenameSets, const string& filenameG, const string& filenameR ) const
{
    ofstream out_set(filenameSets, ofstream::trunc);
    ofstream out_g(filenameG, ofstream::trunc);
    ofstream out_r(filenameR, ofstream::trunc);
    out_set << "#t, T, Ekin, Epot, vS \n\n";
    out_g << "#g \n\n";
    out_r << "#r \n\n";
    MatrixXd M_set(6,datasets.size());
    Dataset set;
    for(int i = 0; i < datasets.size(); ++i)
    {
      set = datasets.at(i);
      M_set.col(i) << set.t, set.T, set.Ekin, set.Epot, set.vS(0), set.vS(1);
    }
    out_set << M << endl;
    out_set.close();

    for(int i = 0; i< g.size(); ++i)
    {
      out_g << g.at(i) << "\n";
    }
    out_g.flush();
    out_g.close();
    
    MatrixXd M_r(2,r.size());
    for(int i = 0; i < r.size(); ++i)
    {
      M_r.col(i) = r.at(i);
    }
    out_r << M_r << endl;
    out_r.close();
    M_set.resize(0,0);
    M_r.resize(0,0);
}

// ------------------------------ Ende Data-Structs ------------------------------------------
// ================================ MD-Klasse ===============================================

class MD
{
    public:
                            MD              ( double L, uint N, uint particlesPerRow, double T,
                                            Potential& potential, Thermostat& thermostat,
                                            uint numBins = 1000 );

        void                equilibrate     ( const double dt, const unsigned int n );
        Data                measure         ( const double dt, const unsigned int n );

    private:
        vector<Vector2d>    r, v;
        double              L;
        uint                N;
        Potential&          potential;
        Thermostat&         thermostat;
        double              t = 0.;

        uint                numBins;
        double              binSize;

        // Teilchen werden in Box [0,L]x[0,L] verschoben
        void                centerParticles ();

        // Berechnungen wichtiger Messgrößen
        double              calcT           () const;
        double              calcEkin        () const;
        double              calcEpot        () const;
        Vector2d            calcvS          () const;
        Dataset             calcDataset     () const;

        // Berechnung der Beschleunigung
        //  Um redundante Rechnungen zu vermeiden, kann es sinnvoll sein, das Histogram
        //  bei der Berechnung der Beschleunigungen zu aktualisieren, daher wird es
        //  hier als Referenz übergeben
        vector<Vector2d>    calcAcc         ( vector<double>& hist ) const;

        // Berechnung des Abstandsvektors zwischen Teilchen r[i] und nähesten Spiegelteilchen von r[j]
        Vector2d            calcDistanceVec ( uint i, uint j ) const;
};

// Initialisierung des Systems per Konstruktor
MD::MD( double L, uint N, uint particlesPerRow, double T,
        Potential& potential, Thermostat& thermostat,
        uint numBins ):
    L(L),
    N(N),
    potential( potential ),
    thermostat( thermostat ),
    numBins( numBins ),
    binSize( L/(2*numBins) )
{
    Vector2d r_vec;
    for(int n = 0; n < particlesPerRow - 1; ++n)
    {
      for(int m = 0; m < particlesPerRow - 1; ++n)
      {
        r_vec << 1 + 2 * n, 1 + 2 * m;
        r.push_back(r_vec);
      }
    }
    Vector2d v_vec;
    for(int i = 0; i < N; ++i)
    {
      v_vec << double(rand()%10), double(rand()%10);
      v.push_back(v_vec);
    }
    Vector2d v_s = calcvS();
    for (Vector2d& n : v)
    {
        n -= v_s;
    }
    r_vec.resize(0);
    v_vec.resize(0);
    v_s.resize(0);
}

// Integration ohne Datenaufnahme zur reinen Äquilibrierung
void MD::equilibrate ( const double dt, const unsigned int n )
{
    /*TODO*/
}

Data MD::measure ( const double dt, const unsigned int n )
{
    /*TODO*/
}

void MD::centerParticles()
{
    for(uint i = 0; i < r.size(); ++i)
    {
      for(uint j = 0; j < r.at(i); ++j)
      {
        if(r.at(i)(j) > L) r.at(i)(j) = fmod(r.at(i)(j), L);
        if(r.at(i)(j) < 0) r.at(i)(j) = fmod(r.at(i)(j), L) +L;
      }
    }
}

double MD::calcT() const
{
//    return calcEkin()*(2./((2*N-2)*1.3806503*pow(10.,-23.)));// mit k_B
    return calcEkin()*(2./((2*N-2));
}

double MD::calcEkin() const
{
    double m = 1.;
    double ekin = 0.;
    for(uint i = 0; i < v.size(); ++i)
    {
      ekin += m*v.at(i).squaredNorm();
    }
    return ekin/2.;
}

double MD::calcEpot() const
{
    double epot = 0.;
    double r2;
    Vector2d dr;
    for(uint i = 0; i < r.size()-1; ++i)
    {
        for(uint j = i+1; j < r.size(); ++i)
        dr = calcDistanceVec(i,j);
        r2 = dr.squaredNorm();
        epot += potential(r2);
    }
    return epot;
}

Vector2d MD::calcvS() const
{
    Vector2d vs = Vector2d::Zero();
    for(int i = 0; i < v.size(); ++i)
    {
      vs += v.at(i);
    }
    return vs;
}

Dataset MD::calcDataset() const
{
    Dataset set;
    set.t = t;
    set.T = calcT();
    set.Ekin = calcEkin();
    set.Epot = calcEpot();
    set.vS = calcvS();
    return set;
}

Vector2d MD::calcDistanceVec( uint i, uint j ) const
{
    double r_cut = L/2.;
    MatrixXd M(2,4);
    M << 0, L, L, L,
         L, 0, L, -L;
    Vector2d vec_shift;
    Vector2d vec_dr;
    bool br = false;
    for(uint i = 0; i < M.cols(); ++i)
    {
      vec_shift = M.col(i);
      vec_dr = r.at(i)-(r.at(j)+vec_shift)
      if(vec_dr < r_cut)
      {
        br = true;
        break;
      }
      vec_dr = r.at(i)-(r.at(j)-vec_shift)
      if(vec_dr < r_cut)
      {
        br = true;
        break;
      }
    }
    if(br) return vec_dr;
    else return r.at(i) - r.at(j);
}

vector<Vector2d> MD::calcAcc( vector<double>& hist ) const
{
    /*TODO*/
}

// ------------------------------ Ende MD-Klasse ------------------------------------------


int main(void)
{
    PotentialLJ      LJ;
    NoThermostat     noThermo;
    IsokinThermostat isoThermo;

    const uint partPerRow       = 4;
    const uint N                = 16;
    const double L              = 8;
    const int numBins           = /*TODO*/;

    // b) Äquilibrierungstest
    {
        const double T          = 1.;
        const double dt         = /*TODO*/;
        const uint steps        = /*TODO*/;

        MD md( L, N, partPerRow, T, LJ, noThermo, numBins );
        md.measure( dt, steps ).save( "b)set.dat", "b)g.dat", "b)r.dat" );
    }

    // c) Paarkorrelationsfunktion
    string TstringVec[3] = { "0.01", "1", "100" };
    for ( auto& Tstring: TstringVec )
    {
        const double T          = stod(Tstring);
        const double dt         = /*TODO*/;
        const uint equiSteps    = /*TODO*/;
        const uint steps        = /*TODO*/;

        MD md( L, N, partPerRow, T, LJ, noThermo, numBins );
        md.equilibrate( dt, equiSteps );
        md.measure( dt, steps ).save( "c)set" + Tstring + ".dat", "c)g" + Tstring + ".dat", "c)r" + Tstring + ".dat" );
    }

    // d) Thermostat
    {
        /*TODO*/
    }

    // e) Animation
    {
        /*TODO*/
    }

    return 0;
}
