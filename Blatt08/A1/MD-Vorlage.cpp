#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

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
    double L = 8.;
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
    /*TODO*/
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
    /*TODO*/
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
    binSize( /*TODO*/ )
{
    /*TODO*/
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
    /*TODO*/
}

double MD::calcT() const
{
    /*TODO*/
}

double MD::calcEkin() const
{
    /*TODO*/
}

double MD::calcEpot() const
{
    /*TODO*/
}

Vector2d MD::calcvS() const
{
    /*TODO*/
}

Dataset MD::calcDataset() const
{
    /*TODO*/
}

Vector2d MD::calcDistanceVec( uint i, uint j ) const
{
    /*TODO*/
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

    const uint partPerRow       = /*TODO*/;
    const uint N                = /*TODO*/;
    const double L              = /*TODO*/;
    const int numBins           = /*TODO*/;

    // b) Äquilibrierungstest
    {
        const double T          = /*TODO*/;
        const double dt         = /*TODO*/;
        const uint steps        = /*TODO*/;

        MD md( L, N, partPerRow, T, LJ, noThermo, numBins );
        md.measure( dt, steps ).save( "b)set.dat", "b)g.dat", "b)r.dat" );
    }

    // c) Paarkorrelationsfunktion
    string TstringVec[3] = { "0.01", "1", "100" };
    for ( auto& Tstring: TstringVec )
    {
        const double T          = stoi(Tstring);
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
