#include<iostream>
#include<random>
#include<time.h>
#include<fstream>
using namespace std;

int spin[100][100];		//Spin-Array

double summe() {
	double summe = 0;
	for (int i = 0; i < 100; i++) {
		for (int j = 0; j < 100; j++) {
			summe += spin[i][j];
		}
	}
	return summe;
}

void print(string s) {
	ofstream stream;
	stream.open(s + ".txt");
	for (int i = 0; i < 100; i++) {		//Array speichern
		for (int j = 0; j < 100; j++) {
			stream << spin[i][j] << "\t";
		}
		stream << endl;
	}
	stream.close();

}

void iniFest() {
	for (int i = 0; i < 100; i++) {
		for (int j = 0; j < 100; j++) {
			spin[i][j] = -1;
		}
	}

}

void iniZufall(mt19937 rng, uniform_real_distribution<> dist) {
	for (int i = 0; i < 100; i++) {		//Array initialisieren mit -1=down
		for (int j = 0; j < 100; j++) {
			double zufallsSpin = dist(rng);

			if (zufallsSpin < 0.5) {
				spin[i][j] = -1;
			} else {
				spin[i][j] = 1;
			}
		}
	}

}

double energie() {
	double e = 0;
	int spinOben = 0, spinUnten = 0, spinRechts = 0, spinLinks = 0, spinMitte =
			0;
	for (int i = 0; i < 100; i++) {			//äußere Summe über I
		for (int j = 0; j < 100; j++) {

			spinOben = spin[i - 1][j];		//nächste Nachbarn
			spinUnten = spin[i + 1][j];
			spinLinks = spin[i][j - 1];
			spinRechts = spin[i][j + 1];

			switch (i) {
			case 0:
				spinOben = spin[99][j];
				break;
			case 99:
				spinUnten = spin[0][j];
				break;
			}
			switch (j) {
			case 0:
				spinLinks = spin[i][99];
				break;
			case 99:
				spinRechts = spin[i][0];
				break;
			}

			spinMitte = spin[i][j];
			e += spinMitte * spinOben + spinMitte * spinUnten
					+ spinMitte * spinLinks + spinMitte * spinRechts; //Summe über nächste Nachbarn
		}
	}
	return -e;
}

double deltaEnergie(int i, int j, int spin_old) {
	double deltaE = 0;
	int spinOben = spin[i - 1][j], spinUnten = spin[i + 1][j], spinRechts =
			spin[i][j + 1], spinLinks = spin[i][j - 1], spinMitte = spin[i][j];
	switch (i) {
	case 0:
		spinOben = spin[99][j];
		break;
	case 99:
		spinUnten = spin[0][j];
		break;
	}
	switch (j) {
	case 0:
		spinLinks = spin[i][99];
		break;
	case 99:
		spinRechts = spin[i][0];
		break;
	}

	deltaE = (spin_old - spinMitte)
			* (spinOben + spinUnten + spinLinks + spinRechts);

	return deltaE;
}

void warmUp(mt19937 rng, uniform_real_distribution<> dist, double beta) {
	for (int n = 0; n < pow(10,4); n++) {
		int zX = dist(rng) * 100;
		int zY = dist(rng) * 100;

		//zufälligen Spin flippen

		int spin_old = spin[zX][zY];
		spin[zX][zY] *= -1;   //Spin Flip

		double deltaE = deltaEnergie(zX, zY, spin_old);

		if (deltaE < 0 || dist(rng) < exp(-beta * deltaE)) {

		} else {
			spin[zX][zY] = spin_old;
		}

	}
}

int main() {
	mt19937 rng;				//Zufallszahlengenerator
	rng.seed(time(NULL));
	uniform_real_distribution<> dist(0, 1);
	iniZufall(rng, dist);
	print("ini");

	int zX = 0, zY = 0; 		//Zufallszahlenspeicher für x,y

	int spin_old = 0;

	double deltaE = 0;
	double beta = 1;

	//Aufgabentteil a)

	for (int b = 1; b <= 3; b += 2) {
		beta = 1 / b;

		//Sweeps
		for (int s = 0; s < 100000; s++) {
			for (int n = 0; n < pow(10, 4); n++) {
				zX = dist(rng) * 100;
				zY = dist(rng) * 100;

				//zufälligen Spin flippen

				spin_old = spin[zX][zY];
				spin[zX][zY] *= -1;   //Spin Flip

				deltaE = deltaEnergie(zX, zY, spin_old);

				if (deltaE < 0 || dist(rng) < exp(-beta * deltaE)) {

				} else {
					spin[zX][zY] = spin_old;
				}

			}

		}
		if (b == 1) {
			print("1kbt");
			cout << "1kbt" << endl;
		} else {
			print("3kbt");
			cout << "3kbt" << endl;
		}

		iniZufall(rng, dist);
	}
	cout << "a) fertig" << endl;

	//Aufgabenteil b)
	ofstream stream;
	for (int s = 0; s < 2; s++) {
		for (int b = 1; b <= 3; b++) {
			if (s == 0) {
				iniZufall(rng, dist);
			} else {
				iniFest();
			}

			beta = 1 / b;
			double E = energie();
			double eSum=0;


			if (s == 0) {
				switch (b) {
				case 1:
					stream.open("1kbt-b-zufall.txt");
					break;
				case 2:
					stream.open("2kbt-b-zufall.txt");
					break;
				case 3:
					stream.open("3kbt-b-zufall.txt");
				}
			} else {
				switch (b) {
				case 1:
					stream.open("1kbt-b-fest.txt");
					break;
				case 2:
					stream.open("2kbt-b-fest.txt");
					break;
				case 3:
					stream.open("3kbt-b-fest.txt");
				}
			}

			for (int n = 0; n < 4*pow(10, 5); n++) {
				zX = dist(rng) * 100;
				zY = dist(rng) * 100;

				//zufälligen Spin flippen

				spin_old = spin[zX][zY];
				spin[zX][zY] *= -1;   //Spin Flip

				deltaE = deltaEnergie(zX, zY, spin_old);

				if (deltaE < 0 || dist(rng) < exp(-beta * deltaE)) {
					E += deltaE;
				} else {
					spin[zX][zY] = spin_old;
				}
				eSum+=E;
				stream << eSum / (10000 * (n + 1)) << "\t" << n + 1 << endl;

			}
			stream.close();

		}
	}
	cout << "b) fertig" << endl;

	//Aufgabenteil c)
	double temperatur[] = { 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.21, 2.22, 2.23,
			2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.31, 2.32, 2.33, 2.34,
			2.35, 2.35, 2.36, 2.36, 2.37, 2.38, 2.39, 2.4, 2.6, 2.8, 3 };
	ofstream streamMT, streamcT;
	streamMT.open("m(T).txt");
	streamcT.open("c(T).txt");
	for (unsigned int b = 0; b < 32; b++) {
		iniZufall(rng, dist);
		//Aufwärmphase
		double beta = 1 / temperatur[b];
		warmUp(rng, dist, beta);
		double E = 0;
		double eSum=0;
		double sum=0;
		ofstream streamMt, streamE;
		if (temperatur[b] == 1) {
			streamMt.open("m(t)-1kbt.txt");
			streamE.open("e(t)-1kbt.txt");
			E = energie();
		} else if (temperatur[b] == 2) {
			streamMt.open("m(t)-2kbt.txt");
			streamE.open("e(t)-2kbt.txt");
			E = energie();
		} else if (temperatur[b] == 3) {
			streamMt.open("m(t)-3kbt.txt");
			streamE.open("e(t)-3kbt.txt");
			E = energie();
		}
		//Sweeps
		for (int s = 0; s < 100000; s++) {
			for (int n = 0; n < pow(10, 4); n++) {
				zX = dist(rng) * 100;
				zY = dist(rng) * 100;

				//zufälligen Spin flippen

				spin_old = spin[zX][zY];
				spin[zX][zY] *= -1;   //Spin Flip

				deltaE = deltaEnergie(zX, zY, spin_old);

				if (deltaE < 0 || dist(rng) < exp(-beta * deltaE)) {
					E += deltaE;
				} else {
					spin[zX][zY] = spin_old;
				}

			}
			eSum+=E;
			sum+=summe();
			if (temperatur[b] == 1 || temperatur[b] == 2
					|| temperatur[b] == 3) { //Magnetisierung/zeit und Energie/Zeit

				streamMt << sum / (10000 * (s + 1)) << "\t"
						<< abs(sum / (10000 * (s + 1))) << "\t" << s + 1
						<< endl;

				streamE << eSum / (10000 * (s + 1)) << "\t" << s + 1 << endl;
			}
		}
		streamMt.close();
		streamE.close();
		streamMT << sum / (10000 * 100000) << "\t"
				<< abs(sum / (10000 * 100000)) << "\t" << temperatur[b]
				<< endl;
		streamcT << pow(eSum, 2) / (100000) << "\t" << eSum / (100000) << "\t"
				<< temperatur[b] << endl;		//Aufgabenteil d)
		cout << temperatur[b] << endl;
	}
	streamMT.close();
	streamcT.close();

}

