#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main(){

    // Teilaufgabe b)

    Matrix3f A, L, U, P;
    Vector3f x, b;

    b << 2,0,2;
    A << 0.5, -0.5, 0,          sqrt(3)/2, sqrt(3)/2, 0,            0, 0, 1;

    // Bestimmung von L, U und P mittels Eeigen

    x = A.partialPivLu().solve(b);
    L = A.partialPivLu().matrixLU().triangularView<UpLoType::UnitLower>();
    U = A.partialPivLu().matrixLU().triangularView<UpLoType::Upper>();
    P = A.partialPivLu().permutationP();

    cout << "x =\n" << x << "\n\n";
    cout << "L =\n" << L << "\n\n";
    cout << "U =\n" << U << "\n\n";
    cout << "P =\n" << P << "\n\n";

    // Teilaufgabe c)

        // A ändert sich nicht, also können L, U und P wiederverwendet werden und es muss nur y' bestimmt werden.

    return 0;
}