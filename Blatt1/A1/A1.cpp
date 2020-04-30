#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;




Vector3f solveLGS(Matrix3f L, Matrix3f U, Matrix3f P, Vector3f b)
{
    Vector3f x, y, z;
    double xx = 0;


    z = P.transpose() * b;


    y(0) = z(0)/L(0,0);

    for (int i=1; i<b.size(); i++){
        
        for (int j=0; j<i; j++){
            xx += L(i, j) * y(j);
        }

        y(i) = 1/L(i, i) * (z(i) - xx);
    }

    x(b.size()-1) = y(b.size()-1)/U(b.size()-1,b.size()-1);

    xx = 0;

    for (int i=b.size()-2; i>=0; i--){
        for (int j=i+1; j<b.size(); j++){
            xx += U(i, j) * x(j);
        }

        x(i) = 1/U(i, i) * (y(i) - xx);
    }

    return x;
}




int main()
{

    // Teilaufgabe b)
    cout << "___________b)____________\n\n";

    Matrix3f A, L, U, P;
    Vector3f b, x;

    b << 2,0,2;
    A << 0.5, -0.5, 0,          sqrt(3)/2, sqrt(3)/2, 0,            0, 0, 1;

    // Bestimmung von L, U und P mittels Eigen

    // x = A.partialPivLu().solve(b);
    L = A.partialPivLu().matrixLU().triangularView<UpLoType::UnitLower>();
    U = A.partialPivLu().matrixLU().triangularView<UpLoType::Upper>();
    P = A.partialPivLu().permutationP();

    x = solveLGS(L, U, P, b);

    cout << "L =\n" << L << "\n\n";
    cout << "U =\n" << U << "\n\n";
    cout << "P =\n" << P << "\n\n";

    cout << "x' =\n" << x << "\n\n";


    // Teilaufgabe c)
    cout << "___________c)____________\n\n";


    //// A ändert sich nicht, also können L, U und P wiederverwendet werden und es muss nur y' bestimmt werden.

    b << 1, 2*sqrt(3), 3;

    x = solveLGS(L, U, P, b);

    cout << "y' =\n" << x << "\n\n";



    // Teilaufgabe d)
    cout << "___________d)____________\n\n";

    A << 0, -0.5, 0.5,          0, sqrt(3)/2, sqrt(3)/2,            1, 0, 0;

    L = A.partialPivLu().matrixLU().triangularView<UpLoType::UnitLower>();
    U = A.partialPivLu().matrixLU().triangularView<UpLoType::Upper>();
    P = A.partialPivLu().permutationP();

    cout << "L =\n" << L << "\n\n";
    cout << "U =\n" << U << "\n\n";
    cout << "P =\n" << P << "\n\n";

    return 0;
}
