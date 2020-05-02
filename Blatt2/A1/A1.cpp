#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include "service.cpp"

using namespace std;
using namespace Eigen;

void saveMatrix(string filename, MatrixXd mat)
{
    ofstream file;
    file.open("build/"+filename+".txt");
    file << mat;
    file.close();
}

void approximate(MatrixXd U, VectorXd S, MatrixXd V, int k)
{
    MatrixXd A(512,512);
    A = A.setZero();
    
    for(int i = 0; i<k; i++)
    {
        A += S(i) * (U.col(i) * V.col(i).transpose());
    }

    saveMatrix("A"+to_string(k), A);
}


int main()
{

    //Bild laden
    MatrixXd bild(512,512), U(512,512), V(512,512);
    VectorXd S(512);

    loadData(bild, "Bild", 512, 512);
    saveMatrix("bild", bild);

    JacobiSVD<MatrixXd> svd(bild, ComputeFullV | ComputeFullU);

    U << svd.matrixU();
    V << svd.matrixV();
    S = svd.singularValues();

    approximate(U, S, V, 10);
    approximate(U, S, V, 20);
    approximate(U, S, V, 50);

    return 0;
}