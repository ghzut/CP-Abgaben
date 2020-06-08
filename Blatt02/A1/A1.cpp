#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include "service.cpp"

using namespace std;
using namespace Eigen;

//Funktion, die die Matrix in einer .txt Datei abspeichert
void saveMatrix(string filename, MatrixXd mat)
{
    ofstream file;
    file.open("build/"+filename+".txt");
    file << mat;
    file.close();
}


int main()
{
    MatrixXd bild(512,512), U(512,512), V(512,512);
    VectorXd S(512);

    //Bild laden
    loadData(bild, "Bild", 512, 512);
    saveMatrix("bild", bild);

    //SVD mithilfe von eigen bestimmen
    JacobiSVD<MatrixXd> svd(bild, ComputeFullV | ComputeFullU);

    U << svd.matrixU();
    V << svd.matrixV();
    S = svd.singularValues();
    
    //Approximationen für k=10, 20 und 50 bestimmen und als .txt abspeichern
    MatrixXd A(512,512);
    A = A.setZero();
    
    for(int i = 0; i<50; i++)
    {
        A += S(i) * (U.col(i) * V.col(i).transpose());

        if(i==9 || i==19 || i==49)
        {
            saveMatrix("A"+to_string(i+1), A);
        }
    }

    

    return 0;
}