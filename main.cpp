#include <iostream>
#include <Eigen/Dense>
#include <cmath>

using namespace std;
using namespace Eigen;


int main
{
    MatrixXd A(3, 3);
    A << 4, 12, -16,
         12, 37, -43,
         -16, -43, 98;

    int n = A.rows();
    MatrixXd L(n, n);

    chol(A, L);

    // Afficher la matrice L
    cout << "Matrice L*Transpose(L) (décomposition de Cholesky) :" << endl;
    cout << L.Transpose()*L << endl;

    return 0;
}