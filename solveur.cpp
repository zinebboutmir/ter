#include <iostream>
#include <Eigen/Dense>
#include <cmath>

#include"solveur.h"

using namespace std;
using namespace Eigen;

void chol(const MatrixXd& A, MatrixXd& L) {
    int n = A.rows();
    bool symmetric = true;

    // Vérifier si la matrice est symétrique
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (A(i, j) != A(j, i)) {
                symmetric = false;
                break;
            }
        }
        if (!symmetric) break;
    }

    if (!symmetric) {
        cout << "Matrice non symétrique" << endl;
        return;
    }

    
    L = MatrixXd::Zero(n, n);

    // Décomposition de Cholesky
    for (int i = 0; i < n; i++) {
        L(i, i) = A(i, i);
        
        // Calcul des éléments de la diagonale de L
        for (int k = 0; k < i; k++) {
            L(i, i) -= L(i, k) * L(i, k);
        }

        // Vérifier si on peut calculer la racine carrée
        if (L(i, i) < 0) {
            cout << "Problème de racine carrée" << endl;
            return;
        } else {
            L(i, i) = sqrt(L(i, i));
        }

        // Calcul des éléments en dessous de la diagonale de L
        for (int j = i + 1; j < n; j++) {
            L(j, i) = A(j, i);
            for (int k = 0; k < i; k++) {
                L(j, i) -= L(i, k) * L(j, k);
            }

            // Vérifier la division par zéro
            if (L(i, i) == 0) {
                cout << "Problème de division par zéro" << endl;
                return;
            } else {
                L(j, i) /= L(i, i);
            }
        }

        // Remettre les éléments au-dessus de la diagonale à zéro
        for (int j = 0; j < i; j++) {
            L(j, i) = 0;
        }
    }
}





