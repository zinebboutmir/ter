#ifndef SOLVEUR_H
#define SOLVEUR_H

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Fonction de décomposition de Cholesky
void chol(const Eigen::MatrixXd& A, Eigen::MatrixXd& L);
void reschol(MatrixXd& L, const VectorXd& b, VectorXd& x);
void gradientConjugue(const MatrixXd& A, const VectorXd& b, VectorXd& x, const VectorXd& x0, int Nmax = 1000, double eps = 1e-6);

#endif // SOLVEUR_H

