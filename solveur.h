#ifndef SOLVEUR_H
#define SOLVEUR_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

// Fonction de décomposition de Cholesky
void chol(const Eigen::MatrixXd& A, Eigen::MatrixXd& L);
void reschol(MatrixXd& L, const VectorXd& b, VectorXd& x);
void gradientConjugue(const MatrixXd& A, const VectorXd& b, VectorXd& x,VectorXd& x_0, int Nmax , double eps );

#endif // SOLVEUR_H

