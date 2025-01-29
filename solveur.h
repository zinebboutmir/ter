#ifndef SOLVEUR_H
#define SOLVEUR_H

#include <Eigen/Dense>

// Fonction de décomposition de Cholesky
void chol(const Eigen::MatrixXd& A, Eigen::MatrixXd& L);
void reschol(MatrixXd& L, const VectorXd& b, VectorXd& x);

#endif // SOLVEUR_H

