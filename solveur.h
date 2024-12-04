#ifndef SOLVEUR_H
#define SOLVEUR_H

#include <Eigen/Dense>

// Fonction de décomposition de Cholesky
void chol(const Eigen::MatrixXd& A, Eigen::MatrixXd& L);

#endif // SOLVEUR_H

