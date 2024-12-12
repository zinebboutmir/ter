#include <iostream>
#include <vector>
//#include <Eigen/Dense>
#include <cmath>

// g++ -std=c++11 -o run main.cpp

// étape 1 : Fonction pour lire le fichier de maillage
void readMsh(const std::string& filename, 
             std::vector<double>& nodeX, std::vector<double>& nodeY, 
             std::vector<std::vector<int>>& elements) {
    // Simule la lecture, adapter selon le format .msh
    nodeX = {0.0, 1.0, 0.0};
    nodeY = {0.0, 0.0, 1.0};
    elements = {{0, 1, 2}};  // Triangle formé par les nœuds 0, 1, 2
}

//étape 2 : Fonction pour calculer la matrice D (propriétés du matériau)
std::vector<std::vector<double>> computeD(double E, double nu) {
    double factor = E / (1 - nu * nu);
    return {
        {factor, factor * nu, 0},
        {factor * nu, factor, 0},
        {0, 0, factor * (1 - nu) / 2}
    };
}

// Fonction pour calculer la matrice B et l'aire d'un élément triangulaire
std::pair<std::vector<std::vector<double>>, double> computeB(const std::vector<std::pair<double, double>>& nodes) {
    // Coordonnées des sommets du triangle
    double x1 = nodes[0].first, y1 = nodes[0].second;
    double x2 = nodes[1].first, y2 = nodes[1].second;
    double x3 = nodes[2].first, y3 = nodes[2].second;
    // nodes[1] fournit (x2,y2)

    // Calcul de l'aire du triangle
    double area = 0.5 * std::abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));

    // Coefficients pour les dérivées des fonctions de forme
    double b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
    double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;

    // Construction de la matrice B
    std::vector<std::vector<double>> B = {
        {b1 / (2 * area), 0, b2 / (2 * area), 0, b3 / (2 * area), 0},
        {0, c1 / (2 * area), 0, c2 / (2 * area), 0, c3 / (2 * area)},
        {c1 / (2 * area), b1 / (2 * area), c2 / (2 * area), b2 / (2 * area), c3 / (2 * area), b3 / (2 * area)}
    };

    return {B, area};
}

// Fonction pour afficher une matrice
void printMatrix(const std::vector<std::vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double val : row) {
            std::cout << val << "\t";
        }
        std::cout << std::endl;
    }
}

std::vector<std::vector<double>> computeKe(
    const std::vector<std::vector<double>>& B, 
    const std::vector<std::vector<double>>& D, 
    double area) 
{
    size_t rows = B.size();
    size_t cols = B[0].size();
    size_t Dsize = D.size();

    // Vérifier que les dimensions sont compatibles
    if (rows != Dsize) {
        throw std::runtime_error("Dimension mismatch between B and D matrices.");
    }

    // Calcul de Bt * D
    std::vector<std::vector<double>> BtD(cols, std::vector<double>(Dsize, 0.0));
    for (size_t i = 0; i < cols; ++i) {
        for (size_t j = 0; j < Dsize; ++j) {
            for (size_t k = 0; k < rows; ++k) {
                BtD[i][j] += B[k][i] * D[k][j];
            }
        }
    }

    // Calcul de BtD * B
    std::vector<std::vector<double>> Ke(cols, std::vector<double>(cols, 0.0));
    for (size_t i = 0; i < cols; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            for (size_t k = 0; k < rows; ++k) {
                Ke[i][j] += BtD[i][k] * B[k][j];
            }
            // Multiplier par l'aire
            Ke[i][j] *= area;
        }
    }

    return Ke;
}

int main() {
    // Exemple : propriétés du matériau
    double E = 210e9; // Module de Young en Pascals
    double nu = 0.3;  // Coefficient de Poisson

    // Calcul de la matrice D
    std::vector<std::vector<double>> D = computeD(E, nu);
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Matrice D:" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    printMatrix(D);

    // Coordonnées des nœuds d'un élément triangulaire
    std::vector<std::pair<double, double>> nodes = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};

    // Calcul de la matrice B et de l'aire
   auto result = computeB(nodes);
    std::vector<std::vector<double>> B = result.first;
    double area = result.second;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Matrice B:" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    printMatrix(B);
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Aire du triangle: " << area << std::endl;
    std::cout << "-------------------------------------" << std::endl;

    // Calcul de la matrice de rigidité élémentaire Ke
    std::vector<std::vector<double>> Ke = computeKe(B, D, area);
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Matrice de rigidité élémentaire Ke:" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    printMatrix(Ke);

    return 0;
}

// =======
// #include <Eigen/Dense>
// #include <cmath>

// #include"solveur.h"

// using namespace std;
// using namespace Eigen;


// int main()
// {
//     MatrixXd A(3, 3);
//     A << 4, 12, -16,
//          12, 37, -43,
//          -16, -43, 98;

//     int n = A.rows();
//     MatrixXd L(n, n);

//     chol(A, L);

//     // Afficher la matrice L
//     cout << "Matrice L*Transpose(L) (décomposition de Cholesky) :" << endl;
//     cout << L*L.transpose() << endl;

//     return 0;
// }