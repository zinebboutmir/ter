#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include "Mesh2D.h"
#include "data_file.h"


using namespace std;
using namespace Eigen;

// g++ -std=c++11 -o run main.cpp


// étape 1 : Fonction pour lire le fichier de maillage

void readMsh(const string& filename, 
             vector<double>& nodeX, vector<double>& nodeY, 
             vector<vector<int>>& elements) {
    // Simule la lecture, adapter selon le format .msh
    nodeX = {0.0, 1.0, 0.0};
    nodeY = {0.0, 0.0, 1.0};
    elements = {{1, 2, 3}};  // Triangle formé par les nœuds 1, 2, 3
}


//étape 2 : Fonction pour calculer la matrice D (propriétés du matériau)
MatrixXd computeD(double E, double nu) {
    double factor = E /( (1 - 2*nu)*(1+nu));
    Matrix3d D;
    D << (1-nu),  nu, 0,
        nu, 1-nu, 0,
        0, 0, (1 - 2*nu) / 2;
    D= D*factor;
    return  D;
}

// Fonction pour calculer la matrice B et l'aire d'un élément triangulaire
pair<MatrixXd , double> computeB(const vector<pair<double, double>>& nodes) {
    // Coordonnées des sommets du triangle
    double x1 = nodes[0].first, y1 = nodes[0].second;
    double x2 = nodes[1].first, y2 = nodes[1].second;
    double x3 = nodes[2].first, y3 = nodes[2].second;
    // nodes[1] fournit (x2,y2)

    // Calcul de l'aire du triangle
    double area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));

    // Coefficients pour les dérivées des fonctions de forme
    double b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
    double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;

    // Construction de la matrice B
    MatrixXd B (3,3);
    B<< b1, 0, b2 , 0, b3 , 0,
        0, c1 , 0, c2, 0, c3 ,
        c1 , b1 , c2 , b2 , c3, b3 ;

    B= B*2*area;
    
    return {B, area};
}

// Fonction pour afficher une matrice
void printMatrix(const vector<vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double val : row) {
            cout << val << "\t";
        }
        cout << endl;
    }
}

MatrixXd computeKe(
    const MatrixXd& B, 
    const MatrixXd& D, 
    double area) 
{
    size_t rows = B.rows();
    size_t cols = B.cols();
    size_t Dsize = D.rows();

    // Vérifier que les dimensions sont compatibles
    if (rows != Dsize) {
        throw runtime_error("Dimension mismatch between B and D matrices.");
    }

    // Calcul de Bt * D
    MatrixXd BtD(cols, Dsize);
    BtD = B.transpose()*D;
    // for (size_t i = 0; i < cols; ++i) {
    //     for (size_t j = 0; j < Dsize; ++j) {
    //         for (size_t k = 0; k < rows; ++k) {
    //             BtD[i][j] += B[k][i] * D[k][j];
    //         }
    //     }
    // }

    // Calcul de BtD * B
    MatrixXd Ke(cols,rows );
    Ke = BtD*B;
    // for (size_t i = 0; i < cols; ++i) {
    //     for (size_t j = 0; j < cols; ++j) {
    //         for (size_t k = 0; k < rows; ++k) {
    //             Ke[i][j] += BtD[i][k] * B[k][j];
    //         }
    //         // Multiplier par l'aire
    //         Ke[i][j] *= area;
    //     }
    // }
    Ke= Ke*area;
    return Ke;
}

int main(int argc, char** argv) {

    if (argc < 2)
   {
      cout << "Please, enter the name of your data file." << endl;
      cout << "Usage: " << argv[0] << " <file.toml>" << endl;
      exit(0);
   }

   const string data_file_name = argv[1];

   // ----------------------- Fichier de données --------------------------------
   DataFile* data_file = new DataFile(data_file_name);



    Mesh2D* mesh = new Mesh2D(data_file->Get_BC_ref(),data_file->Get_BC_type());

    // vector<double> sommet_droite=data_file->Get_vertices()(i,0);
    // vector<double> sommet_gauche=data_file->Get_vertices()(i,0)

    // Exemple : propriétés du matériau
    double E = 15e9; // Module de Young en Pascals
    double nu = 0.25;  // Coefficient de Poisson
    double g=9.81;


    mesh->Read_mesh(data_file->Get_mesh_name());
    // readMsh(fichier,sommet_droite,sommet_gauche,elements)

    MatrixXd D = computeD(E, nu);
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Matrice D:" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    //printMatrix(D);

    // Coordonnées des nœuds d'un élément triangulaire
    vector<pair<double, double>> nodes = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};

    // Calcul de la matrice B et de l'aire
   auto result = computeB(nodes);
    MatrixXd B = result.first;
    double area = result.second;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Matrice B:" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    //printMatrix(B);
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Aire du triangle: " << area << std::endl;
    std::cout << "-------------------------------------" << std::endl;

    // Calcul de la matrice de rigidité élémentaire Ke
    MatrixXd Ke = computeKe(B, D, area);
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Matrice de rigidité élémentaire Ke:" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    //printMatrix(Ke);

    return 0;
}

