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


// void readMsh(const string& filename, 
//              vector<double>& nodeX, vector<double>& nodeY, 
//              vector<vector<int>>& elements) {
//     // Simule la lecture, adapter selon le format .msh
//     nodeX = {0.0, 1.0, 0.0};
//     nodeY = {0.0, 0.0, 1.0};
//     elements = {{1, 2, 3}};  // Triangle formé par les nœuds 1, 2, 3
// }


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
pair <MatrixXd,double> computeB(MatrixXd& nodes, MatrixXd& N) {

    //calcul de J

    MatrixXd J= N*nodes;
    double detJ= J(0,0)*J(1,1)-J(0,1)*J(1,0);

    //construction de T
    MatrixXd J_1= J.inverse();
    MatrixXd T= J_1*N;

    // Coefficients pour les dérivées des fonctions de forme
    double b1 = T(0,0), b2 = T(0,1), b3 = T(0,2);
    double c1 = T(1,0), c2 = T(1,1), c3 = T(1,2);

    // Construction de la matrice B
    MatrixXd B (3,6);
    B<< b1, 0, b2 , 0, b3 , 0,
        0, c1 , 0, c2, 0, c3 ,
        c1 , b1 , c2 , b2 , c3, b3 ;

    
    return {B,detJ};
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

pair <MatrixXd,VectorXd> computeKe_Fe( const MatrixXd& B,  const MatrixXd& D,  double area, double detJ, VectorXi tri, double g ,double rho) 
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

    // Calcul de BtD * B
    MatrixXd Ke(cols,rows );
    Ke = BtD*B;

    Ke= Ke*area*detJ;

    VectorXd Fe(6);
    Fe << 0,1,0,1,0,1;

    Fe=-rho*g*area*Fe/3;

    return {Ke,Fe};
}

VectorXd computeFe( const MatrixXd& B, Mesh2D* _msh)
{
    size_t rows=B.rows();
    VectorXd Fe(rows);
    VectorXd I(rows);
    double alpha(0);
    string BC;
    
	for (unsigned int i = 0; i < _msh->Get_edges().size(); i++)
    {

        int t1=_msh->Get_edges()[i].Get_T1();
		int t2=_msh->Get_edges()[i].Get_T2();  

        if (t2!=-1)     
        {
            Fe(i)=0;
        }

        else
        {   

            double g=9.81;
            double rho=2000.;
            double mu=0.25;
            double h=25.;
            double L=1;
            double w=1000;
            BC=_msh->Get_edges()[i].Get_BC();

            if (BC=="Neumann_homogene")
            {
                alpha=-rho*g*pow(h,2)*L/24.;
                for (int i=0;i<6;i+=2)
                {
                    I(i)=0;
                    I(i+1)=1;
                }

                Fe=alpha*I;


            }
            else if (BC=="Neumann")
            {	
                alpha=-rho*g*pow(h,2)*L/24.+w*g*pow(h,2)*L;
            }

            else if (BC=="Dirichlet")
            {

            }
        }
    }

    return Fe;
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

    // Exemple : propriétés du matériau
    double E = 15e9; // Module de Young en Pascals
    double nu = 0.25;  // Coefficient de Poisson3,84 euros bru
    double g=9.81;
    double rho=2000;

    // recuperation du maillage
    mesh->Read_mesh(data_file->Get_mesh_name());

    const vector<Triangle>& triangles = mesh->Get_triangles();
    const vector<Vertex>& vertices =mesh->Get_vertices();
    MatrixXi table_corresp (mesh->Get_triangles().size(),3);
    MatrixXd N(2,3);
    N << -1, 1,0,
        -1,0,1;
 
    // Coordonnées des nœuds d'un élément triangulaire
    MatrixXd nodes (3,2);

    // Calcul de la matrice D
    MatrixXd D = computeD(E, nu);
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Matrice D:" << std::endl;
    std::cout << "-------------------------------------" << std::endl;

    //Definition de la table de correspondance
    mesh->Build_Bool();
    mesh->Build_Table();
    MatrixXi Table= mesh->Get_Table_degre();
        cout <<"help"<< endl;
    cout << Table << endl;
    int taille=Table.maxCoeff();
        cout <<"help"<< endl;
    //definition de la taille de K et F
    MatrixXd K(taille,taille);
    VectorXd F (taille);
    VectorXi Table_correspondance_locale(6);

    for (long unsigned int i=0;i<mesh->Get_triangles().size();i++)
    {
        Vector3i tri = triangles[i].Get_vertices();
        cout <<tri<<endl;
        table_corresp.row(i) =tri;


        //coordonnées réel des noeuds 

        nodes(0,0) = vertices[i].Get_coor()(0), nodes(0,1)=vertices[i].Get_coor()(1);
        nodes(1,0) = vertices[i].Get_coor()(0), nodes(1,1)=vertices[i].Get_coor()(1);
        nodes(2,0) = vertices[i].Get_coor()(0), nodes(2,1)=vertices[i].Get_coor()(1);
        // Calcul de la matrice B et de l'aire
        cout << "i" << endl;
        auto results =computeB(nodes,N);
        MatrixXd B =  results.first;
        double detJ= results.second;
        std::cout << "-------------------------------------" << std::endl;
        std::cout << "Matrice B:" << std::endl;
        std::cout << "-------------------------------------" << std::endl;

        double area=mesh->Get_triangles_area()(i);
        std::cout << "Aire du triangle: " << area << std::endl;
        std::cout << "-------------------------------------" << std::endl;

        // Calcul de la matrice de rigidité élémentaire Ke

        auto res = computeKe_Fe(B, D, area,detJ,tri,g,rho);
        MatrixXd Ke = res.first;
        VectorXd Fe = res.second;

        //construction table de correspondance locale

        Table_correspondance_locale <<Table(tri[0],0),Table(tri[0],1),Table(tri[1],0),Table(tri[1],1),Table(tri[2],0),Table(tri[2],1);

        //calcul de la matrice de rigidité globale
         for (int k=0; k<3;k++){
            for (int l=0; l<3; l++){

                if (Table_correspondance_locale(k) && Table_correspondance_locale(l) != -1){
                    K(Table_correspondance_locale[k],Table_correspondance_locale(l))+= Ke(k,l);
                    F(Table_correspondance_locale(k))+= Fe(k);
                }
            }
         }

        std::cout << "-------------------------------------" << std::endl;
        std::cout << "Matrice de rigidité élémentaire Ke:" << std::endl;
        std::cout << "-------------------------------------" << std::endl;
        //printMatrix(Ke);
        //nodes= Array32d::Zero();

    }
    cout << F << endl;

    delete mesh;
    delete data_file;

    return 0;
}

