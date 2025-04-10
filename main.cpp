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
pair <MatrixXd,double> computeB(MatrixXd nodes, MatrixXd N) {

    //calcul de J

    MatrixXd J= N*nodes;
    double detJ= J(0,0)*J(1,1)-J(0,1)*J(1,0);
    //cout << detJ << endl;

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

pair <MatrixXd,VectorXd> computeKe_Fe( MatrixXd B,  MatrixXd D,  double area, double detJ, VectorXi tri, double g ,double rho, double L) 
{
    size_t rows = B.rows();
    size_t cols = B.cols();
    size_t Dsize = D.rows();


    // Vérifier que les dimensions sont compatibles
    if (rows != Dsize) {
        throw runtime_error("Dimension mismatch between B and D matrices.");
    }

    // Calcul de BtD * B
    MatrixXd Ke(cols,cols );
    Ke =L*B.transpose()*D*B*detJ/2;
    
    VectorXd Fe(6);
    Fe.setZero();
    Fe << 0,1,0,1,0,1;

    Fe=-rho*g*L*Fe*50*50/24;
    //cout << rho*g*50*50/24<< endl;
    //Fe=-rho*g*area*Fe/3;

    return {Ke,Fe};
}


double pression(double x,double rho, double g){
    return rho*g*(50-x)+10132;
}

VectorXd T(VectorXd S1,VectorXd S2, double eta){
    return (1-eta)*S1/2. +(1+eta)*S2/2.;
}

VectorXd Fct_Forme(double eta){
    VectorXd Forme (4);
    Forme.setZero();
    Forme(0)= (1-eta)/2; Forme(1) =(1-eta)/2; Forme(2) =(1+eta)/2 ; Forme(3)=(1+eta)/2;
    return Forme;
}

double Quadrature(double (*f)(double (*pression)(double,double,double),VectorXd (*T)(VectorXd S1,VectorXd S2, double), VectorXd (*Fct_Forme)(double ), double eta, VectorXd S1,VectorXd S2, int i),VectorXd S_1,VectorXd S_2,int i){
    double w1(1),w2(1);
    double x1(-1/(sqrt(3))),x2(1/(sqrt(3)));
    return w1*f(pression,T,Fct_Forme,x1,S_1,S_2,i)+w2*f(pression,T,Fct_Forme,x2,S_1,S_2,i);
    // double w(2), x(0);
    // return w*f(pression,T,Fct_Forme,x,S_1,S_2,i);
}

double f(double (*pression)(double,double,double),VectorXd (*T)(VectorXd S1,VectorXd S2, double), VectorXd (*Fct_Forme)(double ), double eta, VectorXd S1,VectorXd S2, int i){
    double g=9.81;
    double rho=1000.;
    return pression(T(S1,S2,eta)(1),rho,g)*Fct_Forme(eta)[i];

}

VectorXd computeFe_F_surf(  Mesh2D* _msh,const vector<Edge>& arete_bord, MatrixXi Table,const vector<Vertex>& vertices,int taille,double L)
{
    int rows=6;
    VectorXd Fe(rows);
    VectorXd F_surf(taille+1);

    string BC;
    int sommet1,sommet2,ddl1,ddl2,ddl3,ddl4;
    Vector4i tab;
    double x_norm,y_norm,longueur;
    Vector4d xy;

    Eigen::Matrix<double, Eigen::Dynamic, 2> _edg_normal=_msh->Get_edges_normal();
    Eigen::VectorXd length = _msh->Get_edges_length();
    
    VectorXd coor1,coor2,norm(2);
    Fe.setZero();
    F_surf.setZero();

    
	for (unsigned int i = 0; i < arete_bord.size(); i++)
    {

        BC=arete_bord[i].Get_BC();
        if (BC=="Neumann")
        {

            sommet1= arete_bord[i].Get_vertices()(0);
            sommet2= arete_bord[i].Get_vertices()(1);
            longueur= length[i];
            ddl1=Table(sommet1,0) ; ddl2= Table(sommet1,1); ddl3=Table(sommet2,0);ddl4=Table(sommet2,1);
            tab<< ddl1,ddl2,ddl3,ddl4;
            coor1=vertices[sommet1].Get_coor();
            coor2=vertices[sommet2].Get_coor();
            norm=_edg_normal.row(i);
            x_norm=norm(0);
            y_norm=norm(1);
            xy<<x_norm,y_norm,x_norm,y_norm;

        //     //double alpha=-rho*g*pow(h,2)*L/24.;
            for (int k=0;k<4;k++){
                Fe(k)+= -Quadrature(f,coor1,coor2,k)*xy(k)*longueur/2.;

            }
            for (int k(0); k<4;k++)
            {
                if (tab(k)!=-1)
                {
                    F_surf(tab(k))+=Fe(k);
                }
            }
        }
        //cout << "Fe surf: " << Fe << endl;
        //cout << "facteur : " << 1000*9.81*50*50/24 << endl;
        Fe.setZero();
    
    }
    F_surf= F_surf*L;
    return F_surf;
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
    double nu = 0.25;  // Coefficient de Poisson
    double g=9.81;
    double rho=2200;
    double L=100;

    // recuperation du maillage
    mesh->Read_mesh(data_file->Get_mesh_name());
    mesh-> Build_triangles_center_and_area();
    mesh-> Build_edges_normal_length_and_center();
    mesh-> Build_edges_bord();

    const vector<Triangle>& triangles = mesh->Get_triangles();
    const vector<Vertex>& vertices =mesh->Get_vertices();
    const vector<Edge>& arete = mesh->Get_edges();
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

    int taille=Table.maxCoeff();

    //definition de la taille de K et F
    MatrixXd K(taille+1,taille+1),Ke, B;
    VectorXd F (taille+1),Fe;
    VectorXi Table_correspondance_locale(6);
    K.setZero();
    F.setZero();
    Ke.setZero();
    Fe.setZero();

    Vector3i tri;

    for (long unsigned int i=0;i<mesh->Get_triangles().size();i++)
    {
        tri = triangles[i].Get_vertices();
        table_corresp.row(i) =tri;

        //coordonnées réel des noeuds 

        nodes(0,0) = vertices[tri(0)].Get_coor()(0), nodes(0,1)=vertices[tri(0)].Get_coor()(1);
        nodes(1,0) = vertices[tri(1)].Get_coor()(0), nodes(1,1)=vertices[tri(1)].Get_coor()(1);
        nodes(2,0) = vertices[tri(2)].Get_coor()(0), nodes(2,1)=vertices[tri(2)].Get_coor()(1);
        //cout << "nodes: " << nodes <<endl;

        // Calcul de la matrice B et de l'aire
        auto results =computeB(nodes,N);
        B =  results.first;
        double detJ= results.second;
        std::cout << "-------------------------------------" << std::endl;
        std::cout << "Matrice B:" << std::endl;
        std::cout << "-------------------------------------" << std::endl;
        cout << B << endl;

        double area=mesh->Get_triangles_area()(i);
        std::cout << "Aire du triangle: " << area << std::endl;
        std::cout << "-------------------------------------" << std::endl;

        // Calcul de la matrice de rigidité élémentaire Ke

        auto res = computeKe_Fe(B, D, area,detJ,tri,g,rho,L);
        Ke.setZero();
        Fe.setZero();
        Ke = res.first;
        Fe = res.second;
        cout << "Fe; " << Fe << endl;
        cout << "Ke ;" << Ke << endl;

        //construction table de correspondance locale

        Table_correspondance_locale <<Table(tri[0],0),Table(tri[0],1),Table(tri[1],0),Table(tri[1],1),Table(tri[2],0),Table(tri[2],1);
        cout << "table de corr : " << Table_correspondance_locale << endl;
        //calcul de la matrice de rigidité globale
         for (int k=0; k<6;k++){
            for (int l=0; l<6; l++){

                if ((Table_correspondance_locale(k) !=-1) && (Table_correspondance_locale(l) != -1)){
                    K(Table_correspondance_locale[k],Table_correspondance_locale[l])+= Ke(k,l);
                    //cout << " indice : " << k << " " << l<< endl;
                    //cout<<"K//////////"<<K<<endl;
                }
            }
         }

         for (int k(0); k<6; k++){
            if (Table_correspondance_locale(k) != -1){
                F(Table_correspondance_locale[k])+= Fe(k);
                //cout << "F table : "<<F(Table_correspondance_locale[k])<<endl;
            }
         }

        std::cout << "-------------------------------------" << std::endl;
        std::cout << "Matrice de rigidité élémentaire Ke:" << std::endl;
        std::cout << "-------------------------------------" << std::endl;
        //printMatrix(Ke);


    }
    
    cout<< "K: " << K<< endl;
    cout << "taille de K: " << K.size() << endl;
    F= F+computeFe_F_surf( mesh,arete, Table, vertices, taille, L);

    cout << "F: " << F<< endl;
    cout << "Table de corr :" << Table << endl;

    LLT<MatrixXd> cholesky(K);
    VectorXd q= cholesky.solve(F) ;

    cout << "q: " << q <<endl;

    delete mesh;
    delete data_file;

    return 0;
}

