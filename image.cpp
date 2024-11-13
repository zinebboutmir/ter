#include <iostream>
#include <opencv2/opencv.hpp>

void processImage(const std::string& imagePath) {
    // Charger l'image depuis le fichier
    cv::Mat image = cv::imread(imagePath, cv::IMREAD_COLOR);
    
    // Vérifier si l'image a été chargée correctement
    if (image.empty()) {
        std::cerr << "Erreur : Impossible de charger l'image à partir de " << imagePath << std::endl;
        return;
    }

    // Réduire le bruit
    cv::Mat image_lissée;
    cv::GaussianBlur(image, image_lissée, cv::Size(5, 5), 1.5);

    // Segmentation : convertir l'image en niveaux de gris
    cv::Mat image_grise;
    cv::cvtColor(image_lissée, image_grise, cv::COLOR_BGR2GRAY);
    
    // Appliquer un seuillage pour obtenir une image binaire
    cv::Mat image_binaire;
    cv::threshold(image_grise, image_binaire, 128, 255, cv::THRESH_BINARY);

    // Détection des contours
    std::vector<std::vector<cv::Point>> contours;
    cv::findContours(image_binaire, contours, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE);

    // Dessiner les contours sur une nouvelle image
    cv::Mat image_contours = cv::Mat::zeros(image.size(), CV_8UC3);
    for (size_t i = 0; i < contours.size(); i++) {
        cv::drawContours(image_contours, contours, static_cast<int>(i), cv::Scalar(0, 255, 0), 2);
    }

    // Redimensionner les fenêtres pour s'assurer qu'elles ne soient pas trop grandes
    cv::resizeWindow("Image originale", 800, 600);
    cv::resizeWindow("Image binaire", 800, 600);
    cv::resizeWindow("Contours", 800, 600);

    // Afficher l'image originale, l'image binaire et l'image avec contours
    cv::imshow("Image originale", image);
    cv::imshow("Image binaire", image_binaire);
    cv::imshow("Contours", image_contours);

    // Attendre qu'une touche soit pressée avant de fermer les fenêtres
    cv::waitKey(0); // Attendre une touche
    cv::destroyAllWindows(); // Fermer toutes les fenêtres
}

int main() {
    // Spécifie le chemin de l'image star.jpg
    std::string imagePath = "x_essai4_x50_0.07um_x50_00.png";

    // Traiter l'image
    processImage(imagePath);

    return 0;
}

