#include "vector.cpp"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <random>
#include <chrono>
using namespace std::chrono;
std::random_device rd;
std::mt19937 e {rd()};
std::uniform_real_distribution<float> dist(0, 1);
Vector random_direction(){
    double r1 = dist(e);
    double r2 = dist(e) ;
    double x = cos(2*M_PI*r1)*sqrt(r2*(1. - r2));
    double y = sin(2*M_PI*r1)*sqrt(r2*(1. - r2));
    double z = 1-2*r2;
    return Vector(x,y,z);
}
//Sliced Optimal Transport color transfer algorithm
void Color_transfer(std::vector<Vector>& I, std::vector<Vector>& M){
    int n = I.size();
    int nbiter = 30;
    for(int iter = 0; iter < nbiter; iter++){
        Vector v = random_direction();
        std::vector<std::pair<double, int> > projI(n);
        std::vector<std::pair<double, int> > projM(n);
        for(int i = 0; i<n; i++ ){
            projI[i] = {dot(I[i],v),i};
            projM[i] = {dot(M[i],v),i};
        }
        //sort according to dot product
        std::sort(projI.begin(), projI.end());
        std::sort(projM.begin(), projM.end());
        //advert initial point cloud
        for(int i = 0; i<n; i++){
            int j = projI[i].second;
            I[j]= I[j]+(projM[i].first - projI[i].first)*v;
        }
    }
}

int main() {
    int x,y,n;
    unsigned char *image1 = stbi_load("image.jpeg", &x, &y, &n, 3);
    unsigned char *image2 = stbi_load("model.jpeg", &x, &y, &n, 3);
    auto start = high_resolution_clock::now();
    std::vector<Vector> I(x*y), M(x*y);
    for (int i = 0; i < x*y; i++) {
        I[i] = Vector(image1[3*i+0],image1[3*i+1],image1[3*i+2]);
        M[i] = Vector(image2[3*i+0],image2[3*i+1],image2[3*i+2]);
    }
    Color_transfer(I, M);
    for (int i = 0; i < x*y; i++) {
        image1[3*i+0] = I[i][0];
        image1[3*i+1] = I[i][1];
        image1[3*i+2] = I[i][2];
    }

    stbi_write_jpg("out.jpg", x, y, n, image1, 120);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    duration = duration / 1000000;
    std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;
    return 0;
}