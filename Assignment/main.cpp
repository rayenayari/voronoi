#include "classes.cpp" 
#include <iostream>
#include <random>

int main()
{

    Vector poly1[4] = {Vector(0., 0.), Vector(0., 1.), Vector(1., 1.), Vector(1., 0.)};
    Polygon clipPolygon = Polygon();
    for (int j = 0; j < 4; j++)
    {
        clipPolygon.vertices.push_back(poly1[j]);
    }

    clipPolygon.area();
    int n_vertices = 1000;

    Polygon subjectPolygon = Polygon();
    for (int i = 0; i < n_vertices; i++){
        double x = uniform(engine);
        double y = uniform(engine);
        subjectPolygon.vertices.push_back(Vector(x, y));
    }
    
    double lambdas[subjectPolygon.vertices.size()];
    double weights[subjectPolygon.vertices.size()];
    //center
    Vector C = Vector(0.5,0.5);
    double total;
    for (int i = 0; i < n_vertices; i++){
        Vector point = subjectPolygon.vertices[i];
        Vector diff = C-point;
        //same as in the lecture notes
        lambdas[i] = 0;
        total += lambdas[i];
    }
    for (int i = 0; i < n_vertices; i++)
        {
            lambdas[i] /= total;
            weights[i] = 0;
        }
    Optimal_transport optimizer = Optimal_transport(clipPolygon,subjectPolygon.vertices);
    optimizer.run();
    std::vector<Polygon> result = optimizer.polygons;
    gradient_descent(clipPolygon,subjectPolygon.vertices,lambdas,1e-4,1,weights);

    //std::vector<Polygon> result = voronoi(clipPolygon,subjectPolygon.vertices,weights);
    //result.push_back(subjectPolygon);
    //PART GALOUET MERIGOT

     save_svg(result, "result.svg");
}