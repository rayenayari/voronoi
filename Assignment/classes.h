#include <iostream>
#include <vector>

class Vector
{
public:
    explicit Vector(double x = 0., double y = 0.);
    bool &operator==(const Vector &b);
    bool &operator!=(const Vector &b);
    Vector &operator+=(const Vector &b);
    Vector &operator-=(const Vector &b);
    Vector &operator*=(const Vector &b);
    Vector &operator/=(const Vector &b);

    Vector operator+(const Vector &a);
    Vector operator+(const double a);

    Vector operator-(const Vector &a);
    Vector operator-(const double a);

    Vector operator*(const Vector &a);
    Vector operator*(const double a);

    Vector operator/(const Vector &a);
    Vector operator/(const double a);

    const double &operator[](int i) const;
    double &operator[](int i);

    double dot(const Vector &a);
    Vector cross_product(const Vector &a);
    Vector pow(const double a);
    Vector max(const double a);
    Vector min(const double a);
    int argmin();
    int get_longest();
    double norm();
    Vector normalize();

    void print();

private:
    double coords[3];
};
//Simpler with an Edge class : avoids confusion
class Edge{
    public:
        explicit Edge(const Vector &a,const Vector &b);
        Vector point_a;
        Vector point_b;
};

class Polygon {
    public:
        explicit Polygon();
        std::vector<Edge> edges;
        std::vector<Vector> vertices;
        double area();

};