#include "classes.h" 
#include <iostream> 
#include <cmath>     
#include <random>
#include <string>
#include <algorithm>
#include <list>
#include <stdio.h>
#include "liblbfgs/lbfgs.h"
static std::default_random_engine engine(9); 
static std::uniform_real_distribution<double> uniform(0., 1.);

void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none")
{
    FILE *f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (int i = 0; i < polygons.size(); i++)
    {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (int j = 0; j < polygons[i].vertices.size(); j++)
        {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}



Vector::Vector(double x, double y)
{
    coords[0] = x;
    coords[1] = y;
}

bool &Vector::operator==(const Vector &b)
{
    bool cond1 = coords[0] == b[0];
    bool cond2 = coords[1] == b[1];
    static bool res = cond1 && cond2;
    return res;
}

bool &Vector::operator!=(const Vector &b)
{
    bool cond1 = coords[0] != b[0];
    bool cond2 = coords[1] != b[1];
    static bool res = cond1 || cond2;
    return res;
}

Vector &Vector::operator+=(const Vector &b)
{
    coords[0] += b[0];
    coords[1] += b[1];
    return *this;
}

Vector &Vector::operator*=(const Vector &b)
{
    coords[0] *= b[0];
    coords[1] *= b[1];
    return *this;
}

Vector &Vector::operator/=(const Vector &b)
{
    coords[0] /= b[0];
    coords[1] /= b[1];
    return *this;
}

Vector &Vector::operator-=(const Vector &b)
{
    coords[0] -= b[0];
    coords[1] -= b[1];
    return *this;
}

const double &Vector::operator[](int i) const { return coords[i]; }
double &Vector::operator[](int i) { return coords[i]; }

Vector Vector::operator+(const Vector &a)
{
    return Vector(a[0] + coords[0], a[1] + coords[1]);
}

Vector Vector::operator+(const double a)
{
    return Vector(a + coords[0], a + coords[1]);
}

Vector Vector::operator-(const Vector &a)
{
    return Vector(coords[0] - a[0], coords[1] - a[1]);
}

Vector Vector::operator-(const double a)
{
    return Vector(coords[0] - a, coords[1] - a);
}

Vector Vector::operator*(const Vector &a)
{
    return Vector(a[0] * coords[0], a[1] * coords[1]);
}

Vector Vector::operator*(const double a)
{
    return Vector(a * coords[0], a * coords[1]);
}

Vector Vector::operator/(const Vector &a)
{
    return Vector(coords[0] / a[0], coords[1] / a[1]);
}

Vector Vector::operator/(const double a)
{
    return Vector(coords[0] / a, coords[1] / a);
}

double Vector::dot(const Vector &a)
{
    return a[0] * coords[0] + a[1] * coords[1];
}

Vector Vector::min(double s)
{
    return Vector(std::min(coords[0], s), std::min(coords[1], s));
}

Vector Vector::max(double s)
{
    return Vector(std::max(coords[0], s), std::max(coords[1], s));
}

Vector Vector::pow(double s)
{
    return Vector(std::pow(coords[0], s), std::pow(coords[1], s));
}

void Vector::print()
{
    std::cout << "[" << coords[0] << "," << coords[1] << "]" << std::endl;
}

double Vector::norm()
{
    return sqrt(dot(*this));
}

Vector Vector::normalize()
{
    return *this / norm();
}

/* ================================================================= 

                    Vectors

=================================================================== */

Edge::Edge(const Vector &a, const Vector &b)
{
    this->point_a = a;
    this->point_b = b;
}

Vector intersect(Vector &prevVertex, Vector &curVertex, Edge &clipEdge)
{
    Vector N = Vector(clipEdge.point_b[1] - clipEdge.point_a[1], clipEdge.point_a[0] - clipEdge.point_b[0]);
    Vector a_b = curVertex - prevVertex;
    double t = N.dot(clipEdge.point_a - prevVertex) / N.dot(a_b);
    Vector P = prevVertex + a_b * t;
    if (t < 0 || t > 1)
        return Vector(0., 0.);
    return P;
}

Polygon::Polygon()
{
}

double Polygon::area()
{
    double area = 0;
    int n = vertices.size();
    if (n == 0)
    {
        return 0;
    }

    for (int i = 0; i < n; i++)
    {
        Vector point1 = vertices[i];
        Vector point2 = (i < n - 1) ? vertices[i + 1] : vertices[0];
        area += point1[0] * point2[1] - point1[1] * point2[0];
    }
    
    return std::abs(0.5 * area);
}




bool inside(Vector &vertex, Edge &clipEdge)
{
    Vector N = Vector(clipEdge.point_b[1] - clipEdge.point_a[1], clipEdge.point_a[0] - clipEdge.point_b[0]) * (-1);
    return N.dot(vertex - clipEdge.point_a) <= 0;
}

Polygon clipPolygon(Polygon &subjectPolygon, Polygon &clipPolygon)
{
    Polygon outPolygon;
    for(int i = 0; i < clipPolygon.vertices.size(); i++) {
        Edge clipEdge(clipPolygon.vertices[i], clipPolygon.vertices[(i != 0) ? i - 1 : clipPolygon.vertices.size() - 1]);
        outPolygon = Polygon();
        for (size_t i = 0; i < subjectPolygon.vertices.size(); i++){
            Vector curVertex = subjectPolygon.vertices[i];
            Vector prevVertex = subjectPolygon.vertices[(i>0)?(i-1):subjectPolygon.vertices.size()-1];
            Vector intersection = intersect(prevVertex,curVertex,clipEdge);
            if (inside(curVertex,clipEdge)){
                if (!inside(prevVertex,clipEdge)) outPolygon.vertices.push_back(intersection);
                outPolygon.vertices.push_back(curVertex);
            }
            else if (inside(prevVertex,clipEdge)) outPolygon.vertices.push_back(intersection);
        }
        subjectPolygon = outPolygon;
    }
    return outPolygon;
}
double area(std::vector<Vector> vertices)
{

    double a = 0;
    int n = vertices.size();
    if (n == 0)
    { return a;}
    for (int i = 0; i < n; i++)
    {
        Vector point1 = vertices[i];
        Vector point2 = (i < n - 1) ? vertices[i + 1] : vertices[0];
        a += point1[0] * point2[1] - point1[1] * point2[0];
    }
    a/=2;
    return std::abs(a);
}
//computes triangulation
double triangul(std::vector<Vector> vertices, Vector point)
{
    if (vertices.size() == 0)
        return 0;
    Vector c0 = vertices[0];
    double res = 0;
    for (int i = 0; i < vertices.size() - 2; i++)
    {
        Vector c1 = vertices[i + 1];
        Vector c2 = vertices[i + 2];
        std::vector<Vector> triangle;
        triangle.push_back(c0);
        triangle.push_back(c1);
        triangle.push_back(c2);
        double a = std::abs(area(triangle));
        double b=((c0 - point).dot(c0 - point) + (c0 - point).dot(c1 - point) + (c0 - point).dot(c2 - point) +
                              (c1 - point).dot(c1 - point) + (c1 - point).dot(c2 - point) +
                              (c2 - point).dot(c2 - point));
        res += (a/ 6.) * b;
    }
    return std::abs(res);
}
Polygon clipPolygonLine(Polygon &subjectPolygon, Vector M, Vector Pi_Pj)
{
    int prevIndex;

    Polygon outPolygon = Polygon();
    Vector Pj_Pi = Pi_Pj * (-1);

    for (int j = 0; j < subjectPolygon.vertices.size(); j++)
    {
        Vector curVertex = subjectPolygon.vertices[j];
        if (j > 0)
            prevIndex = j - 1;
        else
            prevIndex = subjectPolygon.vertices.size() - 1;
        Vector prevVertex = subjectPolygon.vertices[prevIndex];
        Vector ab = curVertex - prevVertex;
        double t = Pi_Pj.dot(M - prevVertex) / Pi_Pj.dot(ab);
        Vector intersection = (t >= 0 && t <= 1) ? prevVertex + ab * t : Vector(0., 0.);
        if ((curVertex-M).dot(Pi_Pj) < 0)
        {
            if (!((prevVertex-M).dot(Pi_Pj) < 0))
            {

                outPolygon.vertices.push_back(intersection);
            }
            outPolygon.vertices.push_back(curVertex);
        }
        else if ((prevVertex-M).dot(Pi_Pj) < 0)
        {

            outPolygon.vertices.push_back(intersection);
        }
    }
    return outPolygon;
}

std::vector<Polygon> voronoi(Polygon &clipPolygon, std::vector<Vector> &points, const double *weights)
{
    Vector point1, point2, M, PiPj;
    double weight1, weight2;
    std::vector<Polygon> result;
    for (int i = 0; i < points.size(); i++)
    {
        point1 = points[i];
        weight1 = weights[i];
        Polygon outPolygon = clipPolygon;
        for (int j = 0; j < points.size(); j++)
        {
            point2 = points[j];
            weight2 = weights[j];
            if (i == j)
                continue;
            M = (point1 + point2) * 0.5;
            PiPj = (point2 - point1);
            M = M + PiPj * (weight1 - weight2) / (2 * pow(PiPj.norm(), 2.));
            outPolygon = clipPolygonLine(outPolygon, M, PiPj);
        }
        result.push_back(outPolygon);
    }
    return result;
}



void gradient_descent(Polygon &clipPolygon, std::vector<Vector> &points, const double *lambdas, double eps, double step, double *weights)
{
    double tmp, fx, grad_norm;
    Vector point;
    int N = points.size();
    double grad[N], area[N];
    std::vector<Polygon> polygons = voronoi(clipPolygon, points, lambdas);
    save_svg(polygons, "before_opti.svg");
    double error = 1;
    int k = 1;
    while (error > eps)
    {
        error = 0;
        fx = 0;
        grad_norm = 0;
        polygons = voronoi(clipPolygon, points, weights);

        for (int i = 0; i < N; i++)
        {
            area[i] = polygons[i].area();
            grad[i] = (lambdas[i] - area[i]);
            grad_norm += grad[i] * grad[i];
        }
        for (int i = 0; i < N; i++)
        {
            std::vector<Vector> vertices = polygons[i].vertices;
            point = points[i];
            tmp = (step / sqrt(grad_norm)) * grad[i];
            error += pow(tmp,2);
            weights[i] += tmp;
            fx += (triangul(vertices, point) - weights[i] * area[i] + lambdas[i] * weights[i]);
        }
        error = sqrt(error);
        printf("Iteration %d:\n", k);
        printf("  f_x = %f, x_[0] = %f, x_[1] = %f\n", fx, weights[0], weights[1]);
        printf("  error = %f,\n", error);
        printf("\n");
        k += 1;
        if (k % 5 == 0)
        {
            std::string filename = "optimization_" + std::to_string(k) + ".svg";
            save_svg(polygons, filename);
        }
    }
}
/* =============================================================

================================================================ */

class Optimal_transport
{
protected:
    lbfgsfloatval_t *m_x;

public:
    std::vector<Vector> points;
    Polygon clipPolygon;
    double *weights;
    double *lambdas;
    int max_iterations;
    std::vector<Polygon> polygons;

    Optimal_transport(Polygon &clipPolygon, std::vector<Vector> &points) : m_x(NULL)
    {
        this->clipPolygon = clipPolygon;
        this->points = points;
        max_iterations = 1000;
    }

    virtual ~Optimal_transport()
    {
        if (m_x != NULL)
        {
            lbfgs_free(m_x);
            m_x = NULL;
        }
    }

    int run()
    {
        lbfgsfloatval_t fx;
        int N = this->points.size();
        this->lambdas = (double *)malloc((N) * sizeof(double));

        lbfgsfloatval_t *m_x = lbfgs_malloc(N);
        if (m_x == NULL)
        {
            printf("ERROR: Failed to allocate a memory block for variables.\n");
            return 1;
        }

        /* Initialize the variables. */
        Vector C = Vector(0.5, 0.5);
        Vector diff;
        double total = 0;
        double maxi = 0;
        for (int i = 0; i < N; i++)
        {
            diff = C - this->points[i];
            lambdas[i] = std::exp(-pow(diff.norm(), 2.) / 0.02);
            total += lambdas[i];
        }
        for (int i = 0; i < N; i++)
        {
            lambdas[i] /= total;
            // m_x[i] = lambdas[i];
            if (lambdas[i] > maxi)
            {
                maxi = lambdas[i];
            }
            m_x[i] = 1;
        }
        std::cout << "max dirac:" << maxi << std::endl;
        this->polygons = voronoi(clipPolygon, this->points, lambdas);
        save_svg(polygons, "without_optimization.svg");

        //perform lbfgs
        lbfgs_parameter_t param;
        lbfgs_parameter_init(&param);
        param.max_iterations = max_iterations;

        int ret = lbfgs(N, m_x, &fx, _evaluate, _progress, this, &param);
        this->polygons = voronoi(clipPolygon, this->points, m_x);
        save_svg(polygons, "with_optimization.svg");

        /* result */
        printf("L-BFGS optimization terminated with status code = %d\n", ret);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, m_x[0], m_x[1]);
        free(lambdas);
        return ret;
    }

protected:
    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step)
    {
        return reinterpret_cast<Optimal_transport *>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step)
    {

        lbfgsfloatval_t fx = 0.0;
        this->polygons = voronoi(clipPolygon, this->points, x);

        for (int i = 0; i < n; i++)
        {

            std::vector<Vector> vertices = this->polygons[i].vertices;

            Vector point = this->points[i];
            double area = this->polygons[i].area();
            double tmp = triangul(vertices, point);
            fx += tmp - x[i] * area + this->lambdas[i] * x[i];
            g[i] = area - this->lambdas[i];
        }
        
        return -1*fx;
    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls)
    {
        return reinterpret_cast<Optimal_transport *>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls)
    {
        printf("Iteration %d:\n", k);
        printf("  f_x = %f, x_[0] = %f, x_[1] = %f\n", fx, x[0], x[1]);
        printf("  x_norm = %f, g_norm = %f, step = %f\n", xnorm, gnorm, step);
        printf("\n");
        return 0;
    }
};