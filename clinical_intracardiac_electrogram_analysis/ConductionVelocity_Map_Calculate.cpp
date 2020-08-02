/*
 * 3D CV map calculator
 * Developed by Jun-Seop Song (2017.02.28)
 * CV = 1/|grad(LAT)| (LAT: local activation time)
 */

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

class Node
{
public:
    double xyz[3];
    double lat;
    double grad_lat[3];
    double cv;
    double cv_vec[3];
    double area_sum; // sum of area of adjacent elements
    bool invalid;

    Node()
    {
        cv = 0.0;
        cv_vec[0] = cv_vec[1] = cv_vec[2] = 0.0;
        grad_lat[0] = grad_lat[1] = grad_lat[2] = 0.0;
        area_sum = 0.0;
        invalid = 0;
    }
};


int num_node, num_ele, **ele;
Node *node;


void inputLatMap()
{
    FILE *fid = fopen("LAT.dat", "r");
    int i, j, tmp;

    // Input num_node, num_ele
    fscanf(fid, "%d %d", &num_node, &num_ele);

    // Allocate
    node = new Node[num_node];
    ele = new int *[num_ele];
    for (i = 0; i < num_ele; ++i)
        ele[i] = new int[3];

    // Input xyz, lat
    for (i = 0; i < num_node; ++i)
        fscanf(fid, "%lf %lf %lf %lf", &node[i].xyz[0], &node[i].xyz[1], &node[i].xyz[2], &node[i].lat);

    // Input ele
    for (i = 0; i < num_ele; ++i)
    {
        for (j = 0; j < 3; ++j)
        {
            fscanf(fid, "%d", &tmp);
            ele[i][j] = tmp - 1;
        }
    }

    fclose(fid);
}


void outputCvMap()
{
    FILE *fid;
    int i;

    fid = fopen("CV_scalar.plt", "w");
    // num_node, num_ele
    fprintf(fid, "VARIABLES = \"X\", \"Y\", \"Z\", \"CV (m/s)\"\n");
    fprintf(fid, "ZONE F=FEPOINT, ET=triangle, N=%d , E=%d\n", num_node, num_ele);

    for (i = 0; i < num_node; ++i) // x, y, z, CV (scalar)
        fprintf(fid, "%lf %lf %lf %.3lf\n", node[i].xyz[0], node[i].xyz[1], node[i].xyz[2], node[i].cv);

    for (i = 0; i < num_ele; ++i) // element
        fprintf(fid, "%d %d %d\n", ele[i][0] + 1, ele[i][1] + 1, ele[i][2] + 1);

    fclose(fid);


    fid = fopen("CV_scalar.txt", "w");
    for (i = 0; i < num_node; ++i) // x, y, z, CV (scalar)
        fprintf(fid, "%lf %lf %lf %.3lf\n", node[i].xyz[0], node[i].xyz[1], node[i].xyz[2], node[i].cv);
    fclose(fid);


    fid = fopen("CV_vector.txt", "w");
    for (i = 0; i < num_node; ++i) // x, y, z, CV_x, CV_y, CV_z (vector)
        fprintf(fid, "%lf %lf %lf %lf %lf %lf\n", node[i].xyz[0], node[i].xyz[1], node[i].xyz[2],
                node[i].cv_vec[0], node[i].cv_vec[1], node[i].cv_vec[2]);
    fclose(fid);
}


// Calculate area and normal vector of element
std::pair<double, double *> calcAreaNormal(const double *xyz1, const double *xyz2, const double *xyz3)
{
    int i;
    double a[3], b[3], axb[3], axb_size = 0.0;

    for (i = 0; i < 3; ++i)
    {
        a[i] = xyz2[i] - xyz1[i];
        b[i] = xyz3[i] - xyz1[i];
    }

    for (i = 0; i < 3; ++i)
    {
        axb[i] = a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3];
        axb_size += (axb[i] * axb[i]);
    }
    axb_size = sqrt(axb_size);

    for (i = 0; i < 3; ++i)
        axb[i] /= axb_size;

    return std::make_pair(axb_size / 2.0, axb);
}


// grad(LAT) at each element = sum(LAT(i) * edgevector(i+1, i+2)) / (2 * area)
void calcCV()
{
    int e, i, j;
    std::pair<double, double *> area_normal; // area and normal vector of ele[e]
    double grad_lat_[3], grad_lat[3]; // grad(LAT) at ele[e]
    double grad_lat_abs; // |grad(LAT)| at ele[e]
    Node *v[3];

    // Calculate grad(LAT) for each ele[e]
    for (e = 0; e < num_ele; ++e)
    {
        // v[0], v[1], v[2]: nodes of ele[e] (for convenience)
        for (i = 0; i < 3; ++i)
            v[i] = &node[ele[e][i]];

        // Calculate area and normal vector of ele[e]
        area_normal = calcAreaNormal(v[0]->xyz, v[1]->xyz, v[2]->xyz);

        // Calculate grad(LAT)
        for (i = 0; i < 3; ++i) grad_lat_[i] = 0.0;
        for (i = 0; i < 3; ++i) // for each edge
        {
            for (j = 0; j < 3; ++j) // for each coordinate
            {
                grad_lat_[j] +=
                        v[(i + 2) % 3]->lat * (v[i]->xyz[j] - v[(i + 1) % 3]->xyz[j]) / (2.0 * area_normal.first);
            }
        }

        for (i = 0; i < 3; ++i)
            grad_lat[i] = grad_lat_[(i + 1) % 3] * area_normal.second[(i + 2) % 3] -
                          grad_lat_[(i + 2) % 3] * area_normal.second[(i + 1) % 3];

        // Calculate |grad(LAT)|
        grad_lat_abs = sqrt(grad_lat[0] * grad_lat[0] + grad_lat[1] * grad_lat[1] + grad_lat[2] * grad_lat[2]);


        // Interpolate grad(LAT) for each vertex (Step 1: weighted summation)
        for (i = 0; i < 3; ++i)
        {
            if (v[i]->invalid || grad_lat_abs < 0.0001) // invalid CV if grad(LAT) is too small
            {
                v[i]->invalid = 1;
                continue;
            }

            for (j = 0; j < 3; ++j) // for each coordinate
                v[i]->grad_lat[j] += area_normal.first * grad_lat[j];

            v[i]->area_sum += area_normal.first;
        }
    }

    // Interpolate grad(LAT) for each vertex (Step 2: average)
    for (i = 0; i < num_node; ++i)
    {
        if (!node[i].invalid)
        {
            for (j = 0; j < 3; ++j) // for each coordinate
                node[i].grad_lat[j] /= node[i].area_sum;

            grad_lat_abs = sqrt(node[i].grad_lat[0] * node[i].grad_lat[0] + node[i].grad_lat[1] * node[i].grad_lat[1] +
                                node[i].grad_lat[2] * node[i].grad_lat[2]);

            node[i].cv = 10.0 / grad_lat_abs; // (unit: m/s)

            for (j = 0; j < 3; ++j) // for each coordinate
                node[i].cv_vec[j] = node[i].grad_lat[j] / (grad_lat_abs * grad_lat_abs) * 10.0; // (unit: m/s)
        }
    }
}


int main()
{
    inputLatMap();

    calcCV();

    outputCvMap();

    return 0;
}