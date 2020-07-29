// Calculate curvature map of triangular mesh
// Jun-Seop Song

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

#define PI 3.141592653589793238

typedef struct Node {
	int idx;
	double xyz[3]; // coordinate (x,y,z)
	double curvature; // Gaussian curvature
	double sumArea; // sum of areas of neighboring triangles
} Node;

typedef struct Element {
	Node* v[3];
} Element;


int numVertex, numEle;
vector<Node> vertex;
vector<Element> ele;


void importMesh()
{
	ifstream fin;
	char buf[100];
	int i, t, u[3];

	fin.open("input.ply");

	fin.getline(buf, 100);
	fin.getline(buf, 100);
	fin >> buf >> buf >> numVertex >> buf;
	fin.getline(buf, 100);
	fin.getline(buf, 100);
	fin.getline(buf, 100);
	fin >> buf >> buf >> numEle >> buf;
	fin.getline(buf, 100);
	fin.getline(buf, 100);

	cout << "Vertex:" << numVertex << "   Element:" << numEle << endl;

	vertex.assign(numVertex, Node());
	for (i = 0; i < numVertex; ++i)
	{
		vertex[i].idx = i;
		fin >> vertex[i].xyz[0] >> vertex[i].xyz[1] >> vertex[i].xyz[2];
	}

	ele.assign(numEle, Element());
	for (i = 0; i < numEle; ++i)
	{
		fin >> t >> u[0] >> u[1] >> u[2];
		ele[i].v[0] = &vertex[u[0]];
		ele[i].v[1] = &vertex[u[1]];
		ele[i].v[2] = &vertex[u[2]];
	}

	fin.close();
}


double length(Node* node1, Node* node2)
{
	double len = 0.0;
	for (int i = 0; i < 3; ++i) len += (node1->xyz[i] - node2->xyz[i])*(node1->xyz[i] - node2->xyz[i]);
	return sqrt(len);
}


void calcGaussianCurvature()
{
	int i, j;
	Node* v[3];
	double len[3], s, area;

	for (i = 0; i < numVertex; ++i)
	{
		vertex[i].curvature = 2 * PI;
		vertex[i].sumArea = 0.0;
	}

	for (i = 0; i < numEle; ++i)
	{
		for (j = 0; j < 3; ++j) v[j] = ele[i].v[j];

		for (j = 0; j < 3; ++j) len[j] = length(v[(j + 1) % 3], v[(j + 2) % 3]);
		s = (len[0] + len[1] + len[2]) / 2.0;
		area = sqrt(s*(s - len[0])*(s - len[1])*(s - len[2]));

		for (j = 0; j < 3; ++j)
			v[j]->curvature -= acos((len[(j + 1) % 3] * len[(j + 1) % 3] + len[(j + 2) % 3] * len[(j + 2) % 3] - len[j] * len[j]) / (2.0*len[(j + 1) % 3] * len[(j + 2) % 3]));
		for (j = 0; j < 3; ++j) v[j]->sumArea += area;
	}

	// By Gauss-Bonnet Theorem
	for (i = 0; i < numVertex; ++i)
		vertex[i].curvature /= (vertex[i].sumArea / 3.0);
}


void exportResult()
{
	ofstream fout;
	int i;

	fout.open("output.ply");
	fout << "ply" << endl
		<< "format ascii 1.0" << endl
		<< "element vertex " << numVertex << endl
		<< "property float x" << endl
		<< "property float y" << endl
		<< "property float z" << endl
		<< "property float quality" << endl
		<< "element face " << numEle << endl
		<< "property list uchar int vertex_indices" << endl
		<< "end_header" << endl;
	for (i = 0; i < numVertex; ++i)
		fout << vertex[i].xyz[0] << " " << vertex[i].xyz[1] << " " << vertex[i].xyz[2] << " " << vertex[i].curvature << endl;
	for (i = 0; i < numEle; ++i)
		fout << "3 " << ele[i].v[0]->idx << " " << ele[i].v[1]->idx << " " << ele[i].v[2]->idx << endl;
	fout.close();
}


int main()
{
	importMesh();
	calcGaussianCurvature();
	exportResult();

	return 0;
}
