/*
Calculate Gaussian curvature of LA geometry
Jun-Seop Song 2015.07.14
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_DEGREE 15
#define PI 3.141592653589793238


typedef struct Node{
	bool isBlock;
	double xyz[3];

	int degree;
	int link[MAX_DEGREE];

	double Gaussian_Curvature;
	double sumArea; // sum of areas of neighboring triangles
} Node;


int nde, nel;
int *ele[3];
Node *mesh;


void importGSR()
{
	int i, j, tmpInt;
	
	FILE *fin = fopen("Input.dat", "r");
	fscanf(fin, "%d %d", &nde, &nel);
	
	mesh = new Node[nde];
	for(i=0; i<3; i++) ele[i] = new int[nel];
	
	for(i=0; i<nde; i++)
	{
		mesh[i].isBlock = false;
		mesh[i].degree = 0;
	}

	// Input vertex
	for(i=0; i<nde; i++)
	{
		for(j=0; j<3; j++) fscanf(fin, "%lf", &mesh[i].xyz[j]);
	}

	// Input element
	for(i=0; i<nel; i++)
	{
		for(j=0; j<3; j++) fscanf(fin, "%d", &ele[j][i]);
		fscanf(fin, "%d", &tmpInt);

		for(j=0; j<3; j++) ele[j][i]--;

		if(tmpInt == 4)
		{
			for(j=0; j<3; j++) mesh[ele[j][i]].isBlock = true;
		}
	}

	fclose(fin);
}


void addLink(int from, int to)
{
	int i;

	for(i=0; i<mesh[from].degree; i++)
	{
		if(mesh[from].link[i] == to) return;
	}

	mesh[from].link[mesh[from].degree++] = to;
}


void setLink()
{
	int i, j;

	for(i=0; i<nel; i++)
	{
		for(j=0; j<3; j++)
		{
			addLink(ele[j][i], ele[(j+1)%3][i]);
			addLink(ele[(j+1)%3][i], ele[j][i]);
		}
	}
}


double length(int node1, int node2)
{
	int i;
	double sum = 0.0;

	for(i=0; i<3; i++) sum += (mesh[node1].xyz[i]-mesh[node2].xyz[i])*(mesh[node1].xyz[i]-mesh[node2].xyz[i]);

	return sqrt(sum);
}


void calcGaussianCurvature()
{
	int i, j;
	int v1, v2, v3;
	double len12, len13, len23, s;

	for(i=0; i<nde; i++)
	{
		mesh[i].Gaussian_Curvature = 2*PI;
		mesh[i].sumArea = 0.0;
	}

	for(i=0; i<nel; i++)
	{
		for(j=0; j<3; j++)
		{
			v1 = ele[j][i];
			v2 = ele[(j+1)%3][i];
			v3 = ele[(j+2)%3][i];

			len12 = length(v1, v2);
			len13 = length(v1, v3);
			len23 = length(v2, v3);
			s = (len12 + len13 + len23) / 2.0;

			mesh[v1].Gaussian_Curvature -= acos((len12*len12 + len13*len13 - len23*len23)/(2.0*len12*len13));
			mesh[v1].sumArea += sqrt(s*(s-len12)*(s-len13)*(s-len23));
		}
	}

	// By Gauss-Bonnet Theorem
	for(i=0; i<nde; i++)
	{
		mesh[i].Gaussian_Curvature = mesh[i].Gaussian_Curvature / (mesh[i].sumArea/3.0);
	}
}


void exportResult()
{
	int i;
	FILE *fout;
	
	fout = fopen("Gaussian_Curvature.plt", "w");

	fprintf(fout, "VARIABLES = \"X\", \"Y\", \"Z\", \"Gaussian_Curvature\"\n");
	fprintf(fout, "ZONE F=FEPOINT, ET=triangle, N=%d , E=%d\n", nde, nel);

	for(i=0; i<nde; i++)
	{
		fprintf(fout, "%lf %lf %lf %lf\n", mesh[i].xyz[0], mesh[i].xyz[1], mesh[i].xyz[2], mesh[i].Gaussian_Curvature);
	}

	for(i=0; i<nel; i++) fprintf(fout, "%d %d %d\n", ele[0][i]+1, ele[1][i]+1, ele[2][i]+1);

	fclose(fout);


	fout = fopen("Gaussian_Curvature.txt", "w");
	for(i=0; i<nde; i++) fprintf(fout, "%lf\n", mesh[i].Gaussian_Curvature);
	fclose(fout);
}


void curvatureSmoothing()
{
	double *tmp = new double[nde];
	int i, j;

	for(i=0; i<nde; i++)
	{
		tmp[i] = mesh[i].Gaussian_Curvature;
		
		for(j=0; j<mesh[i].degree; j++) tmp[i] += mesh[mesh[i].link[j]].Gaussian_Curvature;

		tmp[i] = tmp[i] / (mesh[i].degree + 1);
	}

	for(i=0; i<nde; i++) mesh[i].Gaussian_Curvature = tmp[i];

	delete[] tmp;
}


int main()
{
	printf("Curvature Calculator (2015.07.14)\n");
	importGSR();
	setLink();
	calcGaussianCurvature();
	curvatureSmoothing();
	exportResult();

	return 0;
}