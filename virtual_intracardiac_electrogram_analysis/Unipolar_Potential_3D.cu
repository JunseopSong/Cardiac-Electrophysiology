/*
Calculate unipolar electric potential from 3D simulation data
Jun-Seop Song
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>

#define MAX_DEG 15

int nde, nel, numElectrode, *ie[3], *degree, *linkNode;
double *xy[3], *weightEdge, *volt;
bool *ifscar;
std::vector<double> electrodeX, electrodeY, electrodeZ;


// Import nde, nel, xy, ie, ifscar
void importMesh()
{
	int i, j;
	double temp;
	FILE *datain;
	printf("Importing mesh data...\n");

	datain = fopen("cuda_aa.dat", "rb");
	fscanf(datain, "%d", &nde);
	fscanf(datain, "%d", &nde);
	printf("    Number of nodes: %d\n", nde);
	fclose(datain);

	volt = (double *)malloc(sizeof(double) * nde);

	datain = fopen("cuda_xy.dat", "rb");
	for (i = 0; i < 3; ++i) xy[i] = (double *)malloc(sizeof(double) * nde);
	for (i = 0; i < nde; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			fscanf(datain, "%lf", &temp);
			xy[j][i] = temp;
		}
	}
	fclose(datain);

	datain = fopen("cuda_ie.dat", "rb");
	fscanf(datain, "%lf", &temp);
	nel = (int)temp;
	for (i = 0; i < 3; ++i) ie[i] = (int *)malloc(sizeof(int) * nel);
	for (i = 0; i < nel; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			fscanf(datain, "%lf", &temp);
			ie[j][i] = (int)temp;
		}
	}
	fclose(datain);

	datain = fopen("cuda_ifscar.dat", "rb");
	ifscar = (bool *)malloc(sizeof(bool) * nde);
	for (i = 0; i < nde; ++i)
	{
		fscanf(datain, "%lf", &temp);
		ifscar[i] = (bool)((int)temp);
	}
	fclose(datain);
}


double getDouble(char *strBuffer)
{
	char buffer[20];
	int cnt = 0;
	while (*strBuffer)
	{
		if (isdigit(*strBuffer) || *strBuffer == '.' || *strBuffer == '-')
			buffer[cnt++] = *strBuffer;
		else if (cnt != 0) break;
		strBuffer++;
	}
	buffer[cnt] = NULL;
	return atof(buffer);
}


void importElectrodeXYZ()
{
	std::ifstream in;
	const unsigned int maxBufferSize = 100000;
	char buffer[maxBufferSize], fileName[100], *token;
	bool inX, inY, inZ;

	printf("Importing electrode xyz...\n");

	for (int fileId = 1; fileId <= 4; ++fileId)
	{
		inX = inY = inZ = false;
		sprintf(fileName, "DxL_%d.csv", fileId);

		in.open(fileName, std::ifstream::in);
		if (!in.is_open()) continue;

		while (!(inX & inY & inZ))
		{
			in.getline(&buffer[0], maxBufferSize);
			if (buffer[0] == 'r' && buffer[1] == 'o' && buffer[2] == 'v' && buffer[3] == 'i')
			{
				if (!inX && buffer[7] == 'x')
				{
					inX = true;
					token = strtok(&buffer[10], ",");
					while (token != NULL)
					{
						electrodeX.push_back(getDouble(token));
						token = strtok(NULL, ",");
					}
				}
				else if (!inY && buffer[7] == 'y')
				{
					inY = true;
					token = strtok(&buffer[10], ",");
					while (token != NULL)
					{
						electrodeY.push_back(getDouble(token));
						token = strtok(NULL, ",");
					}
				}
				else if (!inZ && buffer[7] == 'z')
				{
					inZ = true;
					token = strtok(&buffer[10], ",");
					while (token != NULL)
					{
						electrodeZ.push_back(getDouble(token));
						token = strtok(NULL, ",");
					}
				}
			}
		}

		in.close();
	}

	numElectrode = electrodeX.size();
	printf("    Number of electrodes: %d\n", numElectrode);
}


void setIFDM()
{
	int i, j, k, k1, k2, v1, v2;
	bool ifNewNode;
	double temp, tmp1, tmp2, tmp3, *medianLine[3];

	int *deg; // degree of each node
	int *link[MAX_DEG]; // link[j][i] : j'th nodes linked to i
	double *weight[MAX_DEG]; // weight[j][i] : weight between i and link[j][i]
	double *lenTri[3]; // length of sides of each triangle(element)

	printf("Processing mesh data...\n");

	// Allocate
	deg = (int *)malloc(sizeof(int) * nde);
	for (i = 0; i < MAX_DEG; ++i) link[i] = (int *)malloc(sizeof(int) * nde);
	for (i = 0; i < MAX_DEG; ++i) weight[i] = (double *)malloc(sizeof(double) * nde);
	for (i = 0; i < 3; ++i) lenTri[i] = (double *)malloc(sizeof(double) * nel);

	for (i = 0; i < nde; ++i) deg[i] = 0;
	for (i = 0; i < 3; ++i)
		medianLine[i] = (double *)malloc(sizeof(double) * nel); // length of median lines of each triangle


	// ------------------ Geometry ------------------
	for (i = 0; i < nel; ++i)
	{
		// Calculate length of sides of triangle(element)
		for (j = 0; j < 3; ++j)
		{
			v1 = ie[(j + 1) % 3][i] - 1;
			v2 = ie[(j + 2) % 3][i] - 1;

			temp = 0.0;
			for (k = 0; k < 3; ++k) temp += (xy[k][v1] - xy[k][v2]) * (xy[k][v1] - xy[k][v2]);

			lenTri[j][i] = sqrt(temp);
		}

		// Calculate length of median lines
		for (j = 0; j < 3; ++j)
		{
			tmp1 = lenTri[j][i];
			tmp2 = lenTri[(j + 1) % 3][i];
			tmp3 = lenTri[(j + 2) % 3][i];
			medianLine[j][i] = sqrt(tmp2 * tmp2 / 2.0 + tmp3 * tmp3 / 2.0 - tmp1 * tmp1 / 4.0);
		}
	}


	// ------------------ Get graph structure ------------------
	for (i = 0; i < nel; ++i)
	{
		// Scan 1
		for (k1 = 0; k1 < 3; ++k1)
		{
			k2 = (k1 + 1) % 3;
			v1 = ie[k1][i] - 1;
			v2 = ie[k2][i] - 1;

			ifNewNode = true;
			for (j = 0; j < deg[v1]; ++j)
			{
				if (link[j][v1] == v2)
				{
					ifNewNode = false;
					break;
				}
			}

			if (ifNewNode)
			{
				link[deg[v1]][v1] = v2;
				deg[v1]++;
			}
		}

		// Scan 2
		for (k1 = 0; k1 < 3; ++k1)
		{
			k2 = (k1 + 2) % 3;
			v1 = ie[k1][i] - 1;
			v2 = ie[k2][i] - 1;

			ifNewNode = true;
			for (j = 0; j < deg[v1]; ++j)
			{
				if (link[j][v1] == v2)
				{
					ifNewNode = false;
					break;
				}
			}

			if (ifNewNode)
			{
				link[deg[v1]][v1] = v2;
				deg[v1]++;
			}
		}
	}


	// ------------------ Calculate weight ------------------
	for (i = 0; i < nde; ++i)
		for (j = 0; j < deg[i]; ++j)
			weight[j][i] = 0.0;

	for (i = 0; i < nel; ++i)
	{
		for (k1 = 0; k1 < 3; ++k1)
		{
			k2 = (k1 + 1) % 3;
			k = (k1 + 2) % 3;
			v1 = ie[k1][i] - 1;
			v2 = ie[k2][i] - 1;

			// find index for v1
			for (j = 0; j < deg[v1]; ++j)
			{
				if (link[j][v1] == v2)
				{
					weight[j][v1] += (medianLine[k][i] / lenTri[k][i] / 3.0);
					break;
				}
			}

			// find index for v2
			for (j = 0; j < deg[v2]; ++j)
			{
				if (link[j][v2] == v1)
				{
					weight[j][v2] += (medianLine[k][i] / lenTri[k][i] / 3.0);
					break;
				}
			}
		}
	}


	// ------------------ Transform 2D array into 1D array ------------------
	degree = (int *)malloc(sizeof(int) * (nde + 1));
	degree[0] = 0;
	for (i = 1; i <= nde; ++i) degree[i] = degree[i - 1] + deg[i - 1];

	linkNode = (int *)malloc(sizeof(int) * degree[nde]);
	for (i = 0; i < nde; ++i)
	{
		for (j = degree[i]; j < degree[i + 1]; ++j) linkNode[j] = link[j - degree[i]][i];
	}

	weightEdge = (double *)malloc(sizeof(double) * degree[nde]);
	for (i = 0; i < nde; ++i)
	{
		for (j = degree[i]; j < degree[i + 1]; ++j) weightEdge[j] = weight[j - degree[i]][i];
	}


	free(deg);
	for (i = 0; i < MAX_DEG; ++i) free(link[i]);
	for (i = 0; i < MAX_DEG; ++i) free(weight[i]);
	for (i = 0; i < 3; ++i) free(lenTri[i]);
	for (i = 0; i < 3; ++i) free(medianLine[i]);
}


// Calculate Laplacian(V) * dA
__global__ void calcLapl(double *lapl_, const double *volt_, const bool *nonExcitable,
	const int *linkNode_, const int *degree_, const double *weight_, const int numNde)
{
	const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
	int i;
	double volt_id, sum;

	if (id < numNde)
	{
		if (nonExcitable[id])
		{
			lapl_[id] = 0.0;
			return;
		}

		volt_id = volt_[id];
		sum = 0.0;
		for (i = degree_[id]; i < degree_[id + 1]; ++i)
		{
			if (!nonExcitable[linkNode_[i]]) sum += (volt_[linkNode_[i]] - volt_id) * weight_[i];
		}

		lapl_[id] = sum;
	}
}


__global__ void calcInvDist(double *invDist_, const double eleX, const double eleY, const double eleZ,
	const double *nodeX, const double *nodeY, const double *nodeZ, const int numNde)
{
	const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;

	if (id < numNde)
	{
		invDist_[id] = rsqrt((eleX - nodeX[id]) * (eleX - nodeX[id])
			+ (eleY - nodeY[id]) * (eleY - nodeY[id])
			+ (eleZ - nodeZ[id]) * (eleZ - nodeZ[id]));
	}
}


int main()
{
	int i, tt, tt_start, tt_end;
	double *lapl_d, *invDist_d, tmpPotential, *unipolarPotential, *volt_d, *weight_d, *nodeX_d, *nodeY_d, *nodeZ_d;
	bool *ifscar_d;
	int *degree_d, *linkNode_d;

	printf("Start: "); scanf("%d", &tt_start);
	printf("End: "); scanf("%d", &tt_end);

	importMesh();
	importElectrodeXYZ();
	setIFDM();

	cublasHandle_t cublasHandle = 0;
	cublasCreate(&cublasHandle);
	cudaMalloc((void **)&lapl_d, sizeof(double) * nde);
	cudaMalloc((void **)&invDist_d, sizeof(double) * nde);
	cudaMalloc((void **)&volt_d, sizeof(double) * nde);
	cudaMalloc((void **)&ifscar_d, sizeof(bool) * nde);
	cudaMalloc((void **)&degree_d, sizeof(int) * (nde + 1));
	cudaMalloc((void **)&linkNode_d, sizeof(int) * degree[nde]);
	cudaMalloc((void **)&weight_d, sizeof(double) * degree[nde]);
	cudaMalloc((void **)&nodeX_d, sizeof(double) * nde);
	cudaMalloc((void **)&nodeY_d, sizeof(double) * nde);
	cudaMalloc((void **)&nodeZ_d, sizeof(double) * nde);
	cudaMemcpy(ifscar_d, ifscar, sizeof(bool) * nde, cudaMemcpyHostToDevice);
	cudaMemcpy(degree_d, degree, sizeof(int) * (nde + 1), cudaMemcpyHostToDevice);
	cudaMemcpy(linkNode_d, linkNode, sizeof(int) * degree[nde], cudaMemcpyHostToDevice);
	cudaMemcpy(weight_d, weightEdge, sizeof(double) * degree[nde], cudaMemcpyHostToDevice);
	cudaMemcpy(nodeX_d, xy[0], sizeof(double) * nde, cudaMemcpyHostToDevice);
	cudaMemcpy(nodeY_d, xy[1], sizeof(double) * nde, cudaMemcpyHostToDevice);
	cudaMemcpy(nodeZ_d, xy[2], sizeof(double) * nde, cudaMemcpyHostToDevice);


	// Calculate electrical potential
	printf("Start!!!\n");
	unipolarPotential = (double *)malloc(sizeof(double) * numElectrode);
	FILE *fin;
	FILE *fout = fopen("unipolar_potential.dat", "wb");
	char finName[100];

	for (tt = tt_start; tt <= tt_end; ++tt)
	{
		printf("    %d\n", tt);

		sprintf(finName, "VoltBinary_%d.txt", tt);
		fin = fopen(finName, "rb");
		fread(volt, sizeof(double), nde, fin);
		fclose(fin);

		cudaMemcpy(volt_d, volt, sizeof(double) * nde, cudaMemcpyHostToDevice);
		calcLapl<<< (int)(nde / 256) + 1, 256 >>>(lapl_d, volt_d, ifscar_d, linkNode_d, degree_d, weight_d, nde);

		for (i = 0; i < numElectrode; ++i)
		{
			calcInvDist<<< (int)(nde / 256) + 1, 256 >>>(invDist_d, electrodeX[i], electrodeY[i], electrodeZ[i], nodeX_d, nodeY_d, nodeZ_d, nde);
			cublasDdot(cublasHandle, nde, lapl_d, 1, invDist_d, 1, &tmpPotential);
			//unipolarPotential[i] = tmpPotential;
			fprintf(fout, "%.10lf ", tmpPotential);
		}

		fprintf(fout, "\n");
		//fwrite(unipolarPotential, sizeof(double), numElectrode, fout);
	}

	fclose(fout);


	free(unipolarPotential);
	cublasDestroy(cublasHandle);
	cudaFree(nodeX_d);
	cudaFree(nodeY_d);
	cudaFree(nodeZ_d);
	cudaFree(lapl_d);
	cudaFree(invDist_d);
	cudaFree(volt_d);
	cudaFree(ifscar_d);
	cudaFree(degree_d);
	cudaFree(linkNode_d);
	cudaFree(weight_d);

	return 0;
}
