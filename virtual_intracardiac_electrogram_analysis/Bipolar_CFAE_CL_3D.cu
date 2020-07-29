/*
Bipolar EGM, CFAE-CL Calculator for 3D Simulation
Jun-Seop Song 2014.12.18
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>

float REFRACTORY = 40.0, WIDTH = 15.0, PPSENSITIVITY = 2.0;


int nde, nel, start, end, sampleNum;
double *vm;
float **egm, *CL;
int *num1, *num2, *node1, *node2;
bool *ifScar;


// MATLAB: fir1(50, [30/500 200/500], hanning(51))
__constant__ double filter_b_[51];
const double filter_b[51] = {4.63918326210317e-05,6.01487136220455e-06,0.000153701779613429,0.00118607963150578,0.00225266327355902,0.00117569604807060,-0.00148147289905911,-0.00129435957227171,0.00325847465547368,0.00529818550591183,-0.00249277177208264,-0.0143174375226661,-0.0149918024814627,-0.00271416877295211,0.00133963485720405,-0.0204930171505778,-0.0503049882244718,-0.0494438908010180,
	-0.0143915554396815,0.00214300336155007,-0.0469199353599694,-0.122636972326032,-0.115283376273340,0.0344336539584967,0.242069171072009,0.339817488087583,0.242069171072009,0.0344336539584967,-0.115283376273340,-0.122636972326032,-0.0469199353599694,0.00214300336155007,-0.0143915554396815,-0.0494438908010180,-0.0503049882244718,-0.0204930171505778,0.00133963485720405,-0.00271416877295211,
	-0.0149918024814627,-0.0143174375226661,-0.00249277177208264,0.00529818550591183,0.00325847465547368,-0.00129435957227171,-0.00148147289905911,0.00117569604807060,0.00225266327355902,0.00118607963150578,0.000153701779613429,6.01487136220455e-06,4.63918326210317e-05};



void inputFiles()
{
	int i, numAbl, nodeIdx, tempIn;
	FILE *fin;

	printf("Importing V_out.dat...\n");
	fin = fopen("V_out.dat", "rb");
	fscanf(fin, "%d", &nde);
	fclose(fin);

	CL = (float*)malloc(sizeof(float)*nde);
	num1 = (int*)malloc(sizeof(int)*(nde+1));
	num2 = (int*)malloc(sizeof(int)*(nde+1));
	ifScar = (bool*)malloc(sizeof(bool)*nde);

	for(i=0; i<nde; i++) ifScar[i] = false;

	printf("Importing electrode1_num.txt...\n");
	fin = fopen("electrode1_num.txt", "rb");
	for(i=1; i<=nde; i++) fscanf(fin, "%d", &num1[i]);
	num1[0] = 0;
	for(i=1; i<=nde; i++) num1[i] = num1[i] + num1[i-1];
	fclose(fin);

	printf("Importing electrode2_num.txt...\n");
	fin = fopen("electrode2_num.txt", "rb");
	for(i=1; i<=nde; i++) fscanf(fin, "%d", &num2[i]);
	num2[0] = 0;
	for(i=1; i<=nde; i++) num2[i] = num2[i] + num2[i-1];
	fclose(fin);

	node1 = (int*)malloc(sizeof(int)*num1[nde]);
	node2 = (int*)malloc(sizeof(int)*num2[nde]);

	printf("Importing electrode1_ele.txt...\n");
	fin = fopen("electrode1_ele.txt", "rb");
	for(i=0; i<num1[nde]; i++)
	{
		fscanf(fin, "%d", &node1[i]);
		node1[i]--;
	}
	fclose(fin);

	printf("Importing electrode2_ele.txt...\n");
	fin = fopen("electrode2_ele.txt", "rb");
	for(i=0; i<num2[nde]; i++)
	{
		fscanf(fin, "%d", &node2[i]);
		node2[i]--;
	}
	fclose(fin);

	printf("Importing cuda_ifscar.dat...\n");
	fin = fopen("cuda_ifscar.dat", "rb");
	if(fin != NULL)
	{
		for(i=0; i<nde; i++)
		{
			fscanf(fin, "%d", &tempIn);
			ifScar[i] = (bool)tempIn;
		}
		fclose(fin);
	}
	else printf("    (No such file...)\n");

	printf("Importing AblationNode.dat...\n");
	fin = fopen("AblationNode.dat", "rb");
	if(fin != NULL)
	{
		fscanf(fin, "%d", &numAbl);
		for(i=0; i<numAbl; i++)
		{
			fscanf(fin, "%d", &nodeIdx);
			if(nodeIdx <= nde) ifScar[nodeIdx-1] = true;
		}
		fclose(fin);
	}
	else printf("    (No such file...)\n");
}



__global__ void noiseNodeKernel(bool *noiseNode_, const bool *ifScar_, const int *num1_, const int *num2_, const int *node1_, const int *node2_, const int N)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int bdim = blockDim.x;
	const unsigned int gdim = gridDim.x;

	int step = bdim*gdim;
	int id, i, i1, i2;

	for(id = bid*bdim + tid; id<N; id+=step)
	{
		if(noiseNode_[id]) continue;

		// Electrode 1
		i1 = num1_[id]; i2 = num1_[id+1];
		for(i=i1; i<i2; i++)
		{
			if(ifScar_[node1_[i]]) noiseNode_[id] = true;
		}

		// Electrode 2
		i1 = num2_[id]; i2 = num2_[id+1];
		for(i=i1; i<i2; i++)
		{
			if(ifScar_[node2_[i]]) noiseNode_[id] = true;
		}
	}
}



__global__ void egmKernel(float *egm_, double *vm_, const int *num1_, const int *num2_, const int *node1_, const int *node2_, const int N)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int bdim = blockDim.x;
	const unsigned int gdim = gridDim.x;

	int step = bdim*gdim;
	int id, i, i1, i2;
	double volt1, volt2;

	for(id = bid*bdim + tid; id<N; id+=step)
	{
		volt1 = 0.0;
		volt2 = 0.0;

		// Electrode 1
		i1 = num1_[id]; i2 = num1_[id+1];
		for(i=i1; i<i2; i++) volt1 += vm_[node1_[i]];

		if(i1 == i2) volt1 = 0.0;
		else volt1 = volt1 / (double)(i2-i1);

		// Electrode 2
		i1 = num2_[id]; i2 = num2_[id+1];
		for(i=i1; i<i2; i++) volt2 += vm_[node2_[i]];

		if(i1 == i2) volt2 = 0.0;
		else volt2 = volt2 / (double)(i2-i1);

		egm_[id] = (float)((int)((volt1 - volt2)*1000.0+0.5))/1000.0;
	}
}



__global__ void filterKernel(float *egm_raw_, float *egm_filtered_, const int M)
{
	int id = blockIdx.x*blockDim.x + threadIdx.x;
	int i, imax;
	double temp = 0.0;

	if(id < M)
	{
		if(id < 50) imax = id;
		else imax = 50;

		for(i=0; i<=imax; i++) temp += filter_b_[i]*(double)egm_raw_[id-i];

		// egm_filtered_[id] = (float)((int)(temp*1000.0+0.5))/1000.0; // filter
		egm_filtered_[M-1-id] = (float)((int)(temp*1000.0+0.5))/1000.0; // filtfilt
	}
}



float calcCL(float *egm_)
{
	float diff[8000], min;
	int i, j, sg[8000], sgCnt = 0, FI[8000], FICnt = 0;
	int REF_FI, CLCnt = 0;

	for(i=0; i<sampleNum-1; i++) diff[i] = egm_[i+1]-egm_[i];
	diff[sampleNum-1] = 0.0;

	for(i=1; i<sampleNum; i++)
	{
		if((diff[i]<0 && diff[i-1]>0) || (diff[i]>0 && diff[i-1]<0)) sg[sgCnt++] = i;
	}

	for(i=0; i<sgCnt-1; i++)
	{
		if(diff[sg[i]-1]>0 && diff[sg[i+1]-1]<0 && sg[i+1]-sg[i]<WIDTH && egm_[sg[i]]-egm_[sg[i+1]]>=PPSENSITIVITY && sg[i]>0)
		{
			min = 10000.0;
			for(j=sg[i]; j<=sg[i+1]; j++)
			{
				if(diff[j] < min) min = diff[j];
			}

			for(j=sg[i]; j<=sg[i+1]; j++)
			{
				if(fabs(diff[j]-min) <= 0.000001)
				{
					FI[FICnt++] = j;
					break;
				}
			}

			if(j < REFRACTORY-1) FICnt--;
		}
	}

	REF_FI = FI[0];
	for(i=1; i<FICnt; i++)
	{
		if(FI[i]-REF_FI >= REFRACTORY)
		{
			REF_FI = FI[i];
			CLCnt++;
		}
	}

	if(CLCnt == 0) return 10000.0;
	else return (float)(REF_FI-FI[0])/(float)CLCnt;
}



void printCL()
{
	FILE *fin, *CLdata, *CLmap;
	int i, node1, node2, node3;
	float x, y, z, temp;

	fin = fopen("V_out.dat", "rb");
	CLdata = fopen("CFAE_CL_data.txt", "wb");
	CLmap = fopen("CFAE_CL_map.plt", "wb");

	fscanf(fin, "%d", &nde);
	fscanf(fin, "%d", &nel);

	fprintf(CLmap, "VARIABLES = \"X\", \"Y\", \"Z\", \"CL\"\n");
	fprintf(CLmap, "ZONE F=FEPOINT, ET=triangle, N=%d , E=%d\n", nde, nel);

	for(i=0; i<nde; i++)
	{
		fprintf(CLdata, "%.1f\n", CL[i]);

		fscanf(fin, "%f %f %f %f", &x, &y, &z, &temp);
		fprintf(CLmap, "%f %f %f %.1f\n", x, y, z, CL[i]);
	}

	for(i=0; i<nel; i++)
	{
		fscanf(fin, "%d %d %d", &node1, &node2, &node3);
		fprintf(CLmap, "%d %d %d\n", node1, node2, node3);
	}


	fclose(fin);
	fclose(CLdata);
	fclose(CLmap);

	printf("Output files:\n     CFAE_CL_data.txt\n     CFAE_CL_map.plt\n");
}



int main()
{
	int filenum, i, j;
	FILE *datafile;
	char datafilename[100];
	float *egmTmp;

	printf("<CFAE 3D 2.2>   (2014.12.18)\n\n");

	printf("Input files:\n     V_out.dat\n     electrode1_num.txt\n     electrode2_num.txt\n     electrode1_ele.txt\n     electrode2_ele.txt\n     cuda_ifscar.dat (if existed)\n     AblationNode.dat (if existed)\n     VoltBinary_#.txt\n\n");

	printf("Refractory (default: 40 ms): "); scanf("%f", &REFRACTORY);
	printf("Width (default: 15 ms): "); scanf("%f", &WIDTH);
	printf("Peak to peak sensitivity (default: 2.0 mV): "); scanf("%f", &PPSENSITIVITY);

	printf("Start # (eg. 1): ");
	scanf("%d", &start);
	printf("Etart # (eg. 6000): ");
	scanf("%d", &end);
	sampleNum = end-start+1;
	printf("\n");

	inputFiles();

	vm = (double*)malloc(sizeof(double)*nde);
	egm = (float**)malloc(sizeof(float*)*nde);
	for(i=0; i<nde; i++) egm[i] = (float*)malloc(sizeof(float)*sampleNum);
	egmTmp = (float*)malloc(sizeof(float)*nde);



	// -------------------- CUDA setting --------------------
	int *num1_d, *num2_d, *node1_d, *node2_d;
	double *vm_d;
	float *egm_d;
	bool *ifScar_d, *noiseNode;

	cudaMalloc((void**)&num1_d, sizeof(int)*(nde+1));
	cudaMalloc((void**)&num2_d, sizeof(int)*(nde+1));
	cudaMalloc((void**)&node1_d, sizeof(int)*num1[nde]);
	cudaMalloc((void**)&node2_d, sizeof(int)*num2[nde]);
	cudaMalloc((void**)&vm_d, sizeof(double)*nde);
	cudaMalloc((void**)&egm_d, sizeof(float)*nde);
	cudaMalloc((void**)&ifScar_d, sizeof(bool)*nde);
	cudaMalloc((void**)&noiseNode, sizeof(bool)*nde);

	cudaMemcpy(num1_d, num1, sizeof(int)*(nde+1), cudaMemcpyHostToDevice);
	cudaMemcpy(num2_d, num2, sizeof(int)*(nde+1), cudaMemcpyHostToDevice);
	cudaMemcpy(node1_d, node1, sizeof(int)*num1[nde], cudaMemcpyHostToDevice);
	cudaMemcpy(node2_d, node2, sizeof(int)*num2[nde], cudaMemcpyHostToDevice);
	cudaMemcpy(ifScar_d, ifScar, sizeof(bool)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(noiseNode, ifScar, sizeof(bool)*nde, cudaMemcpyHostToDevice);

	free(num1);
	free(num2);
	free(node1);
	free(node2);



	// -------------------- Indentify Noise Node --------------------
	noiseNodeKernel<<<256, 128>>>(noiseNode, ifScar_d, num1_d, num2_d, node1_d, node2_d, nde);
	cudaMemcpy(ifScar, noiseNode, sizeof(bool)*nde, cudaMemcpyDeviceToHost);
	cudaFree(ifScar_d);
	cudaFree(noiseNode);



	// -------------------- Calculating EGM --------------------
	printf("\nCalculating EGM...\n");
	for(filenum=start, j=0; filenum<=end; filenum++, j++)
	{
		sprintf(datafilename, "VoltBinary_%d.txt", filenum);
		datafile = fopen(datafilename, "rb");
		fread(vm, sizeof(double), nde, datafile);
		cudaMemcpy(vm_d, vm, sizeof(double)*nde, cudaMemcpyHostToDevice);
		fclose(datafile);
		printf("     %s\n", datafilename);

		egmKernel<<<256, 128>>>(egm_d, vm_d, num1_d, num2_d, node1_d, node2_d, nde);

		cudaMemcpy(egmTmp, egm_d, sizeof(float)*nde, cudaMemcpyDeviceToHost);
		for(i=0; i<nde; i++) egm[i][j] = egmTmp[i];
	}

	free(vm);
	free(egmTmp);
	cudaFree(vm_d);
	cudaFree(egm_d);
	cudaFree(num1_d);
	cudaFree(num2_d);
	cudaFree(node1_d);
	cudaFree(node2_d);



	// -------------------- FIR filter (all-zero filter) --------------------
	//printf("\nFIR filter (all-zero filter)...\n");
	printf("\n30-200Hz FIR forward and reverse filter...\n");

	float *egm_raw, *egm_filtered;
	cudaMalloc((void**)&egm_raw, sizeof(float)*sampleNum);
	cudaMalloc((void**)&egm_filtered, sizeof(float)*sampleNum);

	cudaMemcpyToSymbol(filter_b_, filter_b, sizeof(double)*51);

	for(i=0; i<nde; i++)
	{
		if(i%10000 == 0) printf("     %d / %d\n", i, nde);
		if(ifScar[i]) continue;

		/*
		cudaMemcpy(egm_raw, egm[i], sizeof(float)*sampleNum, cudaMemcpyHostToDevice);
		filterKernel<<<48, 128>>>(egm_raw, egm_filtered, sampleNum);
		cudaMemcpy(egm[i], egm_filtered, sizeof(float)*sampleNum, cudaMemcpyDeviceToHost);
		*/
		cudaMemcpy(egm_raw, egm[i], sizeof(float)*sampleNum, cudaMemcpyHostToDevice);
		filterKernel<<<48, 128>>>(egm_raw, egm_filtered, sampleNum);
		filterKernel<<<48, 128>>>(egm_filtered, egm_raw, sampleNum);
		cudaMemcpy(egm[i], egm_raw, sizeof(float)*sampleNum, cudaMemcpyDeviceToHost);
	}

	cudaFree(egm_raw);
	cudaFree(egm_filtered);



	// -------------------- Calculating CFAE-CL --------------------
	printf("\nCalculating CFAE-CL...\n");

	for(i=0; i<nde; i++)
	{
		CL[i] = ifScar[i] ? 10000.0 : calcCL(egm[i]);
		if(i%10000 == 0) printf("     %d / %d\n", i, nde);
	}

	free(ifScar);



	// -------------------- Print --------------------
	printf("\nPrint...\n\n");
	printCL();


	printf("\n-------------------- END --------------------\n\n");

	int nodenumber;
	while(1)
	{
		printf("Type node number to print EGM (index base 1): ");
		scanf("%d", &nodenumber);
		nodenumber--;

		sprintf(datafilename, "EGM_%d.txt", nodenumber+1);
		datafile = fopen(datafilename, "wb");
		for(i=0; i<sampleNum; i++) fprintf(datafile, "%.2f\n", egm[nodenumber][i]);
		fclose(datafile);

		printf("\n");
	}

	return 0;
}