/*
Bipolar EGM, CFAE-CL Calculator for 2D Simulation
Jun-Seop Song 2014.11.24
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


#define M 600 // M by M tissue
#define sampleNum 6000 // CFAE duration

#define REFRACTORY 40 // ms
#define WIDTH 15 // ms
#define PPSENSITIVITY 4.0 // mV


#define ERow1 3.5
#define ECol1 2.0
#define ERow2 1.5
#define ECol2 2.0
#define STEP 1

#define index(a,b) a*M+b

int ROW, COL, vmSize, egmSize;
int start, end;
float *vm, **egm;


__constant__ double filter_b_[51];
const double filter_b[51] = {4.63918326210317e-05,6.01487136220455e-06,0.000153701779613429,0.00118607963150578,0.00225266327355902,0.00117569604807060,-0.00148147289905911,-0.00129435957227171,0.00325847465547368,0.00529818550591183,-0.00249277177208264,-0.0143174375226661,-0.0149918024814627,-0.00271416877295211,0.00133963485720405,-0.0204930171505778,-0.0503049882244718,-0.0494438908010180,
	-0.0143915554396815,0.00214300336155007,-0.0469199353599694,-0.122636972326032,-0.115283376273340,0.0344336539584967,0.242069171072009,0.339817488087583,0.242069171072009,0.0344336539584967,-0.115283376273340,-0.122636972326032,-0.0469199353599694,0.00214300336155007,-0.0143915554396815,-0.0494438908010180,-0.0503049882244718,-0.0204930171505778,0.00133963485720405,-0.00271416877295211,
	-0.0149918024814627,-0.0143174375226661,-0.00249277177208264,0.00529818550591183,0.00325847465547368,-0.00129435957227171,-0.00148147289905911,0.00117569604807060,0.00225266327355902,0.00118607963150578,0.000153701779613429,6.01487136220455e-06,4.63918326210317e-05};



__global__ void egmKernel(float *egm_, float *vm_, const int colSize, const int N)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int bdim = blockDim.x;
	const unsigned int gdim = gridDim.x;

	int step = bdim*gdim;
	int id, rowId, colId, i, j;
	float volt1, volt2;

	for(id = bid*bdim + tid; id<N; id+=step)
	{
		rowId = (int)(id/colSize);
		colId = id%colSize;

		volt1 = 0.0;
		volt2 = 0.0;

		// Electrode 1
		for(i = colId; i <= colId+(ECol1*4-1); i++)
		{
			for(j = rowId; j <= rowId+(ERow1*4-1); j++)
			{
				volt1 += vm_[index(i, j)];
			}
		}
		volt1 = volt1 / ((ERow1*4)*(ECol1*4));

		// Electrode 2
		for(i = colId; i <= colId+(ECol2*4-1); i++)
		{
			for(j = rowId+(ERow1*4+STEP*4); j <= rowId+(ERow1*4+STEP*4)+(ERow2*4-1); j++)
			{
				volt2 += vm_[index(i, j)];
			}
		}
		volt2 = volt2 / ((ERow2*4)*(ECol2*4));

		// EGM
		egm_[id] = volt1 - volt2;
	}
}



__global__ void filterKernel(float *egm_raw_, float *egm_filtered_, const int NN)
{
	int id = blockIdx.x*blockDim.x + threadIdx.x;
	int i, imax;
	double temp = 0.0;

	if(id < NN)
	{
		if(id < 50) imax = id;
		else imax = 50;

		for(i=0; i<=imax; i++) temp += filter_b_[i]*(double)egm_raw_[id-i];

		// egm_filtered_[id] = (float)((int)(temp*1000.0+0.5))/1000.0; // filter
		egm_filtered_[NN-1-id] = (float)((int)(temp*1000.0+0.5))/1000.0; // filtfilt
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
				if(fabs(diff[j]-min) < 0.000001)
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



int main()
{
	int filenum, i, j;
	FILE *datafile;
	char datafilename[100];
	float *egmTmp;

	printf("<Bipolar EGM, CFAE-CL Calculator for 2D Simulation>\n");
	printf("JSSONG 2014.11.24\n\n");

	printf("Start time (ms): ");
	scanf("%d", &start);


	ROW = M-(4*ERow1+4*STEP+4*ERow2)+1;
	if(ECol1 > ECol2) COL = M-4*ECol1+1;
	else COL = M-4*ECol2+1;
	vmSize = M*M;
	egmSize = ROW*COL;
	end = sampleNum+start-1;
	printf("\nstart: %d   end: %d\n", start, end);


	printf("\nAllocating memory...\n");
	vm = (float*)malloc(sizeof(float)*vmSize);
	egm = (float**)malloc(sizeof(float*)*egmSize);
	for(i=0; i<egmSize; i++) egm[i] = (float*)malloc(sizeof(float)*sampleNum);
	egmTmp = (float*)malloc(sizeof(float)*egmSize);


	// -------------------- CUDA setting --------------------
	float *vm_d, *egm_d;

	cudaMalloc((void**)&vm_d, sizeof(float)*vmSize);
	cudaMalloc((void**)&egm_d, sizeof(float)*egmSize);


	// -------------------- Calculating EGM --------------------
	printf("\nCalculating EGM...\n");
	for(filenum=start, j=0; filenum<=end; filenum++, j++)
	{
		sprintf(datafilename, "vm%d.txt", filenum);
		datafile = fopen(datafilename, "rb");
		for(i=0; i<vmSize; i++) fscanf(datafile, "%f", &vm[i]);
		cudaMemcpy(vm_d, vm, sizeof(float)*vmSize, cudaMemcpyHostToDevice);
		fclose(datafile);
		printf("     %s\n", datafilename);

		egmKernel<<<256, 128>>>(egm_d, vm_d, COL, egmSize);

		cudaMemcpy(egmTmp, egm_d, sizeof(float)*egmSize, cudaMemcpyDeviceToHost);

		for(i=0; i<egmSize; i++) egm[i][j] = egmTmp[i];
	}

	free(vm);
	free(egmTmp);
	cudaFree(vm_d);
	cudaFree(egm_d);


	// -------------------- FIR filter (all-zero filter) --------------------
	printf("\nFIR filter (all-zero filter)...\n");

	float *egm_raw, *egm_filtered;
	cudaMalloc((void**)&egm_raw, sizeof(float)*sampleNum);
	cudaMalloc((void**)&egm_filtered, sizeof(float)*sampleNum);

	cudaMemcpyToSymbol(filter_b_, filter_b, sizeof(double)*51);

	for(i=0; i<egmSize; i++)
	{
		/*
		cudaMemcpy(egm_raw, egm[i], sizeof(float)*sampleNum, cudaMemcpyHostToDevice);
		filterKernel<<<48, 128>>>(egm_raw, egm_filtered, sampleNum);
		cudaMemcpy(egm[i], egm_filtered, sizeof(float)*sampleNum, cudaMemcpyDeviceToHost);
		*/
		cudaMemcpy(egm_raw, egm[i], sizeof(float)*sampleNum, cudaMemcpyHostToDevice);
		filterKernel<<<48, 128>>>(egm_raw, egm_filtered, sampleNum);
		filterKernel<<<48, 128>>>(egm_filtered, egm_raw, sampleNum);
		cudaMemcpy(egm[i], egm_raw, sizeof(float)*sampleNum, cudaMemcpyDeviceToHost);

		if(i%10000 == 0) printf("     %d / %d\n", i, egmSize);
	}

	cudaFree(egm_raw);
	cudaFree(egm_filtered);


	// -------------------- Calculating CFAE-CL --------------------
	printf("\nCalculating CFAE-CL...\n");

	sprintf(datafilename, "CFAE_CL_%d_%d.txt", start, end);
	datafile = fopen(datafilename, "wb");

	for(i=0; i<egmSize; i++)
	{
		fprintf(datafile, "%.1f\n", calcCL(egm[i]));
		if(i%10000 == 0) printf("     %d / %d\n", i, egmSize);
	}

	fclose(datafile);


	// ---------------- Print EGM ----------------
	printf("\n------------------------------------------------\n");
	int nodeNumber;
	while(1)
	{
		printf("\nType node number to export EGM: ");
		scanf("%d", &nodeNumber);

		sprintf(datafilename, "egm%d.txt", nodeNumber);
		datafile = fopen(datafilename, "wb");
		for(i=0; i<sampleNum; i++) fprintf(datafile, "%.3f\n", egm[nodeNumber][i]);
		fclose(datafile);
	}


	return 0;
}