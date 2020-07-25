/*
Phase Singularity Detection (from 3D simulation)
+++ Location-Centric Method +++
(developed by Jun-Seop Song & Young-Seon Lee)

Reference:
Lee, Y. S., Song, J. S., Hwang, M., Lim, B., Joung, B., & Pak, H. N. (2016). A new efficient method for detecting phase singularity in cardiac fibrillation. PLoS one, 11(12), e0167567.
https://doi.org/10.1371/journal.pone.0167567
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define tau 30 // ms


int nde, nel, start, end, timeNum;
float **vm;
bool *ifPhaseSingular;
int *PStime;



void printPS()
{
	FILE *fin, *PSdata, *PSmap;
	int i, node1, node2, node3;
	float x, y, z, temp;

	fin = fopen("V_out.dat", "rb");
	PSdata = fopen("PS_data.txt", "wb");
	PSmap = fopen("PS_map.plt", "wb");

	fscanf(fin, "%d", &nde);
	fscanf(fin, "%d", &nel);

	fprintf(PSmap, "VARIABLES = \"X\", \"Y\", \"Z\", \"Time\"\n");
	fprintf(PSmap, "ZONE F=FEPOINT, ET=triangle, N=%d , E=%d\n", nde, nel);

	for(i=0; i<nde; i++)
	{
		if(ifPhaseSingular[i]) fprintf(PSdata, "1\n");
		else fprintf(PSdata, "0\n");

		fscanf(fin, "%f %f %f %f", &x, &y, &z, &temp);
		fprintf(PSmap, "%f %f %f %d\n", x, y, z, PStime[i]);
		//if(ifPhaseSingular[i]) fprintf(PSmap, "%f %f %f 1\n", x, y, z);
		//else fprintf(PSmap, "%f %f %f 0\n", x, y, z);
	}

	for(i=0; i<nel; i++)
	{
		fscanf(fin, "%d %d %d", &node1, &node2, &node3);
		fprintf(PSmap, "%d %d %d\n", node1, node2, node3);
	}


	fclose(fin);
	fclose(PSdata);
	fclose(PSmap);
}



int main()
{
	int filenum, i, j;
	FILE *datafile;
	char datafilename[100];
	float *vm_smooth, F_mean, theta, theta_pre;
	double *vm_input;


	printf("<Phase Singularity 3D>   (2014.12.01)\n\n");

	printf("Start # (eg. 1): ");
	scanf("%d", &start);
	printf("Etart # (eg. 6000): ");
	scanf("%d", &end);
	timeNum = end-start+1;
	printf("\n");

	datafile = fopen("V_out.dat", "rb");
	fscanf(datafile, "%d", &nde);
	fclose(datafile);

	vm = (float**)malloc(sizeof(float*)*nde);
	vm_input = (double*)malloc(sizeof(double)*nde);
	for(i=0; i<nde; i++) vm[i] = (float*)malloc(sizeof(float)*(timeNum+tau));
	ifPhaseSingular = (bool*)malloc(sizeof(bool)*nde);
	vm_smooth = (float*)malloc(sizeof(float)*(timeNum+tau));
	PStime = (int*)malloc(sizeof(int)*nde);

	for(i=0; i<nde; i++)
	{
		ifPhaseSingular[i] = false;
		PStime[i] = -100;
	}


	// -------------------- Calculating Phase Singularity --------------------
	printf("\nImporting voltage data...\n");
	for(filenum=start, j=0; filenum<=end+tau; filenum++, j++)
	{
		sprintf(datafilename, "VoltBinary_%d.txt", filenum);
		datafile = fopen(datafilename, "rb");
		fread(vm_input, sizeof(double), nde, datafile);
		for(i=0; i<nde; i++) vm[i][j] = (float)vm_input[i];
		fclose(datafile);
		printf("     %s\n", datafilename);
	}
	free(vm_input);

	// -------------------- Calculating Phase Singularity --------------------
	printf("\nCalculating phase singularity...\n");
	for(i=0; i<nde; i++)
	{
		for(j=1; j<timeNum+tau-1; j++) vm_smooth[j] = (vm[i][j-1]+vm[i][j]+vm[i][j+1])/3.0;
		vm_smooth[0] = vm[i][0];
		vm_smooth[timeNum+tau-1] = vm[i][timeNum+tau-1];

		F_mean = 0.0;
		for(j=tau; j<timeNum+tau; j++) F_mean += vm_smooth[j];
		F_mean = F_mean/(float)timeNum;

		theta_pre = (atan2(vm_smooth[tau]-F_mean, vm_smooth[0]-F_mean)*1000.0 + 0.5)/1000.0;
		for(j=1; j<timeNum; j++)
		{
			theta = (atan2(vm_smooth[j+tau]-F_mean, vm_smooth[j]-F_mean)*1000.0 + 0.5)/1000.0;

			if(theta-theta_pre < -3.14159) // PS detection criteria
			{
				ifPhaseSingular[i] = true;
				PStime[i] = j;
			}

			theta_pre = theta;
		}

		if(i%10000 == 0) printf("     %d / %d\n", i, nde);
	}

	printPS();

	return 0;
}