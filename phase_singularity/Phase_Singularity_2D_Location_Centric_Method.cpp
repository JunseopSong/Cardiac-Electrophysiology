/*
Phase Singularity Detection (from 2D simulation)
+++ Location-Centric Method +++
(developed by Jun-Seop Song & Young-Seon Lee)

Reference:
Lee, Y. S., Song, J. S., Hwang, M., Lim, B., Joung, B., & Pak, H. N. (2016). A new efficient method for detecting phase singularity in cardiac fibrillation. PLoS one, 11(12), e0167567.
https://doi.org/10.1371/journal.pone.0167567
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define tau 30 // ms
#define nde 512*512 // size of tissue
#define startTime 1  // ms
#define timeNum 6000  // ms


int nel;
float **vm;
bool *ifPhaseSingular;
int *PStime;


void printPS()
{
	FILE *PSdata = fopen("PS_data.txt", "w");

	for (int i = 0; i < nde; ++i)
	{
		if (ifPhaseSingular[i]) fprintf(PSdata, "1\n");
		else fprintf(PSdata, "0\n");
	}

	fclose(PSdata);
}



int main()
{
	int i, j;
	FILE *datafile;
	char fileName[100];
	float *vm_smooth, F_mean, theta, theta_pre;
	double *vm_input;
	clock_t t1, t2, t3;

	printf("<Phase Singularity 2D>\n\n");

	vm = (float**)malloc(sizeof(float*)*nde);
	vm_input = (double*)malloc(sizeof(double)*nde);
	for (i = 0; i < nde; ++i) vm[i] = (float*)malloc(sizeof(float)*(timeNum + tau));
	ifPhaseSingular = (bool*)malloc(sizeof(bool)*nde);
	vm_smooth = (float*)malloc(sizeof(float)*(timeNum + tau));
	PStime = (int*)malloc(sizeof(int)*nde);

	for (i = 0; i < nde; ++i)
	{
		ifPhaseSingular[i] = false;
		PStime[i] = -100;
	}

	t1 = clock();

	// -------------------- Input 2D simulation data --------------------
	printf("\nImporting voltage data...\n");
	for (j = 0; j < timeNum + tau; ++j)
	{
		sprintf(fileName, "vm%d.txt", (j + startTime));
		datafile = fopen(fileName, "rb");
		fread(vm_input, sizeof(double), nde, datafile);
		for (i = 0; i < nde; ++i) vm[i][j] = (float)vm_input[i];
		fclose(datafile);
		printf("%d\n", j);
	}
	free(vm_input);

	t2 = clock();

	// -------------------- Calculating PS --------------------
	printf("\nCalculating phase singularity...\n");
	for (i = 0; i < nde; ++i)
	{
		for (j = 1; j < timeNum + tau - 1; ++j) vm_smooth[j] = (vm[i][j - 1] + vm[i][j] + vm[i][j + 1]) / 3.0;
		vm_smooth[0] = vm[i][0];
		vm_smooth[timeNum + tau - 1] = vm[i][timeNum + tau - 1];

		F_mean = 0.0;
		for (j = tau; j < timeNum + tau; ++j) F_mean += vm_smooth[j];
		F_mean = F_mean / (float)timeNum;

		theta_pre = (atan2(vm_smooth[tau] - F_mean, vm_smooth[0] - F_mean)*1000.0 + 0.5) / 1000.0;
		for (j = 1; j < timeNum; ++j)
		{
			theta = (atan2(vm_smooth[j + tau] - F_mean, vm_smooth[j] - F_mean)*1000.0 + 0.5) / 1000.0;

			if (theta - theta_pre < -3.14159) // PS detection criteria
			{
				ifPhaseSingular[i] = true;
				PStime[i] = j;
			}

			theta_pre = theta;
		}

		if (i % 10000 == 0) printf("     %d / %d\n", i, nde);
	}

	t3 = clock();

	printPS();

	printf("\n--------------------------------\n");
	printf("Data importing time: %.3f ms\n", (float)(t2 - t1));
	printf("Calculating time: %.3f ms\n", (float)(t3 - t2));
	printf("Total time: %.3f ms\n", (float)(t3 - t1));
	printf("--------------------------------\n");

	return 0;
}