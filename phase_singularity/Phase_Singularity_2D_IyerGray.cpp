/*
Phase Singularity Detection (from 2D simulation)
+++ Iyer-Gray Method +++
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define tau 30  // ms
#define N 512  // size of tissue
#define startTime 1  // ms
#define timeDur 6000  // ms

#define id(x,y) ((x)*N+(y))
#define PI 3.14159265358979
#define EPS 0.000001

const int nde = N*N;
float **vm, **theta;
bool *ifPhaseSingular;
int *PStime;


void printPS()
{
	FILE *PSdata = fopen("PS_data.txt", "w");

	for (int i = 0; i < nde; ++i)
	{
		if (ifPhaseSingular[i]) fprintf(PSdata, "1\n");  // PS point
		else fprintf(PSdata, "0\n");  // Not PS point
	}

	fclose(PSdata);
}


bool mod(float a, float m)
{
	a = fabs(a);
	while (a > 0) a -= m;
	if (a<EPS && a>-EPS) return 1;
	a += m;
	if (a<EPS && a>-EPS) return 1;
	return 0;
}


int main()
{
	int i, j, k;
	FILE *datafile;
	char fileName[100];
	float *vm_smooth, F_mean, K11, K21, K31, K32, K33, K23, K13, K12, K_SUM;
	double *vm_input;
	clock_t t1, t2, t3;

	printf("<Phase Singularity 2D - Iyer-Gray>\n\n");

	vm = (float**)malloc(sizeof(float*)*nde);
	vm_input = (double*)malloc(sizeof(double)*nde);
	for (i = 0; i < nde; ++i) vm[i] = (float*)malloc(sizeof(float)*(timeDur + tau));
	theta = (float**)malloc(sizeof(float*)*nde);
	for (i = 0; i < nde; ++i) theta[i] = (float*)malloc(sizeof(float)*timeDur);
	ifPhaseSingular = (bool*)malloc(sizeof(bool)*nde);
	vm_smooth = (float*)malloc(sizeof(float)*(timeDur + tau));
	PStime = (int*)malloc(sizeof(int)*nde);

	for (i = 0; i < nde; ++i)
	{
		ifPhaseSingular[i] = false;
		PStime[i] = -100;
	}

	t1 = clock();

	// -------------------- Input --------------------
	printf("\nImporting voltage data...\n");
	for (j = 0; j < timeDur + tau; ++j)
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

	// -------------------- Calculating theta --------------------
	printf("\nCalculating phase...\n");
	for (i = 0; i < nde; ++i)
	{
		for (j = 1; j < timeDur + tau - 1; ++j) vm_smooth[j] = (vm[i][j - 1] + vm[i][j] + vm[i][j + 1]) / 3.0;
		vm_smooth[0] = vm[i][0];
		vm_smooth[timeDur + tau - 1] = vm[i][timeDur + tau - 1];

		F_mean = 0.0;
		for (j = tau; j < timeDur + tau; ++j) F_mean += vm_smooth[j];
		F_mean = F_mean / (float)timeDur;

		for (j = 0; j < timeDur; ++j)
		{
			theta[i][j] = (atan2(vm_smooth[j + tau] - F_mean, vm_smooth[j] - F_mean)*1000.0 + 0.5) / 1000.0;
		}

		if (i % 10000 == 0) printf("     %d / %d\n", i, nde);
	}

	// -------------------- Identifying PS --------------------
	printf("\nIdentifying PS...\n");
	for (k = 0; k < timeDur; ++k)
	{
		for (i = 2; i < N - 1; ++i)
		{
			for (j = 2; j < N - 1; ++j)
			{
				K11 = theta[id(j, i - 1)][k] - theta[id(j - 1, i - 1)][k];
				if (!(fabs(K11) <= PI) && (K11 > 0)) K11 -= 2 * PI;
				else if (!(fabs(K11) <= PI) && (K11 < 0)) K11 += 2 * PI;

				K21 = theta[id(j + 1, i - 1)][k] - theta[id(j, i - 1)][k];
				if (!(fabs(K21) <= PI) && (K21 > 0)) K21 -= 2 * PI;
				else if (!(fabs(K21) <= PI) && (K21 < 0)) K21 += 2 * PI;

				K31 = theta[id(j + 1, i)][k] - theta[id(j + 1, i - 1)][k];
				if (!(fabs(K31) <= PI) && (K31 > 0)) K31 -= 2 * PI;
				else if (!(fabs(K31) <= PI) && (K31 < 0)) K31 += 2 * PI;

				K32 = theta[id(j + 1, i + 1)][k] - theta[id(j + 1, i)][k];
				if (!(fabs(K32) <= PI) && (K32 > 0)) K32 -= 2 * PI;
				else if (!(fabs(K32) <= PI) && (K32 < 0)) K32 += 2 * PI;

				K33 = theta[id(j, i + 1)][k] - theta[id(j + 1, i + 1)][k];
				if (!(fabs(K33) <= PI) && (K33 > 0)) K33 -= 2 * PI;
				else if (!(fabs(K33) <= PI) && (K33 < 0)) K33 += 2 * PI;

				K23 = theta[id(j - 1, i + 1)][k] - theta[id(j, i + 1)][k];
				if (!(fabs(K23) <= PI) && (K23 > 0)) K23 -= 2 * PI;
				else if (!(fabs(K23) <= PI) && (K23 < 0)) K23 += 2 * PI;

				K13 = theta[id(j - 1, i)][k] - theta[id(j - 1, i + 1)][k];
				if (!(fabs(K13) <= PI) && (K13 > 0)) K13 -= 2 * PI;
				else if (!(fabs(K13) <= PI) && (K13 < 0)) K13 += 2 * PI;

				K12 = theta[id(j - 1, i - 1)][k] - theta[id(j - 1, i)][k];
				if (!(fabs(K12) <= PI) && (K12 > 0)) K12 -= 2 * PI;
				else if (!(fabs(K12) <= PI) && (K12 < 0)) K12 += 2 * PI;

				K_SUM = K11 + K21 + K31 + K32 + K33 + K23 + K13 + K12;
				ifPhaseSingular[id(i, j)] |= mod(K_SUM, PI)&(fabs(K_SUM) > EPS);
			}
		}

		if (k % 100 == 0) printf("     %d / %d\n", k, timeDur);
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