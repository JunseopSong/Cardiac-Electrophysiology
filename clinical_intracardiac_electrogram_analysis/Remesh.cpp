#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>

using namespace std;


int GetNumbers(char* str)
{
	char buffer[100];
	int cnt=0;
	while(*str){
		if(isdigit(*str))
			buffer[cnt++] = *str;
		str++;
	}

	return atoi(buffer);
}



void xml_to_ply()
{
	int nde, nel, i;
	char tmp[5000];
	double *xyz[3];
	int *ie[3];
	std::ifstream fin;

	fin.open("dif001.xml");

	for(i=0; i<29; i++) fin.getline(tmp, 1000);

	nde = GetNumbers(tmp);
	for(i=0; i<3; i++) xyz[i] = new double[nde];

	for(i=0; i<nde; i++) fin  >> xyz[0][i]  >>  xyz[1][i]  >>  xyz[2][i];

	for(i=0; i<nde+5; i++) fin.getline(tmp, 1000);

	nel = GetNumbers(tmp);
	for(i=0; i<3; i++) ie[i] = new int[nel];

	for(i=0; i<nel; i++) fin  >> ie[0][i]  >>  ie[1][i]  >>  ie[2][i];

	fin.close();


	FILE *fout = fopen("RawMesh.ply", "w");
	fprintf(fout, "ply\n");
	fprintf(fout, "format ascii 1.0\n");
	fprintf(fout, "element vertex %d\n", nde);
	fprintf(fout, "property double x\n");
	fprintf(fout, "property double y\n");
	fprintf(fout, "property double z\n");
	fprintf(fout, "element face %d\n", nel);
	fprintf(fout, "property list uchar int vertex_index\n");
	fprintf(fout, "end_header\n");
	for(i=0; i<nde; i++) fprintf(fout, "%lf %lf %lf\n", xyz[0][i], xyz[1][i], xyz[2][i]);
	for(i=0; i<nel; i++) fprintf(fout, "3 %d %d %d\n", ie[0][i]-1, ie[1][i]-1, ie[2][i]-1);
	fclose(fout);

	
	for(i=0; i<3; i++) free(xyz[i]);
	for(i=0; i<3; i++) free(ie[i]);
}


void extractNdeNel(int *nde, int *nel, char *str)
{
	char bufferNde[20], bufferNel[20];
	int bufferPt=0;

	// Skip string
	while(!isdigit(*str)) str++;

	// Detect nde
	while(isdigit(*str))
	{
		bufferNde[bufferPt++] = *str;
		str++;
	}
	bufferNde[bufferPt] = NULL;

	bufferPt = 0;

	// Skip string
	while(!isdigit(*str)) str++;

	// Detect nel
	while(*str && isdigit(*str))
	{
		bufferNel[bufferPt++] = *str;
		str++;
	}
	bufferNel[bufferPt] = NULL;

	*nde = atoi(bufferNde);
	*nel = atoi(bufferNel);
}


void outputFormatConvert()
{
	char tmp[5000];
	std::ifstream fin;
	int nde, nel, i;
	double *xyz[3];
	int *ie[3];

	fin.open("output.dat");

	for(i=0; i<7; i++) fin.getline(tmp, 1000);

	extractNdeNel(&nde, &nel, tmp);
	for(i=0; i<3; i++) xyz[i] = new double[nde];
	for(i=0; i<3; i++) ie[i] = new int[nel];

	for(i=0; i<5; i++) fin.getline(tmp, 1000);

	for(i=0; i<nde; i++) fin  >> xyz[0][i]  >>  xyz[1][i]  >>  xyz[2][i];
	for(i=0; i<nel; i++) fin  >> ie[0][i]  >>  ie[1][i]  >>  ie[2][i];

	fin.close();


	FILE *fout = fopen("output_GS.dat", "w");
	fprintf(fout, "%d %d\n", nde, nel);
	for(i=0; i<nde; i++) fprintf(fout, "%lf %lf %lf 0\n", xyz[0][i], xyz[1][i], xyz[2][i]);
	for(i=0; i<nel; i++) fprintf(fout, "%d %d %d\n", ie[0][i], ie[1][i], ie[2][i]);
	fclose(fout);


	for(i=0; i<3; i++) free(xyz[i]);
	for(i=0; i<3; i++) free(ie[i]);
}


int main()
{
	int type;

	printf("1: dif001.xml  ->  RawMesh.ply\n");
	printf("2: acvd plymcout.ply 500000 0 -d 0\n");
	printf("3: output.dat (ASCII)  ->  output_GS.dat (input file for KNU GUI)\n");

	scanf("%d", &type);

	switch(type)
	{
	case 1:
		xml_to_ply();
		break;

	case 2:
		system("acvd plymcout.ply 500000 0 -d 0");
		remove("simplification.ply");
		break;

	case 3:
		outputFormatConvert();
		break;
	}


	return 0;
}