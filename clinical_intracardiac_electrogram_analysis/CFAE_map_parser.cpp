#include <fstream>
#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#define MAX_DEG 20

using namespace std;

double xyz[5000000][5], CL[5000000];
int ie[5000000][8], nde, nel, realNde, realNel;
char tmp[5000];
char* str = "";

int *deg, *link[MAX_DEG], *degree, *linkNode;
bool *isLAnode;

int GetNumbers(char* str);



void makeLinkArray()
{
	int i, j, k, k1, k2, v1, v2;
	bool ifNewNode;
	double temp;

	deg = (int*)malloc(sizeof(int)*nde);
	for(i=0; i<MAX_DEG; i++) link[i] = (int*)malloc(sizeof(int)*nde);


	for(i=0; i<nde; i++) deg[i] = 0;


	// ------------------ Get linking information ------------------
	for(i=0; i<nel; i++)
	{
		// Scan 1-1
		for(k1=0; k1<3; k1++)
		{
			k2 = (k1+1)%3;
			v1 = ie[i][k1]-1;
			v2 = ie[i][k2]-1;

			ifNewNode = true;
			for(j=0; j<deg[v1]; j++)
			{
				if(link[j][v1] == v2)
				{
					ifNewNode = false;
					break;
				}
			}

			if(ifNewNode)
			{
				link[deg[v1]][v1] = v2;
				deg[v1]++;
			}
		}

		// Scan 1-2
		for(k1=0; k1<3; k1++)
		{
			k2 = (k1+2)%3;
			v1 = ie[i][k1]-1;
			v2 = ie[i][k2]-1;

			ifNewNode = true;
			for(j=0; j<deg[v1]; j++)
			{
				if(link[j][v1] == v2)
				{
					ifNewNode = false;
					break;
				}
			}

			if(ifNewNode)
			{
				link[deg[v1]][v1] = v2;
				deg[v1]++;
			}
		}
	}



	// ------------------ Transform 2D array into 1D array ------------------
	degree = (int*)malloc(sizeof(int)*(nde+1));
	degree[0] = 0;
	for(i=1; i<=nde; i++) degree[i] = degree[i-1] + deg[i-1];


	linkNode = (int*)malloc(sizeof(int)*degree[nde]);
	for(i=0; i<nde; i++)
	{
		for(j=degree[i]; j<degree[i+1]; j++)
		{
			linkNode[j] = link[j-degree[i]][i];
		}
	}
}


void detectLA(int nodeId)
{
	int i;

	isLAnode[nodeId] = true;

	for(i=degree[nodeId]; i<degree[nodeId+1]; i++){
		if(!isLAnode[linkNode[i]]) detectLA(linkNode[i]);
	}
}


int main()
{
	double maxZ;
	int i, LAnode;

	char tmp1[1000];

	std::ifstream fin;
	fin.open("DxLandmarkGeo.xml");

	for(i=0; i<20; i++) fin.getline(tmp, 1000);

	nde = GetNumbers(tmp);

	maxZ = -1000000;
	for(i=0; i<nde; i++){
		fin  >> xyz[i][0]  >>  xyz[i][1]  >>  xyz[i][2];

		if(xyz[i][2] > maxZ){
			maxZ = xyz[i][2];
			LAnode = i;
		}
	}

	for(i=0; i<4; i++) fin.getline(tmp, 1000);

	for(i=0; i<nde; i++){
		fin  >> CL[i];
	}

	for(i=0; i<(3*nde+15); i++) fin.getline(tmp, 1000);

	nel = GetNumbers(tmp);

	for(i=0; i<nel; i++){
		fin  >> ie[i][0]  >>  ie[i][1]  >>  ie[i][2];
	}

	fin.close();



	// Remove esophagus
	makeLinkArray();

	isLAnode = new bool[nde];
	for(i=0; i<nde; i++) isLAnode[i] = false;

	detectLA(LAnode);


	// Count nde, nel without esophagus
	realNde = 0;
	for(i=0; i<nde; i++){
		if(isLAnode[i]) realNde++;
	}


	if(realNde < nde/2) // fix LA detection error
	{
		for(i=0; i<nde; i++) isLAnode[i] = !isLAnode[i];

		realNde = 0;
		for(i=0; i<nde; i++){
			if(isLAnode[i]) realNde++;
		}
	}


	realNel = 0;
	for(i=0; i<nel; i++){
		if(isLAnode[ie[i][0]-1] & isLAnode[ie[i][1]-1] & isLAnode[ie[i][2]-1]) realNel++;
	}



	// Clinical_CFAE_Data.plt
	std::ofstream fout;
	fout.open("Clinical_CFAE_Data.plt");

	fout << "VARIABLES = \"X\", \"Y\", \"Z\", \"CL\"" << std::endl;
	fout << "ZONE F=FEPOINT, ET=triangle, N=" << realNde << " , E=" << realNel <<std::endl;

	for(i=0; i<nde; i++){
		if(isLAnode[i]) fout << xyz[i][0] << " " << xyz[i][1] << " " << xyz[i][2] << " " << CL[i] << std::endl;
	}

	for(i=0; i<nel; i++){
		if(isLAnode[ie[i][0]-1] & isLAnode[ie[i][1]-1] & isLAnode[ie[i][2]-1]) fout << ie[i][0] << " " << ie[i][1] << " " << ie[i][2] << std::endl;
	}

	fout.close();


	// CFAE_CL_data.txt
	fout.open("CFAE_CL_data.txt");
	for(i=0; i<nde; i++){
		if(isLAnode[i]) fout << CL[i] << std::endl;
	}
	fout.close();


	// GSR file
	fout.open("Clinical_GSR(Rename_This_File_To_Input_dot_dat).dat");

	fout << realNde << " " << realNel <<std::endl;

	for(i=0; i<nde; i++){
		if(isLAnode[i]) fout << xyz[i][0] << " " << xyz[i][1] << " " << xyz[i][2] << std::endl;
	}

	for(i=0; i<nel; i++){
		if(isLAnode[ie[i][0]-1] & isLAnode[ie[i][1]-1] & isLAnode[ie[i][2]-1]) fout << ie[i][0] << " " << ie[i][1] << " " << ie[i][2] << " 0" << std::endl;
	}

	fout.close();


	/*
	// Vtk File
	fout.open("RawMesh.vtk");

	fout << "# vtk DataFile Version 1.0" << std::endl << "Atrium" << std::endl << "ASCII" << std::endl << std::endl;

	fout << "DATASET POLYDATA" << std::endl << "POINTS " << realNde << " double" <<std::endl;

	for(i=0; i<nde; i++){
	if(isLAnode[i]) fout << xyz[i][0] << " " << xyz[i][1] << " " << xyz[i][2] << std::endl;
	}

	fout << std::endl << "POLYGONS " << realNel << " " << realNel*4 << std::endl;

	for(i=0; i<nel; i++){
	if(isLAnode[ie[i][0]-1] & isLAnode[ie[i][1]-1] & isLAnode[ie[i][2]-1]) fout << "3 " << ie[i][0]-1 << " " << ie[i][1]-1 << " " << ie[i][2]-1 << std::endl;
	}
	fout.close();


	// Ply File
	fout.open("RawMesh.ply");

	fout << "ply" << std::endl << "format ascii 1.0" << std::endl;
	fout << "element vertex " << realNde << std::endl;
	fout << "property double x" << std::endl << "property double y" << std::endl << "property double z" << std::endl;
	fout << "element face " << realNel << std::endl;
	fout << "property list uchar int vertex_index" << std::endl << "end_header" << std::endl;

	for(i=0; i<nde; i++){
	if(isLAnode[i]) fout << xyz[i][0] << " " << xyz[i][1] << " " << xyz[i][2] << std::endl;
	}

	for(i=0; i<nel; i++){
	if(isLAnode[ie[i][0]-1] & isLAnode[ie[i][1]-1] & isLAnode[ie[i][2]-1]) fout << "3 " << ie[i][0]-1 << " " << ie[i][1]-1 << " " << ie[i][2]-1 << std::endl;
	}
	fout.close();
	*/


	return 0;
}


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

