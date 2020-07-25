#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#define MAXSTR 100000

char *strBuffer;


int getInt()
{
	char buffer[20];
	int cnt = 0;
	while (*strBuffer) {
		if (isdigit(*strBuffer))
			buffer[cnt++] = *strBuffer;
		else if (cnt != 0) break;
		strBuffer++;
	}
	buffer[cnt] = NULL;

	return atoi(buffer);
}


double getDouble()
{
	char buffer[20];
	int cnt = 0;
	while (*strBuffer) {
		if (isdigit(*strBuffer) || *strBuffer == '.' || *strBuffer == '-')
			buffer[cnt++] = *strBuffer;
		else if (cnt != 0) break;
		strBuffer++;
	}
	buffer[cnt] = NULL;

	return atof(buffer);
}


int main()
{
	FILE *fid, *fout;
	int fileId, i, numData;
	char fileName[100], fileIdString[2], foutName[100];

	strBuffer = new char[MAXSTR];


	for (fileId = 1; fileId <= 4; fileId++)
	{
		fileIdString[0] = '0' + fileId;
		fileIdString[1] = NULL;
		strcpy(fileName, "DxL_");
		strcat(fileName, fileIdString);
		strcat(fileName, ".csv");
		fid = fopen(fileName, "r");
		if (fid == NULL) continue;

		printf("Processing %s:\n", fileName);


		// number of data
		for (i = 0; i < 13; i++) fgets(strBuffer, MAXSTR, fid);
		numData = getInt();
		printf("    Number of data: %d\n", numData);


		// LAT
		printf("    LAT\n");
		strcpy(foutName, "LAT");
		strcat(foutName, fileIdString);
		strcat(foutName, ".csv");
		fout = fopen(foutName, "w");

		for (i = 0; i < 31; i++) fgets(strBuffer, MAXSTR, fid);
		fputs(strBuffer + 9, fout);
		fclose(fout);


		// peaktopeak
		printf("    peaktopeak\n");
		strcpy(foutName, "peaktopeak");
		strcat(foutName, fileIdString);
		strcat(foutName, ".csv");
		fout = fopen(foutName, "w");

		for (i = 0; i < 1; i++) fgets(strBuffer, MAXSTR, fid);
		fputs(strBuffer + 11, fout);
		fclose(fout);


		// CFAE
		printf("    CFAE\n");
		strcpy(foutName, "CFAE");
		strcat(foutName, fileIdString);
		strcat(foutName, ".csv");
		fout = fopen(foutName, "w");

		for (i = 0; i < 2; i++) fgets(strBuffer, MAXSTR, fid);
		fputs(strBuffer + 10, fout);
		fclose(fout);


		// xyz
		printf("    xyz\n");
		strcpy(foutName, "xyz");
		strcat(foutName, fileIdString);
		strcat(foutName, ".csv");
		fout = fopen(foutName, "w");

		for (i = 0; i < 5; i++) fgets(strBuffer, MAXSTR, fid);
		fputs(strBuffer + 10, fout);
		fgets(strBuffer, MAXSTR, fid);
		fputs(strBuffer + 10, fout);
		fgets(strBuffer, MAXSTR, fid);
		fputs(strBuffer + 10, fout);
		fclose(fout);


		// utilized
		printf("    utilized\n");
		strcpy(foutName, "utilized");
		strcat(foutName, fileIdString);
		strcat(foutName, ".csv");
		fout = fopen(foutName, "w");

		for (i = 0; i < 3; i++) fgets(strBuffer, MAXSTR, fid);
		fputs(strBuffer + 10, fout);
		fclose(fout);


		// egm
		printf("    egm\n");
		strcpy(foutName, "egm");
		strcat(foutName, fileIdString);
		strcat(foutName, ".csv");
		fout = fopen(foutName, "w");

		for (i = 0; i < 13; i++) fgets(strBuffer, MAXSTR, fid);
		for (i = 0; i < 10173; i++)
		{
			fgets(strBuffer, MAXSTR, fid);
			fputs(strBuffer, fout);
		}
		fclose(fout);


		fclose(fid);
	}
}