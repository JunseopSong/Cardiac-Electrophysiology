// ------------------------------ PROGRAM SETTING START ------------------------------
/*
couplingType
0: Dpar1 (SN)
1: Dpar2 (RA)
2: Exit Channel
-1: Fibrosis
Dpar_ATRIUM_SN: Atrial cell - SN cell Junction

cellModelType
-10: Atrial cell
0: SN cell
1: Perinodal cell
*/

#define RADIUS_CENTER 0.0
#define RADIUS_HETERO 1.0 // cm

#define RADIUS_EXIT_CHANNEL 3 // number of cell
//#define EXIT_CHANNEL_START_POINT 20

#define Dpar00 0.0002
#define Dpar11 0.002

/*
#define Dpar1 0.02 // mm^2/ms
#define Dpar2 0.20 // mm^2/ms
double Dpar_INTRA_CHANNEL = 0.25; // mm^2/ms
double Dpar_ATRIUM_SN = 0.25; // mm^2/ms
*/

const int probe[] = {10910,130994,37303};
const int PROBE_NUMBER = sizeof(probe) / sizeof(probe[0]);
//#define probe_x -5.5903372
//#define probe_y 90.1180948
//#define probe_z -7.2320996


//#define FIBROSIS_RATIO 0.3
//#define FIBROSIS
//#define CONDUCTION_PATHWAY

// ------------------------------ PROGRAM SETTING END --------------------------------


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <Windows.h>
#include <process.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


#define dt 0.02
#define itmax 150100 // time limit
#define tol 0.00001

#define MAX_DEG 15

////////////////////// CRN Constants Start //////////////////////
#define NUM_CRN_VAR 25

#define g_Na  7.8   // nanoS_per_picoF (in fast_sodium_current)
#define g_K1  0.09   // nanoS_per_picoF (in time_independent_potassium_current)
#define g_to  0.1652   // nanoS_per_picoF (in transient_outward_K_current)
#define g_Kr  0.029411765   // nanoS_per_picoF (in rapid_delayed_rectifier_K_current)
#define g_Ca_L  0.12375   // nanoS_per_picoF (in L_type_Ca_channel)
#define GKur 1.0 // dimensionless

#define G_f 1.25 // 3.76
#define E_f -22.0
#define G_CaT 11.0
#define E_CaT 45.0

#define CMDN_max  0.05   // millimolar (in Ca_buffers)
#define CSQN_max  10.0   // millimolar (in Ca_buffers)
#define Km_CMDN  0.00238   // millimolar (in Ca_buffers)
#define Km_CSQN  0.8   // millimolar (in Ca_buffers)
#define Km_TRPN  0.0005   // millimolar (in Ca_buffers)
#define TRPN_max  0.07   // millimolar (in Ca_buffers)
#define Ca_up_max  15.0   // millimolar (in Ca_leak_current_by_the_NSR)
#define K_rel  30.0   // per_millisecond (in Ca_release_current_from_JSR)
#define I_up_max  0.005   // millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)
#define K_up  0.00092   // millimolar (in Ca_uptake_current_by_the_NSR)
#define I_NaCa_max  1600.0   // picoA_per_picoF (in Na_Ca_exchanger_current)
#define K_mCa  1.38   // millimolar (in Na_Ca_exchanger_current)
#define K_mNa  87.5   // millimolar (in Na_Ca_exchanger_current)
#define K_sat  0.1   // dimensionless (in Na_Ca_exchanger_current)
#define gamma  0.35   // dimensionless (in Na_Ca_exchanger_current)
#define g_B_Ca  0.001131   // nanoS_per_picoF (in background_currents)
#define g_B_K  0.0   // nanoS_per_picoF (in background_currents)
#define g_B_Na  0.0006744375   // nanoS_per_picoF (in background_currents)
#define V_cell  20100.0   // micrometre_3 (in intracellular_ion_concentrations)
#define Cm  100.0   // picoF (in membrane)
#define F  96.4867   // coulomb_per_millimole (in membrane)
#define R  8.3143   // joule_per_mole_kelvin (in membrane)
#define T  310.0   // kelvin (in membrane)
#define RT_F 26.712832 // R*T/F
#define stim_amplitude  -2800.0   // picoA (in membrane)
#define i_CaP_max  0.275   // picoA_per_picoF (in sarcolemmal_calcium_pump_current)
#define g_Ks  0.12941176   // nanoS_per_picoF (in slow_delayed_rectifier_K_current)
#define Km_K_o  1.5   // millimolar (in sodium_potassium_pump)
#define Km_Na_i  10.0   // millimolar (in sodium_potassium_pump)
#define i_NaK_max  0.59933874   // picoA_per_picoF (in sodium_potassium_pump)
#define Ca_o  1.8   // millimolar (in standard_ionic_concentrations)
#define K_o  5.4   // millimolar (in standard_ionic_concentrations)
#define Na_o  140.0   // millimolar (in standard_ionic_concentrations)
#define tau_tr  180.0   // millisecond (in transfer_current_from_NSR_to_JSR)
#define K_Q10  3.0   // dimensionless (in transient_outward_K_current)

// V_rel = 0.0048*V_cell;
// V_i = V_cell*0.68;
// V_up = 0.0552*V_cell;
// sigma = 1.0/7.0*(exp(Na_o/67.3)-1.0);

#define tau_u  8.0     // millisecond (in Ca_release_current_from_JSR_u_gate)
#define tau_f_Ca  2.0   // millisecond (in L_type_Ca_channel_f_Ca_gate)
#define V_i  13668.0   // micrometre_3 (in intracellular_ion_concentrations)
#define V_rel  96.48   // micrometre_3 (in intracellular_ion_concentrations)
#define V_up  1109.52   // micrometre_3 (in intracellular_ion_concentrations)
#define sigma  1.00091   // dimensionless (in sodium_potassium_pump)
////////////////////// CRN Constants END //////////////////////


// ---------------- CRN Lookup Table ----------------
#define CRNTableVolt 32768
#define vlo (-90.0)
#define dvt 0.005

// If you change dt, you have to change following two parameters:
#define expdt_tau_u 0.997503122 // expdt_tau_u = exp(-dt/tau_u)
#define expdt_tau_f_Ca 0.990049834 // expdt_tau_f_Ca = exp(-dt/tau_f_Ca)
// --------------------------------------------------

char *couplingType;
double *yglob[NUM_CRN_VAR], yinitial[NUM_CRN_VAR], *volt_print, ttt, ttt_print, *cellModelType;

FILE *printVoltageFile;

int nde=1000000, nel=1000000, protocolType, nstim, *init_stim, *ie[3];
bool *Istim, *ifscar, *ifablation, icheck1 = false, icheck2 = false, idid = false, stimEnd = true;

double *volt, *xy[3], *volt_tecplot;
double ttt0, ttt_printVoltage, ttt_printTecplot;

int *degree, *linkNode;
double *weightEdge;


void importData()
{
	int i, j;
	double temp;
	FILE *datain;

	datain = fopen("cuda_aa.dat", "rb");
	fscanf(datain, "%d", &nde);
	fscanf(datain, "%d", &nde);
	printf("   Number of nodes : %d\n", nde);
	fclose(datain);

	for(i=0; i<NUM_CRN_VAR; i++) yglob[i] = (double*)malloc(sizeof(double)*nde);
	volt = (double*)malloc(sizeof(double)*nde);
	volt_tecplot = (double*)malloc(sizeof(double)*nde);
	Istim = (bool*)malloc(sizeof(bool)*nde);
	for(i=0; i<3; i++) xy[i] = (double*)malloc(sizeof(double)*nde);
	ifablation = (bool*)malloc(sizeof(bool)*nde);
	for(i=0; i<nde; i++) ifablation[i] = false;


	datain = fopen("cuda_ifscar.dat", "rb");
	printf("   Importing scar zone info...\n");
	ifscar = (bool*)malloc(sizeof(bool)*nde);
	for(i=0; i<nde; i++)
	{
		fscanf(datain, "%lf", &temp);
		ifscar[i] = (bool)((int)temp);
	}
	fclose(datain);

	/*if(ABLATION)
	{
		int tempin, abnodes;
		datain = fopen("AblationNode.dat", "rb");
		fscanf(datain, "%d", &abnodes);
		printf("   Importing ablation node info...\n");
		for(i=0; i<abnodes; i++)
		{
			fscanf(datain, "%d", &tempin);
			tempin--; // Fix Fortran, C index mismatch

			if(tempin < nde)
			{
				ifscar[tempin] = true;
				ifablation[tempin] = true;
			}
		
		fclose(datain);
	}*/

	datain = fopen("cuda_initstim.dat", "rb");
	fscanf(datain, "%lf", &temp); nstim = (int)temp;
	printf("   Importing stimulation site info...\n");
	init_stim = (int*)malloc(sizeof(int)*nstim);
	for(i=0; i<nstim; i++)
	{
		fscanf(datain, "%lf", &temp);
		init_stim[i] = (int)temp;
	}
	for(i=0; i<nstim; i++) init_stim[i]--; // Fix Fortran, C index mismatch
	fclose(datain);

	datain = fopen("cuda_xy.dat", "rb");
	printf("   Importing coordinate info...\n");
	for(i=0; i<nde; i++)
	{
		for(j=0; j<3; j++)
		{
			fscanf(datain, "%lf", &temp);
			xy[j][i] = temp;
		}
	}
	fclose(datain);

	datain = fopen("cuda_ie.dat", "rb");
	fscanf(datain, "%lf", &temp); nel = (int)temp;
	printf("   Importing mesh info...\n");
	for(i=0; i<3; i++) ie[i] = (int*)malloc(sizeof(int)*nel);
	for(i=0; i<nel; i++)
	{
		for(j=0; j<3; j++)
		{
			fscanf(datain, "%lf", &temp);
			ie[j][i] = (int)temp;
		}
	}
	fclose(datain);
}



void setIFDM()
{
	int i, j, k, k1, k2, v1, v2;
	bool ifNewNode;
	double temp, tmp1, tmp2, tmp3, *area, *area_triangle, *medianLine[3];

	int *deg; // degree of each node
	int *link[MAX_DEG]; // link[j][i] : j'th nodes linked to i
	double *weight[MAX_DEG]; // weight[j][i] : weight between i and link[j][i]
	double *lenTri[3]; // length of sides of each triangle(element)


	// Allocate
	deg = (int*)malloc(sizeof(int)*nde);
	for(i=0; i<MAX_DEG; i++) link[i] = (int*)malloc(sizeof(int)*nde);
	for(i=0; i<MAX_DEG; i++) weight[i] = (double*)malloc(sizeof(double)*nde);
	for(i=0; i<3; i++) lenTri[i] = (double*)malloc(sizeof(double)*nel);

	for(i=0; i<nde; i++) deg[i] = 0;
	area = (double*)malloc(sizeof(double)*nde); // area of Voronoi voxel
	area_triangle = (double*)malloc(sizeof(double)*nel); // area of triangle
	for(i=0; i<3; i++) medianLine[i] = (double*)malloc(sizeof(double)*nel); // length of median lines of each triangle


	double lengthMin;
	lengthMin = 1000000.0;

	// ------------------ Geometry ------------------
	for(i=0; i<nel; i++)
	{
		// Calculate length of sides of triangle(element)
		for(j=0; j<3; j++)
		{
			v1 = ie[(j+1)%3][i]-1;
			v2 = ie[(j+2)%3][i]-1;

			temp = 0.0;
			for(k=0; k<3; k++) temp += (xy[k][v1]-xy[k][v2])*(xy[k][v1]-xy[k][v2]);

			lenTri[j][i] = sqrt(temp);

			if(lenTri[j][i]<lengthMin) lengthMin = lenTri[j][i];
		}

		// Calculate area of triangle
		temp = (lenTri[0][i] + lenTri[1][i] + lenTri[2][i])/2.0;
		area_triangle[i] = sqrt(temp*(temp-lenTri[0][i])*(temp-lenTri[1][i])*(temp-lenTri[2][i]));

		// Calculate length of median lines
		for(j=0; j<3; j++)
		{
			tmp1 = lenTri[j][i];
			tmp2 = lenTri[(j+1)%3][i];
			tmp3 = lenTri[(j+2)%3][i];
			medianLine[j][i] = sqrt(tmp2*tmp2/2.0 + tmp3*tmp3/2.0 - tmp1*tmp1/4.0);
		}
	}

	//printf("%lf\n\n", lengthMin);
	//system("pause");


	// ------------------ Get graph structure ------------------
	for(i=0; i<nel; i++)
	{
		// Scan 1
		for(k1=0; k1<3; k1++)
		{
			k2 = (k1+1)%3;
			v1 = ie[k1][i]-1;
			v2 = ie[k2][i]-1;

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

		// Scan 2
		for(k1=0; k1<3; k1++)
		{
			k2 = (k1+2)%3;
			v1 = ie[k1][i]-1;
			v2 = ie[k2][i]-1;

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


	// ------------------ Calculate weight and area ------------------
	for(i=0; i<nde; i++)
	{
		area[i] = 0.0;
		for(j=0; j<deg[i]; j++) weight[j][i] = 0.0;
	}

	for(i=0; i<nel; i++)
	{
		for(j=0; j<3; j++) area[ie[j][i]-1] += area_triangle[i]/3.0;

		for(k1=0; k1<3; k1++)
		{
			k2 = (k1+1)%3;
			k = (k1+2)%3;
			v1 = ie[k1][i]-1;
			v2 = ie[k2][i]-1;

			// find index for v1
			for(j=0; j<deg[v1]; j++)
			{
				if(link[j][v1] == v2)
				{
					weight[j][v1] += (medianLine[k][i]/lenTri[k][i]/3.0);
					break;
				}
			}

			// find index for v2
			for(j=0; j<deg[v2]; j++)
			{
				if(link[j][v2] == v1)
				{
					weight[j][v2] += (medianLine[k][i]/lenTri[k][i]/3.0);
					break;
				}
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
		for(j=degree[i]; j<degree[i+1]; j++) linkNode[j] = link[j-degree[i]][i];
	}

	weightEdge = (double*)malloc(sizeof(double)*degree[nde]);
	for(i=0; i<nde; i++)
	{
		for(j=degree[i]; j<degree[i+1]; j++) weightEdge[j] = weight[j-degree[i]][i] / area[i];
	}


	// Deallocate
	free(deg);
	for(i=0; i<MAX_DEG; i++) free(link[i]);
	for(i=0; i<MAX_DEG; i++) free(weight[i]);
	for(i=0; i<3; i++) free(lenTri[i]);

	free(area);
	free(area_triangle);
	for(i=0; i<3; i++) free(medianLine[i]);
}


void preprocess()
{
	int i;
	double dist, probe_x, probe_y, probe_z;

	for(i=0; i<NUM_CRN_VAR; i++) yglob[i] = (double*)malloc(sizeof(double)*nde);
	volt = (double*)malloc(sizeof(double)*nde);
	volt_print = (double*)malloc(sizeof(double)*nde);
	couplingType = (char*)malloc(sizeof(char)*nde);
	cellModelType = (double*)malloc(sizeof(double)*nde);


	probe_x = xy[0][probe[0]];
	probe_y = xy[1][probe[0]];
	probe_z = xy[2][probe[0]];

	for(i=0; i<nde; i++) couplingType[i] = 0;


	// Gap junction
	for(i=0; i<nde; i++)
	{
		if((xy[0][i]-probe_x)*(xy[0][i]-probe_x) + (xy[1][i]-probe_y)*(xy[1][i]-probe_y) + (xy[2][i]-probe_z)*(xy[2][i]-probe_z) >= RADIUS_HETERO*RADIUS_HETERO)
		{
			couplingType[i] = 1;
		}
	}


	// Cell EP model hetero
	for(i=0; i<nde; i++)
	{
		dist = sqrt((double)((xy[0][i]-probe_x)*(xy[0][i]-probe_x) + (xy[1][i]-probe_y)*(xy[1][i]-probe_y) + (xy[2][i]-probe_z)*(xy[2][i]-probe_z)));

		if(dist < RADIUS_CENTER)
		{
			cellModelType[i] = 0; // SN cell
		}
		else if(dist < RADIUS_HETERO)
		{
			cellModelType[i] = (dist-RADIUS_CENTER) / (RADIUS_HETERO-RADIUS_CENTER); // D cell
			cellModelType[i] = 1.0309*cellModelType[i] / (1.0 + 0.7745*exp(-(3*cellModelType[i]-2.05)/0.295)); // F cell
		}
		else
		{
			cellModelType[i] = -10; // atrial cell
		}
	}


/*#ifdef FIBROSIS
	sampleNum = populationNum * FIBROSIS_RATIO;

	srand((unsigned)time(NULL));

	for(sampleCnt=0; sampleCnt<sampleNum; sampleCnt++)
	{
		i = rand()%M;
		j = rand()%M;

		distance2 = (xy[0][i]-probe_x)*(xy[0][i]-probe_x) + (xy[1][i]-probe_y)*(xy[1][i]-probe_y) + (xy[2][i]-probe_z)*(xy[2][i]-probe_z);
		if(distance2 <= RADIUS_HETERO*RADIUS_HETERO) couplingType[index(nde)] = -1;
		else sampleCnt--;
	}
#endif

#ifdef CONDUCTION_PATHWAY
	for(i=M/2+EXIT_CHANNEL_START_POINT; i<=M/2+RADIUS_HETERO; i++)
	{
		for(j=M/2-RADIUS_EXIT_CHANNEL; j<=M/2+RADIUS_EXIT_CHANNEL; j++)
		{
			couplingType[index(i, j)] = 2;
			cellModelType[index(i, j)] = -1; // atrial cell
		}
	}

	// Draw hemi-circle
	for(i=M/2+EXIT_CHANNEL_START_POINT-RADIUS_EXIT_CHANNEL; i<=M/2+EXIT_CHANNEL_START_POINT; i++)
	{
		for(j=M/2-RADIUS_EXIT_CHANNEL; j<=M/2+RADIUS_EXIT_CHANNEL; j++)
		{
			if((i-(M/2+EXIT_CHANNEL_START_POINT))*(i-(M/2+EXIT_CHANNEL_START_POINT)) + (j-M/2)*(j-M/2) <= RADIUS_EXIT_CHANNEL*RADIUS_EXIT_CHANNEL)
			{
				couplingType[index(i, j)] = 2;
				cellModelType[index(i, j)] = -1; // atrial cell
			}
		}
	}
#endif*/
}



void CRN_initialValue()
{
	int i, j;


	// SN cell
	yinitial[0] = 0.0;  // u (dimensionless) (in Ca_release_current_from_JSR_u_gate)
	yinitial[1] = 1.0;        // v (dimensionless) (in Ca_release_current_from_JSR_v_gate)
	yinitial[2] = 0.9926;     // w (dimensionless) (in Ca_release_current_from_JSR_w_gate)
	yinitial[3] = 0.0156;   // d (dimensionless) (in L_type_Ca_channel_d_gate)
	yinitial[4] = 0.3899;   // f_Ca (dimensionless) (in L_type_Ca_channel_f_Ca_gate)
	yinitial[5] = 0.4100;   // f (dimensionless) (in L_type_Ca_channel_f_gate)
	yinitial[6] = 0.0024;   // h (dimensionless) (in fast_sodium_current_h_gate)
	yinitial[7] = 0.0017;   // j (dimensionless) (in fast_sodium_current_j_gate)
	yinitial[8] = 0.4848;   // m (dimensionless) (in fast_sodium_current_m_gate)
	yinitial[9] = 0.0005;   // Ca_i (millimolar) (in intracellular_ion_concentrations)
	yinitial[10] = 4.9738;     // Ca_rel (millimolar) (in intracellular_ion_concentrations)
	yinitial[11] = 4.9837;     // Ca_up (millimolar) (in intracellular_ion_concentrations)
	yinitial[12] = 139.0;    // K_i (millimolar) (in intracellular_ion_concentrations)
	yinitial[13] = 9.4330;   // Na_i (millimolar) (in intracellular_ion_concentrations)
	yinitial[14] = -55.3145;    // V (millivolt) (in membrane)
	yinitial[15] = 0.4148;  // xr (dimensionless) (in rapid_delayed_rectifier_K_current_xr_gate)
	yinitial[16] = 0.1769;  // xs (dimensionless) (in slow_delayed_rectifier_K_current_xs_gate)
	yinitial[17] = 0.2238;  // oa (dimensionless) (in transient_outward_K_current_oa_gate)
	yinitial[18] = 0.2092;  // oi (dimensionless) (in transient_outward_K_current_oi_gate)
	yinitial[19] = 0.2233;  // ua (dimensionless) (in ultrarapid_delayed_rectifier_K_current_ua_gate)
	yinitial[20] = 0.9814;  // ui (dimensionless) (in ultrarapid_delayed_rectifier_K_current_ui_gate)
	yinitial[21] = 0.0038;  // fun (in funny current)
	yinitial[22] = 0.0561;  // dt (in T type Ca channel)
	yinitial[23] = 0.0335;  // ft (in T type Ca chennel)
	yinitial[24] = 0.1619;  // dnew (in L type Ca chennel)

	for(i=0; i<nde; i++)
	{
		for(j=0; j<NUM_CRN_VAR; j++) yglob[j][i] = yinitial[j];
		volt[i] = yinitial[14];
	}


	// Atrial cell
	yinitial[0] = 2.35e-112;  // u (dimensionless) (in Ca_release_current_from_JSR_u_gate)
	yinitial[1] = 1.0;        // v (dimensionless) (in Ca_release_current_from_JSR_v_gate)
	yinitial[2] = 0.9992;     // w (dimensionless) (in Ca_release_current_from_JSR_w_gate)
	yinitial[3] = 1.367e-4;   // d (dimensionless) (in L_type_Ca_channel_d_gate)
	yinitial[4] = 7.755e-1;   // f_Ca (dimensionless) (in L_type_Ca_channel_f_Ca_gate)
	yinitial[5] = 9.996e-1;   // f (dimensionless) (in L_type_Ca_channel_f_gate)
	yinitial[6] = 9.649e-1;   // h (dimensionless) (in fast_sodium_current_h_gate)
	yinitial[7] = 9.775e-1;   // j (dimensionless) (in fast_sodium_current_j_gate)
	yinitial[8] = 2.908e-3;   // m (dimensionless) (in fast_sodium_current_m_gate)
	yinitial[9] = 1.013e-4;   // Ca_i (millimolar) (in intracellular_ion_concentrations)
	yinitial[10] = 1.488;     // Ca_rel (millimolar) (in intracellular_ion_concentrations)
	yinitial[11] = 1.488;     // Ca_up (millimolar) (in intracellular_ion_concentrations)
	yinitial[12] = 1.39e2;    // K_i (millimolar) (in intracellular_ion_concentrations)
	yinitial[13] = 1.117e1;   // Na_i (millimolar) (in intracellular_ion_concentrations)
	yinitial[14] = -81.18;    // V (millivolt) (in membrane)
	yinitial[15] = 3.296e-5;  // xr (dimensionless) (in rapid_delayed_rectifier_K_current_xr_gate)
	yinitial[16] = 1.869e-2;  // xs (dimensionless) (in slow_delayed_rectifier_K_current_xs_gate)
	yinitial[17] = 3.043e-2;  // oa (dimensionless) (in transient_outward_K_current_oa_gate)
	yinitial[18] = 9.992e-1;  // oi (dimensionless) (in transient_outward_K_current_oi_gate)
	yinitial[19] = 4.966e-3;  // ua (dimensionless) (in ultrarapid_delayed_rectifier_K_current_ua_gate)
	yinitial[20] = 9.986e-1;  // ui (dimensionless) (in ultrarapid_delayed_rectifier_K_current_ui_gate)

	for(i=0; i<nde; i++)
	{
		if(abs(cellModelType[i] + 10) < tol)
		{
			for(j=0; j<NUM_CRN_VAR; j++) yglob[j][i] = yinitial[j];
			volt[i] = yinitial[14];
		}
	}


	// Perinodal cell
	yinitial[0] = 0.000000;
	yinitial[1] = 1.000000;
	yinitial[2] = 0.999100;
	yinitial[3] = 0.000200;
	yinitial[4] = 0.758600;
	yinitial[5] = 0.999300;
	yinitial[6] = 0.938500;
	yinitial[7] = 0.956900;
	yinitial[8] = 0.004500;
	yinitial[9] = 0.000100;
	yinitial[10] = 1.639300;
	yinitial[11] = 1.637300;
	yinitial[12] = 139.560000;
	yinitial[13] = 10.191100;
	yinitial[14] = -78.571400;
	yinitial[15] = 0.000000;
	yinitial[16] = 0.020700;
	yinitial[17] = 0.035100;
	yinitial[18] = 0.998800;
	yinitial[19] = 0.006500;
	yinitial[20] = 0.998500;
	yinitial[21] = 0.192000;
	yinitial[22] = 0.000200;
	yinitial[23] = 0.698700;
	yinitial[24] = 0.002300;

	for(i=0; i<nde; i++)
	{
		if(cellModelType[i] == 1)
		{
			for(j=0; j<NUM_CRN_VAR; j++) yglob[j][i] = yinitial[j];
			volt[i] = yinitial[14];
		}
	}
}



__global__ void setCRNTableVolt(double *expdt_tau_m, double *m_inf, double *expdt_tau_h, double *h_inf, double *expdt_tau_j, double *j_inf, double *expdt_tau_xr, double *xr_infinity,
	double *expdt_tau_w, double *w_infinity, double *expdt_tau_xs, double *xs_infinity, double *expdt_tau_d, double *d_infinity, double *expdt_tau_oa, double *oa_infinity,
	double *expdt_tau_oi, double *oi_infinity, double *expdt_tau_ua, double *ua_infinity, double *expdt_tau_ui, double *ui_infinity, double *expdt_tau_f, double *f_infinity, double *j_NaK,
	double *j_NaCa1, double *j_NaCa2, double *g_Kur, double *j_K1, double *j_Kr,
	double *expdt_tau_fun, double *fun_infinity, double *expdt_tau_dt, double *dt_infinity, double *expdt_tau_ft, double *ft_infinity, double *d_infinity_new)
{
	double vv;

	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int bdim = blockDim.x;
	const unsigned int gdim = gridDim.x;
	int step = bdim*gdim;

	double alpha, beta;

	for(int i = bid*bdim + tid; i<CRNTableVolt; i+=step)
	{
		vv = vlo + (double)i * dvt;

		// m
		if (fabs(vv+47.13) < tol) alpha = 3.2;
		else alpha = 0.32*(vv+47.13)/(1.0-exp(-0.1*(vv+47.13)));
		beta = 0.08*exp(-vv/11.0);
		m_inf[i] = alpha/(alpha+beta);
		expdt_tau_m[i] = exp(-dt/(1.0/(alpha+beta)));

		// h
		if (vv < -40.0) alpha = 0.135*exp((vv+80.0)/-6.8);
		else alpha = 0.0;
		if (vv < -40.0) beta = 3.56*exp(0.079*vv)+3.1e5*exp(0.35*vv);
		else beta = 1.0/(0.13*(1.0+exp((vv+10.66)/-11.1)));
		h_inf[i] = alpha/(alpha+beta);
		expdt_tau_h[i] = exp(-dt/(1.0/(alpha+beta)));

		// j
		if (vv < -40.0) alpha = (-1.2714e5*exp(0.2444*vv)-3.474e-5*exp(-0.04391*vv))*(vv+37.78)/(1.0+exp(0.311*(vv+79.23)));
		else alpha = 0.0;
		if (vv < -40.0) beta = 0.1212*exp(-0.01052*vv)/(1.0+exp(-0.1378*(vv+40.14)));
		else beta = 0.3*exp(-2.535e-7*vv)/(1.0+exp(-0.1*(vv+32.0)));
		j_inf[i] = alpha/(alpha+beta);
		expdt_tau_j[i] = exp(-dt/(1.0/(alpha+beta)));

		// xr
		if (fabs(vv+14.1) < tol) alpha = 0.0015;
		else alpha = 0.0003*(vv+14.1)/(1.0-exp((vv+14.1)/-5.0));
		if (fabs(vv-3.3328) < tol) beta = 3.7836118e-4;
		else beta = 0.000073898*(vv-3.3328)/(exp((vv-3.3328)/5.1237)-1.0);
		expdt_tau_xr[i] = exp(-dt/(1/(alpha+beta)));
		xr_infinity[i] = 1/(1.0+exp((vv+14.1)/-6.5));

		// w
		if (fabs(vv-7.9) < tol) expdt_tau_w[i] = exp(-dt/(6.0*0.2/1.3));
		else expdt_tau_w[i] = exp(-dt/(6.0*(1.0-exp(-(vv-7.9)/5.0))/((1.0+0.3*exp(-(vv-7.9)/5.0))*1.0*(vv-7.9))));
		w_infinity[i] = 1.0-1.0/(1.0+exp(-(vv-40.0)/17.0));

		// xs
		if (fabs(vv-19.9) < tol) alpha = 0.00068;
		else alpha = 0.00004*(vv-19.9)/(1.0-exp((vv-19.9)/-17.0));
		if (fabs(vv-19.9) < tol) beta = 0.000315;
		else beta = 0.000035*(vv-19.9)/(exp((vv-19.9)/9.0)-1.0);
		expdt_tau_xs[i] = exp(-dt/(0.5*1/(alpha+beta)));
		xs_infinity[i] = pow(1.0+exp((vv-19.9)/-12.7), -0.5);

		// d
		d_infinity[i] = 1.0/(1.0+exp((vv+10.0)/-8.0));
		if (fabs(vv+10.0) < tol) expdt_tau_d[i] = exp(-dt/(4.579/(1.0+exp((vv+10.0)/-6.24))));
		else expdt_tau_d[i] = exp(-dt/((1.0-exp((vv+10.0)/-6.24))/(0.035*(vv+10.0)*(1.0+exp((vv+10.0)/-6.24)))));

		// oa
		alpha = 0.65*1.0/(exp((vv-(-10.0))/-8.5)+exp((vv-(-10.0)-40.0)/-59.0));
		beta = 0.65*1.0/(2.5+exp((vv-(-10.0)+72.0)/17.0));
		expdt_tau_oa[i] = exp(-dt/(1/(alpha+beta)/K_Q10));
		oa_infinity[i] = 1.0/(1.0+exp((vv-(-10.0)+10.47)/-17.54));

		// oi
		alpha = 1.0/(18.53+1.0*exp((vv-(-10.0)+103.7)/10.95));
		beta = 1.0/(35.56+1.0*exp((vv-(-10.0)-8.74)/-7.44));
		expdt_tau_oi[i] = exp(-dt/(1/(alpha+beta)/K_Q10));
		oi_infinity[i] = 1.0/(1.0+exp((vv-(-10.0)+33.1)/5.3));

		// ua
		alpha = 0.65*1.0/(exp((vv-(-10.0))/-8.5)+exp((vv-(-10.0)-40.0)/-59.0));
		beta = 0.65*1.0/(2.5+exp((vv-(-10.0)+72.0)/17.0));
		expdt_tau_ua[i] = exp(-dt/(1/(alpha+beta)/K_Q10));
		ua_infinity[i] = 1.0/(1.0+exp((vv-(-10.0)+20.3)/-9.6));

		// ui
		alpha = 1.0/(21.0+1.0*exp((vv-(-10.0)-195.0)/-28.0));
		beta = 1.0/exp((vv-(-10.0)-168.0)/-16.0);
		expdt_tau_ui[i] = exp(-dt/(1.0/(alpha+beta)/K_Q10));
		ui_infinity[i] = 1.0/(1.0+exp((vv-(-10.0)-109.45)/27.48));

		// f
		f_infinity[i] = exp(-(vv+28.0)/6.9)/(1.0+exp(-(vv+28.0)/6.9));
		expdt_tau_f[i] = exp(-dt/(9.0*1.0/(0.0197*exp(-0.0337*0.0337*(vv+10.0)*(vv+10.0))+0.02)));

		// j_NaK
		j_NaK[i] = Cm*i_NaK_max/(1.0+0.1245*exp(-0.1*F*vv/(R*T))+0.0365*sigma*exp(-F*vv/(R*T)))*K_o/(K_o+Km_K_o);

		// j_NaCa
		j_NaCa1[i] = Cm*I_NaCa_max*exp(gamma*F*vv/(R*T))*Ca_o/((K_mNa*K_mNa*K_mNa + Na_o*Na_o*Na_o)*(K_mCa+Ca_o)*(1.0+K_sat*exp((gamma-1.0)*vv*F/(R*T))));
		j_NaCa2[i] = Cm*I_NaCa_max*exp((gamma-1.0)*F*vv/(R*T))*Na_o*Na_o*Na_o/((K_mNa*K_mNa*K_mNa + Na_o*Na_o*Na_o)*(K_mCa+Ca_o)*(1.0+K_sat*exp((gamma-1.0)*vv*F/(R*T))));

		// g_Kur
		g_Kur[i] = Cm*GKur*(0.005+0.05/(1.0+exp((vv-15.0)/-13.0)));

		// j_K1
		j_K1[i] = Cm*g_K1/(1.0+exp(0.07*(vv+80.0)));

		// j_Kr
		j_Kr[i] = Cm*g_Kr/(1.0+exp((vv+15.0)/22.4));

		// I_f
		alpha = 0.36*(vv+148.8)/(exp(0.066*(vv+148.8))-1.0);
		beta = 0.1*(vv+87.3)/(1.0-exp(-0.21*(vv+87.3)));
		expdt_tau_fun[i] = exp(-dt/(1000.0/(alpha + beta) - 54.0));
		fun_infinity[i] = 1.0/(1.0+exp((vv+95.95)/12.1));

		// I_CaT
		expdt_tau_dt[i] = exp(-dt/(1000.0/(1068.0*(exp((vv+26.3)/30.0) + exp(-(vv+26.3)/30.0)))));
		dt_infinity[i] = 1.0/(1.0 + exp(-(vv+26.3)/6.0));
		expdt_tau_ft[i] = exp(-dt/(1000.0/(15.3*exp(-(vv+61.7)/83.3) + 15.0*exp((vv+61.7)/15.38))));
		ft_infinity[i] = 1.0/(1.0 + exp((vv+71.0)/9.0));

		// d_new
		d_infinity_new[i] = 1.0/(1.0+exp((vv+20.0+10.0)/-8.0));
	}
}



__global__ void CRNKernel(double* u, double* vc, double* w, double* d, double*  f_Ca, double* f, double* h, double* jj, double* m, double* Ca_i, double* Ca_rel, double* Ca_up, double* K_i, double* Na_i, double* v, double* xr, double* xs, double* oa, double* oi, double* ua, double* ui, double *fun, double *dt_, double *ft_, double *dnew,
	double* volt_, const int num, const double *cellModelType_d,
	const double *expdt_tau_m, const double *m_inf, const double *expdt_tau_h, const double *h_inf, const double *expdt_tau_j, const double *j_inf, const double *expdt_tau_xr, const double *xr_infinity,
	const double *expdt_tau_w, const double *w_infinity, const double *expdt_tau_xs, const double *xs_infinity, const double *expdt_tau_d, const double *d_infinity, const double *expdt_tau_oa, const double *oa_infinity,
	const double *expdt_tau_oi, const double *oi_infinity, const double *expdt_tau_ua, const double *ua_infinity, const double *expdt_tau_ui, const double *ui_infinity, const double *expdt_tau_f, const double *f_infinity, const double *j_NaK,
	const double *j_NaCa1, const double *j_NaCa2, const double *g_Kur, const double *j_K1, const double *j_Kr, const double *expdt_tau_fun, const double *fun_infinity, const double *expdt_tau_dt, const double *dt_infinity, const double *expdt_tau_ft, const double *ft_infinity, const double *d_infinity_new)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int bdim = blockDim.x;
	const unsigned int gdim = gridDim.x;

	int step = bdim*gdim;

	double qv1, qv2, temp1, temp2;
	int vindex1, vindex2;

	double Vm, Nai, Cai, Ki, tt;
	double i_up_leak, i_rel, i_Ca_L, i_NaCa, Fn, i_up, E_Ca, E_Na, i_B_Na, i_B_Ca, E_K, i_B_K, i_Na;
	double i_NaK, i_K1, i_to, i_Kur, i_Kr, i_Ks, i_CaP, i_tr;
	double i_f, i_CaT;

	for (int id = bid*bdim + tid; id<num; id+=step)
	{
		Vm = volt_[id];
		Nai = Na_i[id];
		Cai = Ca_i[id];
		Ki = K_i[id];

		v[id] = Vm;

		// Lookup table index
		temp1 = (Vm-vlo)/dvt;
		vindex1 = (int)temp1;
		vindex2 = vindex1 + 1;
		qv1 = temp1 - (double)vindex1;
		qv2 = (double)vindex2 - temp1;

		tt = cellModelType_d[id];
		
		if(abs(tt + 10) < tol) // atrial cell
		{
			i_up_leak = I_up_max*Ca_up[id]/Ca_up_max;

			i_rel = K_rel*u[id]*u[id]*vc[id]*w[id]*(Ca_rel[id]-Cai);

			i_Ca_L = Cm*g_Ca_L*d[id]*f[id]*f_Ca[id]*(Vm-65.0);

			temp1 = j_NaCa1[vindex1]*qv2 + j_NaCa1[vindex2]*qv1;
			temp2 = j_NaCa2[vindex1]*qv2 + j_NaCa2[vindex2]*qv1;
			i_NaCa = temp1*Nai*Nai*Nai - temp2*Cai;

			Fn = 1.0e-12*V_rel*i_rel-5.0e-13/F*(0.5*i_Ca_L-0.2*i_NaCa);

			temp1 = 1.0/(1.0+exp(-(Fn-3.4175e-13)/13.67e-16));
			u[id] = temp1-(temp1-u[id])*expdt_tau_u;

			temp1 = 1.0-1.0/(1.0+exp(-(Fn-6.835e-14)/13.67e-16));
			temp2 = 1.91+2.09*1.0/(1.0+exp(-(Fn-3.4175e-13)/13.67e-16));
			vc[id] = temp1-(temp1-vc[id])*exp(-dt/temp2);

			temp1 = w_infinity[vindex1]*qv2 + w_infinity[vindex2]*qv1;
			temp2 = expdt_tau_w[vindex1]*qv2 + expdt_tau_w[vindex2]*qv1;
			w[id] = temp1-(temp1-w[id])*temp2;

			i_up = I_up_max/(1.0+K_up/Cai);


			temp1 = d_infinity[vindex1]*qv2 + d_infinity[vindex2]*qv1;
			temp2 = expdt_tau_d[vindex1]*qv2 + expdt_tau_d[vindex2]*qv1;
			d[id] = temp1-(temp1-d[id])*temp2;

			temp1 = 1.0/(1.0+Cai/0.00035);
			f_Ca[id] = temp1-(temp1-f_Ca[id])*expdt_tau_f_Ca;

			temp1 = f_infinity[vindex1]*qv2 + f_infinity[vindex2]*qv1;
			temp2 = expdt_tau_f[vindex1]*qv2 + expdt_tau_f[vindex2]*qv1;
			f[id] = temp1-(temp1-f[id])*temp2;

			E_Ca = RT_F/2.0*log(Ca_o/Cai);
			E_Na = RT_F*log(Na_o/Nai);
			i_B_Na = Cm*g_B_Na*(Vm-E_Na);
			i_B_Ca = Cm*g_B_Ca*(Vm-E_Ca);
			E_K = RT_F*log(K_o/Ki);
			i_B_K = Cm*g_B_K*(Vm-E_K);

			i_Na = Cm*g_Na*m[id]*m[id]*m[id]*h[id]*jj[id]*(Vm-E_Na);


			temp1 = h_inf[vindex1]*qv2 + h_inf[vindex2]*qv1;
			temp2 = expdt_tau_h[vindex1]*qv2 + expdt_tau_h[vindex2]*qv1;
			h[id] = temp1-(temp1-h[id])*temp2;

			temp1 = j_inf[vindex1]*qv2 + j_inf[vindex2]*qv1;
			temp2 = expdt_tau_j[vindex1]*qv2 + expdt_tau_j[vindex2]*qv1;
			jj[id] = temp1-(temp1-jj[id])*temp2;

			temp1 = m_inf[vindex1]*qv2 + m_inf[vindex2]*qv1;
			temp2 = expdt_tau_m[vindex1]*qv2 + expdt_tau_m[vindex2]*qv1;
			m[id] = temp1-(temp1-m[id])*temp2;

			temp1 = j_NaK[vindex1]*qv2 + j_NaK[vindex2]*qv1;
			temp2 = Km_Na_i/Nai;
			i_NaK = temp1 / (1.0+temp2*sqrt(temp2));

			temp1 = j_K1[vindex1]*qv2 + j_K1[vindex2]*qv1;
			i_K1 = temp1*(Vm-E_K);

			i_to = Cm*g_to*oa[id]*oa[id]*oa[id]*oi[id]*(Vm-E_K);

			temp1 = g_Kur[vindex1]*qv2 + g_Kur[vindex2]*qv1;
			i_Kur = temp1*ua[id]*ua[id]*ua[id]*ui[id]*(Vm-E_K);

			temp1 = j_Kr[vindex1]*qv2 + j_Kr[vindex2]*qv1;
			i_Kr = temp1*xr[id]*(Vm-E_K);

			i_Ks = Cm*g_Ks*xs[id]*xs[id]*(Vm-E_K);
			i_CaP = Cm*i_CaP_max*Cai/(0.0005+Cai);
			i_tr = (Ca_up[id]-Ca_rel[id])/tau_tr;


			temp1 = xr_infinity[vindex1]*qv2 + xr_infinity[vindex2]*qv1;
			temp2 = expdt_tau_xr[vindex1]*qv2 + expdt_tau_xr[vindex2]*qv1;
			xr[id] = temp1-(temp1-xr[id])*temp2;

			temp1 = xs_infinity[vindex1]*qv2 + xs_infinity[vindex2]*qv1;
			temp2 = expdt_tau_xs[vindex1]*qv2 + expdt_tau_xs[vindex2]*qv1;
			xs[id] = temp1-(temp1-xs[id])*temp2;

			temp1 = oa_infinity[vindex1]*qv2 + oa_infinity[vindex2]*qv1;
			temp2 = expdt_tau_oa[vindex1]*qv2 + expdt_tau_oa[vindex2]*qv1;
			oa[id] = temp1-(temp1-oa[id])*temp2;

			temp1 = oi_infinity[vindex1]*qv2 + oi_infinity[vindex2]*qv1;
			temp2 = expdt_tau_oi[vindex1]*qv2 + expdt_tau_oi[vindex2]*qv1;
			oi[id] = temp1-(temp1-oi[id])*temp2;

			temp1 = ua_infinity[vindex1]*qv2 + ua_infinity[vindex2]*qv1;
			temp2 = expdt_tau_ua[vindex1]*qv2 + expdt_tau_ua[vindex2]*qv1;
			ua[id] = temp1-(temp1-ua[id])*temp2;

			temp1 = ui_infinity[vindex1]*qv2 + ui_infinity[vindex2]*qv1;
			temp2 = expdt_tau_ui[vindex1]*qv2 + expdt_tau_ui[vindex2]*qv1;
			ui[id] = temp1-(temp1-ui[id])*temp2;

			i_f = 0.0;

			i_CaT = 0.0;

			temp1 = (2.0*i_NaCa-(i_CaP+i_Ca_L+i_CaT+i_B_Ca))/(2.0*V_i*F)+(V_up*(i_up_leak-i_up)+i_rel*V_rel)/V_i; // B1
			temp2 = 1.0+TRPN_max*Km_TRPN/((Cai+Km_TRPN)*(Cai+Km_TRPN))+CMDN_max*Km_CMDN/((Cai+Km_CMDN)*(Cai+Km_CMDN)); // B2


			Vm -= (i_Na+i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_Na+i_B_Ca+i_NaK+i_CaP+i_NaCa+i_Ca_L+i_f+i_CaT)/Cm*dt;

			Cai += (temp1/temp2)*dt;

			Ca_up[id] += (i_up-(i_up_leak+i_tr*V_rel/V_up))*dt;

			Ca_rel[id] += ((i_tr-i_rel)*1.0/(1.0+CSQN_max*Km_CSQN/((Ca_rel[id]+Km_CSQN)*(Ca_rel[id]+Km_CSQN))))*dt;

			Nai += (-3.0*i_NaK-(3.0*i_NaCa+i_B_Na+i_Na))/(V_i*F)*dt;

			Ki += (2.0*i_NaK-(i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_K))/(V_i*F)*dt;
		}
		else
		{
			i_up_leak = (0.30*(1-tt) + 0.96*tt)*I_up_max*Ca_up[id]/Ca_up_max;

			i_rel = K_rel*u[id]*u[id]*vc[id]*w[id]*(Ca_rel[id]-Cai);

			i_Ca_L = (0.68*0.5*(1-tt) + 2.584*0.75*tt)*Cm*g_Ca_L*d[id]*f[id]*f_Ca[id]*(Vm-65.0);

			temp1 = j_NaCa1[vindex1]*qv2 + j_NaCa1[vindex2]*qv1;
			temp2 = j_NaCa2[vindex1]*qv2 + j_NaCa2[vindex2]*qv1;
			i_NaCa = 0.5*(0.5*(1-tt) + 0.75*tt)*(temp1*Nai*Nai*Nai - temp2*Cai);

			Fn = 1.0e-12*(V_rel*0.5)*i_rel-5.0e-13/F*(0.5*i_Ca_L-0.2*i_NaCa);

			temp1 = 1.0/(1.0+exp(-(Fn-3.4175e-13)/13.67e-16));
			u[id] = temp1-(temp1-u[id])*expdt_tau_u;

			temp1 = 1.0-1.0/(1.0+exp(-(Fn-6.835e-14)/13.67e-16));
			temp2 = 1.91+2.09*1.0/(1.0+exp(-(Fn-3.4175e-13)/13.67e-16));
			vc[id] = temp1-(temp1-vc[id])*exp(-dt/temp2);

			temp1 = w_infinity[vindex1]*qv2 + w_infinity[vindex2]*qv1;
			temp2 = expdt_tau_w[vindex1]*qv2 + expdt_tau_w[vindex2]*qv1;
			w[id] = temp1-(temp1-w[id])*temp2;

			i_up = I_up_max/(1.0+K_up/Cai);


			temp1 = d_infinity[vindex1]*qv2 + d_infinity[vindex2]*qv1;
			temp2 = expdt_tau_d[vindex1]*qv2 + expdt_tau_d[vindex2]*qv1;
			d[id] = temp1-(temp1-d[id])*temp2;
			temp1 = d_infinity_new[vindex1]*qv2 + d_infinity_new[vindex2]*qv1;
			dnew[id] = temp1-(temp1-dnew[id])*temp2;

			temp1 = 1.0/(1.0+Cai/0.00035);
			f_Ca[id] = temp1-(temp1-f_Ca[id])*expdt_tau_f_Ca;

			temp1 = f_infinity[vindex1]*qv2 + f_infinity[vindex2]*qv1;
			temp2 = expdt_tau_f[vindex1]*qv2 + expdt_tau_f[vindex2]*qv1;
			f[id] = temp1-(temp1-f[id])*temp2;

			E_Ca = RT_F/2.0*log(Ca_o/Cai);
			E_Na = RT_F*log(Na_o/Nai);
			i_B_Na = (Cm*(0.5*(1-tt) + 0.75*tt))*g_B_Na*(Vm-E_Na);
			i_B_Ca = (Cm*(0.5*(1-tt) + 0.75*tt))*g_B_Ca*(Vm-E_Ca);
			E_K = RT_F*log(K_o/Ki);
			i_B_K = (Cm*(0.5*(1-tt) + 0.75*tt))*g_B_K*(Vm-E_K);

			i_Na = (0.06*0.5*(1-tt) + 0.30*0.75*tt)*Cm*g_Na*m[id]*m[id]*m[id]*h[id]*jj[id]*(Vm-E_Na);


			temp1 = h_inf[vindex1]*qv2 + h_inf[vindex2]*qv1;
			temp2 = expdt_tau_h[vindex1]*qv2 + expdt_tau_h[vindex2]*qv1;
			h[id] = temp1-(temp1-h[id])*temp2;

			temp1 = j_inf[vindex1]*qv2 + j_inf[vindex2]*qv1;
			temp2 = expdt_tau_j[vindex1]*qv2 + expdt_tau_j[vindex2]*qv1;
			jj[id] = temp1-(temp1-jj[id])*temp2;

			temp1 = m_inf[vindex1]*qv2 + m_inf[vindex2]*qv1;
			temp2 = expdt_tau_m[vindex1]*qv2 + expdt_tau_m[vindex2]*qv1;
			m[id] = temp1-(temp1-m[id])*temp2;

			temp1 = (0.5*(1-tt) + 0.75*tt)*(j_NaK[vindex1]*qv2 + j_NaK[vindex2]*qv1);
			temp2 = Km_Na_i/Nai;
			i_NaK = temp1 / (1.0+temp2*sqrt(temp2));

			temp1 = (0.5*(1-tt) + 0.75*tt)*(j_K1[vindex1]*qv2 + j_K1[vindex2]*qv1);
			i_K1 = 0.42*temp1*(Vm-E_K);

			i_to = (0.40*0.5*(1-tt) + 0.92*0.75*tt)*Cm*g_to*oa[id]*oa[id]*oa[id]*oi[id]*(Vm-E_K);

			temp1 = (0.5*(1-tt) + 0.75*tt)*(g_Kur[vindex1]*qv2 + g_Kur[vindex2]*qv1);
			i_Kur = 0.26*temp1*ua[id]*ua[id]*ua[id]*ui[id]*(Vm-E_K);

			temp1 = (0.5*(1-tt) + 0.75*tt)*(j_Kr[vindex1]*qv2 + j_Kr[vindex2]*qv1);
			i_Kr = (0.45*(1-tt) + 2.745*tt)*temp1*xr[id]*(Vm-E_K);

			i_Ks = (0.69*0.5*(1-tt) + 4.1883*0.75*tt)*Cm*g_Ks*xs[id]*xs[id]*(Vm-E_K);
			i_CaP = (Cm*(0.5*(1-tt) + 0.75*tt))*i_CaP_max*Cai/(0.0005+Cai);
			i_tr = (Ca_up[id]-Ca_rel[id])/tau_tr;


			temp1 = xr_infinity[vindex1]*qv2 + xr_infinity[vindex2]*qv1;
			temp2 = expdt_tau_xr[vindex1]*qv2 + expdt_tau_xr[vindex2]*qv1;
			xr[id] = temp1-(temp1-xr[id])*temp2;

			temp1 = xs_infinity[vindex1]*qv2 + xs_infinity[vindex2]*qv1;
			temp2 = expdt_tau_xs[vindex1]*qv2 + expdt_tau_xs[vindex2]*qv1;
			xs[id] = temp1-(temp1-xs[id])*temp2;

			temp1 = oa_infinity[vindex1]*qv2 + oa_infinity[vindex2]*qv1;
			temp2 = expdt_tau_oa[vindex1]*qv2 + expdt_tau_oa[vindex2]*qv1;
			oa[id] = temp1-(temp1-oa[id])*temp2;

			temp1 = oi_infinity[vindex1]*qv2 + oi_infinity[vindex2]*qv1;
			temp2 = expdt_tau_oi[vindex1]*qv2 + expdt_tau_oi[vindex2]*qv1;
			oi[id] = temp1-(temp1-oi[id])*temp2;

			temp1 = ua_infinity[vindex1]*qv2 + ua_infinity[vindex2]*qv1;
			temp2 = expdt_tau_ua[vindex1]*qv2 + expdt_tau_ua[vindex2]*qv1;
			ua[id] = temp1-(temp1-ua[id])*temp2;

			temp1 = ui_infinity[vindex1]*qv2 + ui_infinity[vindex2]*qv1;
			temp2 = expdt_tau_ui[vindex1]*qv2 + expdt_tau_ui[vindex2]*qv1;
			ui[id] = temp1-(temp1-ui[id])*temp2;


			temp1 = fun_infinity[vindex1]*qv2 + fun_infinity[vindex2]*qv1;
			temp2 = expdt_tau_fun[vindex1]*qv2 + expdt_tau_fun[vindex2]*qv1;
			fun[id] = temp1-(temp1-fun[id])*temp2;
			i_f = (7.5*(1-tt) + 29.475*tt)*G_f*fun[id]*(Vm-(E_f));


			temp1 = dt_infinity[vindex1-500]*qv2 + dt_infinity[vindex2-500]*qv1;
			temp2 = expdt_tau_dt[vindex1]*qv2 + expdt_tau_dt[vindex2]*qv1;
			dt_[id] = temp1-(temp1-dt_[id])*temp2;

			temp1 = ft_infinity[vindex1]*qv2 + ft_infinity[vindex2]*qv1;
			temp2 = expdt_tau_ft[vindex1]*qv2 + expdt_tau_ft[vindex2]*qv1;
			ft_[id] = temp1-(temp1-ft_[id])*temp2;

			i_CaT = 13.2*G_CaT*dt_[id]*ft_[id]*(Vm-E_CaT);


			temp1 = (2.0*i_NaCa-(i_CaP+i_Ca_L+i_CaT+i_B_Ca))/(2.0*(V_i*(0.5*(1-tt) + 0.75*tt))*F)+((V_up*(0.5*(1-tt) + 0.75*tt))*(i_up_leak-i_up)+i_rel*(V_rel*(0.5*(1-tt) + 0.75*tt)))/(V_i*(0.5*(1-tt) + 0.75*tt)); // B1
			temp2 = 1.0+TRPN_max*Km_TRPN/((Cai+Km_TRPN)*(Cai+Km_TRPN))+CMDN_max*Km_CMDN/((Cai+Km_CMDN)*(Cai+Km_CMDN)); // B2


			Vm -= (i_Na+i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_Na+i_B_Ca+i_NaK+i_CaP+i_NaCa+i_Ca_L+i_f+i_CaT)/(Cm*(0.5*(1-tt) + 0.75*tt))*dt;

			Cai += (temp1/temp2)*dt;

			Ca_up[id] += (i_up-(i_up_leak+i_tr*(V_rel*(0.5*(1-tt) + 0.75*tt))/(V_up*(0.5*(1-tt) + 0.75*tt))))*dt;

			Ca_rel[id] += ((i_tr-i_rel)*1.0/(1.0+CSQN_max*Km_CSQN/((Ca_rel[id]+Km_CSQN)*(Ca_rel[id]+Km_CSQN))))*dt;

			Nai += (-3.0*i_NaK-(3.0*i_NaCa+i_B_Na+i_Na))/((V_i*(0.5*(1-tt) + 0.75*tt))*F)*dt;

			Ki += (2.0*i_NaK-(i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_K))/((V_i*(0.5*(1-tt) + 0.75*tt))*F)*dt;
		}


		volt_[id] = Vm;
		Na_i[id] = Nai;
		Ca_i[id] = Cai;
		K_i[id] = Ki;
	}
}



__global__ void diffKernel(double *volt_after, double *volt_before, const bool *nonExcitable, const double *cellModelType_d, const int *linkNode_, const int *degree_, const double *weight_, const int num)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int bdim = blockDim.x;
	const unsigned int gdim = gridDim.x;
	int step=bdim*gdim, i;

	double volt, sum1, sum2, tt;

	for(int id=bid * bdim + tid; id<num; id+=step)
	{
		volt = volt_before[id];
		volt_after[id] = volt;

		if(nonExcitable[id]) continue;

		sum1 = 0.0; sum2 = 0.0;
		for(i=degree_[id]; i<degree_[id+1]; i++)
		{
			if(nonExcitable[linkNode_[i]]) sum1 += (volt_before[linkNode_[i]] - volt)*weight_[i];
			else 
			{
				tt = (cellModelType_d[id] + cellModelType_d[linkNode_[i]]) / 2.0;

				if(tt >= 0) sum2 += (Dpar00*(1-tt) + Dpar11*tt)*(volt_before[linkNode_[i]] - volt)*weight_[i];
				else sum2 += Dpar11*(volt_before[linkNode_[i]] - volt)*weight_[i];
			}
		}

		volt_after[id] += (sum2*dt/2.0 + sum1*0.00001*dt/2.0);
	}
}



void printVoltage()
{
	fprintf(printVoltageFile, "%.1lf ", ttt_printVoltage);
	for(int i=0; i<PROBE_NUMBER; i++) fprintf(printVoltageFile, "%.2lf ", volt[probe[i]]);
	fprintf(printVoltageFile, "\n");
}


unsigned int __stdcall printTecplot(void* arg)
{
	FILE *fout;
	char stringarr[30];
	int i;

	sprintf(stringarr, "vm%05d.plt", (int)(ttt_printTecplot+0.5));
	fout = fopen(stringarr, "wb");

	fprintf(fout, "VARIABLES = \"X\", \"Y\", \"Z\", \"Vm\"\n");
	fprintf(fout, "ZONE F=FEPOINT, ET=triangle, N=%d , E=%d\n", nde, nel);

	for(i=0; i<nde; i++)
	{
		if(ifablation[i]) volt_tecplot[i] = 100.0;
		fprintf(fout, "%lf %lf %lf %.2lf\n", xy[0][i], xy[1][i], xy[2][i], volt_tecplot[i]);
	}

	for(i=0; i<nel; i++) fprintf(fout, "%d %d %d\n", ie[0][i], ie[1][i], ie[2][i]);

	fclose(fout);

	printf("\n%s is created.\n\n", stringarr);

	return 0;
}


int main()
{
	int i, j, iter, ivmaxn;
	double vmax;

	printf("\n\nImporting data...\n\n");
	importData();
	printf("\n\nPreparing...");
	setIFDM();

	printVoltageFile = fopen("Voltage_LA.dat", "wb");


	///////////////// CUDA Threads, Blocks Setting Start ///////////////
	dim3 grid(128, 1, 1);
	dim3 threads(128, 1, 1);
	///////////////// CUDA Threads, Blocks Setting End /////////////////


	preprocess();

	// Variables for CPU thread
	bool PRINT_TECPLOT_FIRST = true;
	unsigned int printTecplotId;
	HANDLE printTecplotHandle;

	/*unsigned int printVoltageId;
	HANDLE printVoltageHandle;*/


	// ============================ Device Variables ============================
	char *couplingType_d;
	double *volt_d, *volt_d_, *cellModelType_d;
	bool *stimulation, *ifscar_d;
	double *yglob0, *yglob1, *yglob2, *yglob3, *yglob4, *yglob5, *yglob6, *yglob7, *yglob8, *yglob9, *yglob10, *yglob11, *yglob12, *yglob13, *yglob14, *yglob15, *yglob16, *yglob17, *yglob18, *yglob19, *yglob20, *yglob21, *yglob22, *yglob23, *yglob24;
	
	cudaMalloc((void**)&yglob0, sizeof(double)*nde);
	cudaMalloc((void**)&yglob1, sizeof(double)*nde);
	cudaMalloc((void**)&yglob2, sizeof(double)*nde);
	cudaMalloc((void**)&yglob3, sizeof(double)*nde);
	cudaMalloc((void**)&yglob4, sizeof(double)*nde);
	cudaMalloc((void**)&yglob5, sizeof(double)*nde);
	cudaMalloc((void**)&yglob6, sizeof(double)*nde);
	cudaMalloc((void**)&yglob7, sizeof(double)*nde);
	cudaMalloc((void**)&yglob8, sizeof(double)*nde);
	cudaMalloc((void**)&yglob9, sizeof(double)*nde);
	cudaMalloc((void**)&yglob10, sizeof(double)*nde);
	cudaMalloc((void**)&yglob11, sizeof(double)*nde);
	cudaMalloc((void**)&yglob12, sizeof(double)*nde);
	cudaMalloc((void**)&yglob13, sizeof(double)*nde);
	cudaMalloc((void**)&yglob14, sizeof(double)*nde);
	cudaMalloc((void**)&yglob15, sizeof(double)*nde);
	cudaMalloc((void**)&yglob16, sizeof(double)*nde);
	cudaMalloc((void**)&yglob17, sizeof(double)*nde);
	cudaMalloc((void**)&yglob18, sizeof(double)*nde);
	cudaMalloc((void**)&yglob19, sizeof(double)*nde);
	cudaMalloc((void**)&yglob20, sizeof(double)*nde);
	cudaMalloc((void**)&yglob21, sizeof(double)*nde);
	cudaMalloc((void**)&yglob22, sizeof(double)*nde);
	cudaMalloc((void**)&yglob23, sizeof(double)*nde);
	cudaMalloc((void**)&yglob24, sizeof(double)*nde);

	cudaMalloc((void**)&volt_d, sizeof(double)*nde);
	cudaMalloc((void**)&volt_d_, sizeof(double)*nde);
	cudaMalloc((void**)&couplingType_d, sizeof(char)*nde);
	cudaMalloc((void**)&cellModelType_d, sizeof(double)*nde);
	cudaMalloc((void**)&stimulation, sizeof(bool)*nde);
	cudaMalloc((void**)&ifscar_d, sizeof(bool)*nde);

	for(i=0; i<nde; i++) Istim[i] = false;
	cudaMemcpy(stimulation, Istim, sizeof(bool)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(ifscar_d, ifscar, sizeof(bool)*nde, cudaMemcpyHostToDevice);


	printf("...");
	CRN_initialValue();

	cudaMemcpy(volt_d, volt, sizeof(double)*nde, cudaMemcpyHostToDevice);

	cudaMemcpy(yglob0, yglob[0], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob1, yglob[1], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob2, yglob[2], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob3, yglob[3], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob4, yglob[4], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob5, yglob[5], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob6, yglob[6], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob7, yglob[7], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob8, yglob[8], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob9, yglob[9], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob10, yglob[10], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob11, yglob[11], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob12, yglob[12], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob13, yglob[13], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob14, yglob[14], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob15, yglob[15], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob16, yglob[16], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob17, yglob[17], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob18, yglob[18], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob19, yglob[19], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob20, yglob[20], sizeof(double)*nde, cudaMemcpyHostToDevice);

	free(ifscar);

	// ------------------------ Making CRN Lookup Table  Start ------------------------
	double *expdt_tau_m, *m_inf, *expdt_tau_h, *h_inf, *expdt_tau_j, *j_inf, *expdt_tau_xr, *xr_infinity, *expdt_tau_w, *w_infinity, *expdt_tau_xs, *xs_infinity, *expdt_tau_d, *d_infinity, *expdt_tau_oa, *oa_infinity;
	double *expdt_tau_oi, *oi_infinity, *expdt_tau_ua, *ua_infinity, *expdt_tau_ui, *ui_infinity, *expdt_tau_f, *f_infinity, *j_NaK, *j_NaCa1, *j_NaCa2, *g_Kur, *j_K1, *j_Kr;
	double *expdt_tau_fun, *fun_infinity, *expdt_tau_dt, *dt_infinity, *expdt_tau_ft, *ft_infinity, *d_infinity_new;
	cudaMalloc((void**)&expdt_tau_m, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&m_inf, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&expdt_tau_h, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&h_inf, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&expdt_tau_j, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&j_inf, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&expdt_tau_xr, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&xr_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&expdt_tau_w, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&w_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&expdt_tau_xs, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&xs_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&expdt_tau_d, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&d_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&expdt_tau_oa, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&oa_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&expdt_tau_oi, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&oi_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&expdt_tau_ua, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&ua_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&expdt_tau_ui, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&ui_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&expdt_tau_f, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&f_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&j_NaK, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&j_NaCa1, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&j_NaCa2, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&g_Kur, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&j_K1, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&j_Kr, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&expdt_tau_fun, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&fun_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&expdt_tau_dt, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&dt_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&expdt_tau_ft, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&ft_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&d_infinity_new, sizeof(double)*CRNTableVolt);

	setCRNTableVolt<<<grid, threads>>>(expdt_tau_m, m_inf, expdt_tau_h, h_inf, expdt_tau_j, j_inf, expdt_tau_xr, xr_infinity, expdt_tau_w, w_infinity, expdt_tau_xs, xs_infinity, expdt_tau_d, d_infinity, expdt_tau_oa, oa_infinity,
		expdt_tau_oi, oi_infinity, expdt_tau_ua, ua_infinity, expdt_tau_ui, ui_infinity, expdt_tau_f, f_infinity, j_NaK, j_NaCa1, j_NaCa2, g_Kur, j_K1, j_Kr,
		expdt_tau_fun, fun_infinity, expdt_tau_dt, dt_infinity, expdt_tau_ft, ft_infinity, d_infinity_new);
	cudaThreadSynchronize();
	// ------------------------ Making CRN Lookup Table  End --------------------------

	// Device Variables for diffKernel
	int *degree_d, *linkNode_d;
	double *weight_d;
	cudaMalloc((void**)&degree_d, sizeof(int)*(nde+1));
	cudaMalloc((void**)&linkNode_d, sizeof(int)*degree[nde]);
	cudaMalloc((void**)&weight_d, sizeof(double)*degree[nde]);
	cudaMalloc((void**)&volt_d_, sizeof(double)*nde);
	cudaMemcpy(degree_d, degree, sizeof(int)*(nde+1), cudaMemcpyHostToDevice);
	cudaMemcpy(linkNode_d, linkNode, sizeof(int)*degree[nde], cudaMemcpyHostToDevice);
	cudaMemcpy(weight_d, weightEdge, sizeof(double)*degree[nde], cudaMemcpyHostToDevice);

	
	// Initialize cells
	CRN_initialValue();


	// CUDA memcpy
	cudaMemcpy(volt_d, volt, sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(volt_d_, volt, sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(couplingType_d, couplingType, sizeof(char)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(cellModelType_d, cellModelType, sizeof(double)*nde, cudaMemcpyHostToDevice);

	cudaMemcpy(yglob0, yglob[0], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob1, yglob[1], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob2, yglob[2], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob3, yglob[3], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob4, yglob[4], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob5, yglob[5], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob6, yglob[6], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob7, yglob[7], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob8, yglob[8], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob9, yglob[9], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob10, yglob[10], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob11, yglob[11], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob12, yglob[12], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob13, yglob[13], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob14, yglob[14], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob15, yglob[15], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob16, yglob[16], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob17, yglob[17], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob18, yglob[18], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob19, yglob[19], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob20, yglob[20], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob21, yglob[21], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob22, yglob[22], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob23, yglob[23], sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob24, yglob[24], sizeof(double)*nde, cudaMemcpyHostToDevice);

	ttt = ttt0;
	
	printf("\nStart!!\n\n");
	printf("Iter       Time(ms)    Vmax(mV)   Vmax node\n");
	for(iter=1; iter<=itmax; iter++)
	{
		ttt = dt*iter;

		CRNKernel<<<grid, threads>>>(yglob0, yglob1, yglob2, yglob3, yglob4, yglob5, yglob6, yglob7, yglob8, yglob9, yglob10, yglob11, yglob12, yglob13, yglob14, yglob15, yglob16, yglob17, yglob18, yglob19, yglob20, yglob21, yglob22, yglob23, yglob24, volt_d, nde, cellModelType_d,
			expdt_tau_m, m_inf, expdt_tau_h, h_inf, expdt_tau_j, j_inf, expdt_tau_xr, xr_infinity, expdt_tau_w, w_infinity, expdt_tau_xs, xs_infinity, expdt_tau_d, d_infinity, expdt_tau_oa, oa_infinity,
			expdt_tau_oi, oi_infinity, expdt_tau_ua, ua_infinity, expdt_tau_ui, ui_infinity, expdt_tau_f, f_infinity, j_NaK, j_NaCa1, j_NaCa2, g_Kur, j_K1, j_Kr,
			expdt_tau_fun, fun_infinity, expdt_tau_dt, dt_infinity, expdt_tau_ft, ft_infinity, d_infinity_new);

		diffKernel<<<grid, threads>>>(volt_d_, volt_d, ifscar_d, cellModelType_d, linkNode_d, degree_d, weight_d, nde);

		diffKernel<<<grid, threads>>>(volt_d, volt_d_, ifscar_d, cellModelType_d, linkNode_d, degree_d, weight_d, nde);


		if(iter%10 == 0)
		{
			ttt = (double)((int)(ttt*2.0 + 0.25))/2.0; // Fix precision error

			cudaMemcpy(volt, volt_d, sizeof(double)*nde, cudaMemcpyDeviceToHost);

			ivmaxn = 0;
			vmax = -1000.0;
			for(i=0; i<nde; i++)
			{
				if(vmax < volt[i])
				{
					vmax = volt[i];
					ivmaxn = i;
				}
			}
			printf("%10d %10.1lf %10.3lf %10d\n", iter, ttt, vmax, ivmaxn);

			ttt_printVoltage = ttt;
			printVoltage();
		}

		if(iter%500 == 0)
		{
			
			WaitForSingleObject(printTecplotHandle, INFINITE);
			Sleep(5);
			CloseHandle(printTecplotHandle);
			
			//else PRINT_TECPLOT_FIRST = false;

			ttt_printTecplot = ttt;
			cudaMemcpy(volt_tecplot, volt_d, sizeof(double)*nde, cudaMemcpyDeviceToHost);
			printTecplotHandle = (HANDLE)_beginthreadex(NULL, 0, printTecplot, NULL, 0, (unsigned*)&printTecplotId);
		}
	}

	/////////////////// Deallocating Variables ///////////////////
	Sleep(30000);

	fclose(printVoltageFile);

	cudaFree(volt_d);
	cudaFree(stimulation);
	cudaFree(ifscar_d);

	cudaFree(yglob0);
	cudaFree(yglob1);
	cudaFree(yglob2);
	cudaFree(yglob3);
	cudaFree(yglob4);
	cudaFree(yglob5);
	cudaFree(yglob6);
	cudaFree(yglob7);
	cudaFree(yglob8);
	cudaFree(yglob9);
	cudaFree(yglob10);
	cudaFree(yglob11);
	cudaFree(yglob12);
	cudaFree(yglob13);
	cudaFree(yglob14);
	cudaFree(yglob15);
	cudaFree(yglob16);
	cudaFree(yglob17);
	cudaFree(yglob18);
	cudaFree(yglob19);
	cudaFree(yglob20);


	for(i=0; i<NUM_CRN_VAR; i++) free(yglob[i]);
	free(volt);
	free(volt_tecplot);
	free(Istim);
	free(init_stim);
	free(ifablation);
	for(i=0; i<3; i++) free(xy[i]);
	for(i=0; i<3; i++) free(ie[i]);

	printf("\n------------ END ------------\n");
	system("pause");





	return 0;
}