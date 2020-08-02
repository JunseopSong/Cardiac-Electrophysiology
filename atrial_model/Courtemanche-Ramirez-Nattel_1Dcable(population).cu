/*
Generate a population of 1D cable models
This code compute APD restitution curve for each ion channel parameter combination
Jun-Seop Song
*/


// ------------------------------ PROGRAM SETTING START ------------------------------
#define NUM_PARALLEL 64 // number of cables
#define N_ 128 // number of cells in each cable

#define NUM_PARAMETER 6 // number of remodeling parameters

#define Dfu 0.00154*0.3 // cm^2/ms
#define dx_ 0.025 // cm
// ------------------------------ PROGRAM SETTING END --------------------------------


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


int nde = NUM_PARALLEL * N_;
#define dt 0.05
#define dtinv 20.0 // dtinv = 1.0 / dt
#define itmax 360000 // time limit
#define tol 0.000001

double pacingStartTime[200];
int PACING_NUM;


////////////////////// CRN Constants Start //////////////////////
#define NUM_CRN_VAR 21

#define g_Na  7.8   // nanoS_per_picoF (in fast_sodium_current)
#define g_K1  0.09   // nanoS_per_picoF (in time_independent_potassium_current)
#define g_to  0.1652   // nanoS_per_picoF (in transient_outward_K_current)
#define g_Kr  0.029411765   // nanoS_per_picoF (in rapid_delayed_rectifier_K_current)
#define g_Ca_L  0.12375   // nanoS_per_picoF (in L_type_Ca_channel)
#define GKur 1.0 // dimensionless

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
#define stim_amplitude  -2900.0   // picoA (in membrane)
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
#define CRNTableVolt 16384
#define vlo (-90.0)
#define dvt 0.01

#define expdt_tau_u_1 0.993769491
#define expdt_tau_u_10 0.998750781
#define expdt_f_Ca_1 0.975309912
#define expdt_f_Ca_10 0.995012479
// --------------------------------------------------

int nstim, *init_stim;
bool *Istim;
double *yglob[NUM_CRN_VAR], yinitial[NUM_CRN_VAR], *volt, ttt;

FILE *databaseFile;
double *parameter;
int NUM_REPEAT;



/*
You can alternatively use the Sobol sequence.
*/
void parameterSetting()
{
	int i;

	for(i=0; i<NUM_PARALLEL; i++)
	{
		parameter[NUM_PARAMETER*i + 0] = 0.50 + (double)(rand()%101)/100.0; // g_Na: 0.50 ~ 1.50
		parameter[NUM_PARAMETER*i + 1] = 1.00 + (double)(rand()%151)/100.0; // g_K1: 1.00 ~ 2.50
		parameter[NUM_PARAMETER*i + 2] = 0.01 + (double)(rand()%100)/100.0; // g_to: 0.01 ~ 1.00
		parameter[NUM_PARAMETER*i + 3] = 0.50 + (double)(rand()%151)/100.0; // g_Kr: 0.50 ~ 2.00
		parameter[NUM_PARAMETER*i + 4] = 0.01 + (double)(rand()%100)/100.0; // g_CaL: 0.01 ~ 1.00
		parameter[NUM_PARAMETER*i + 5] = 0.01 + (double)(rand()%100)/100.0; // g_Kur: 0.01 ~ 1.00
	}
}


/*
Pacing protocol for computing APD restitution curve (steady-state)
*/
void pacingSetting()
{
	int i;

	pacingStartTime[0] = 1000.0;
	for(i=1; i<8; i++) pacingStartTime[i] = pacingStartTime[i-1] + 600.0;
	for(i=8*1; i<8*2; i++) pacingStartTime[i] = pacingStartTime[i-1] + 500.0;
	for(i=8*2; i<8*3; i++) pacingStartTime[i] = pacingStartTime[i-1] + 400.0;
	for(i=8*3; i<8*4; i++) pacingStartTime[i] = pacingStartTime[i-1] + 300.0;
	for(i=8*4; i<8*5; i++) pacingStartTime[i] = pacingStartTime[i-1] + 250.0;
	PACING_NUM = 8*5;
}


void preprocess()
{
	int i, j;

	for(i=0; i<NUM_CRN_VAR; i++) yglob[i] = (double*)malloc(sizeof(double)*nde);
	volt = (double*)malloc(sizeof(double)*nde);
	Istim = (bool*)malloc(sizeof(bool)*nde);
	parameter = (double*)malloc(sizeof(double)*NUM_PARALLEL*NUM_PARAMETER);

	nstim = 4*NUM_PARALLEL;
	init_stim = (int*)malloc(sizeof(int)*nstim);
	for(i=0; i<NUM_PARALLEL; i++)
	{
		for(j=0; j<4; j++) init_stim[i*4 + j] = i*N_ + j;
	}
}


void CRN_initialValue()
{
	int i, j;

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
		for(j=0; j<NUM_CRN_VAR; j++) yglob[j][i] = yinitial[j];
		volt[i] = yinitial[14];
	}
}


__global__ void setCRNTableVolt(double *tau_m, double *m_inf, double *tau_h, double *h_inf, double *tau_j, double *j_inf, double *tau_xr, double *xr_infinity,
	double *tau_w, double *w_infinity, double *tau_xs, double *xs_infinity, double *tau_d, double *d_infinity, double *tau_oa, double *oa_infinity,
	double *tau_oi, double *oi_infinity, double *tau_ua, double *ua_infinity, double *tau_ui, double *ui_infinity, double *tau_f, double *f_infinity, double *j_NaK,
	double *j_NaCa1, double *j_NaCa2, double *g_Kur, double *j_K1, double *j_Kr)
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
		tau_m[i] = (1.0/(alpha+beta));

		// h
		if (vv < -40.0) alpha = 0.135*exp((vv+80.0)/-6.8);
		else alpha = 0.0;
		if (vv < -40.0) beta = 3.56*exp(0.079*vv)+3.1e5*exp(0.35*vv);
		else beta = 1.0/(0.13*(1.0+exp((vv+10.66)/-11.1)));
		h_inf[i] = alpha/(alpha+beta);
		tau_h[i] = (1.0/(alpha+beta));

		// j
		if (vv < -40.0) alpha = (-1.2714e5*exp(0.2444*vv)-3.474e-5*exp(-0.04391*vv))*(vv+37.78)/(1.0+exp(0.311*(vv+79.23)));
		else alpha = 0.0;
		if (vv < -40.0) beta = 0.1212*exp(-0.01052*vv)/(1.0+exp(-0.1378*(vv+40.14)));
		else beta = 0.3*exp(-2.535e-7*vv)/(1.0+exp(-0.1*(vv+32.0)));
		j_inf[i] = alpha/(alpha+beta);
		tau_j[i] = (1.0/(alpha+beta));

		// xr
		if (fabs(vv+14.1) < tol) alpha = 0.0015;
		else alpha = 0.0003*(vv+14.1)/(1.0-exp((vv+14.1)/-5.0));
		if (fabs(vv-3.3328) < tol) beta = 3.7836118e-4;
		else beta = 0.000073898*(vv-3.3328)/(exp((vv-3.3328)/5.1237)-1.0);
		tau_xr[i] = (1/(alpha+beta));
		xr_infinity[i] = 1/(1.0+exp((vv+14.1)/-6.5));

		// w
		if (fabs(vv-7.9) < tol) tau_w[i] = (6.0*0.2/1.3);
		else tau_w[i] = (6.0*(1.0-exp(-(vv-7.9)/5.0))/((1.0+0.3*exp(-(vv-7.9)/5.0))*1.0*(vv-7.9)));
		w_infinity[i] = 1.0-1.0/(1.0+exp(-(vv-40.0)/17.0));

		// xs
		if (fabs(vv-19.9) < tol) alpha = 0.00068;
		else alpha = 0.00004*(vv-19.9)/(1.0-exp((vv-19.9)/-17.0));
		if (fabs(vv-19.9) < tol) beta = 0.000315;
		else beta = 0.000035*(vv-19.9)/(exp((vv-19.9)/9.0)-1.0);
		tau_xs[i] = (0.5*1/(alpha+beta));
		xs_infinity[i] = pow(1.0+exp((vv-19.9)/-12.7), -0.5);

		// d
		d_infinity[i] = 1.0/(1.0+exp((vv+10.0)/-8.0));
		if (fabs(vv+10.0) < tol) tau_d[i] = (4.579/(1.0+exp((vv+10.0)/-6.24)));
		else tau_d[i] = ((1.0-exp((vv+10.0)/-6.24))/(0.035*(vv+10.0)*(1.0+exp((vv+10.0)/-6.24))));

		// oa
		alpha = 0.65*1.0/(exp((vv-(-10.0))/-8.5)+exp((vv-(-10.0)-40.0)/-59.0));
		beta = 0.65*1.0/(2.5+exp((vv-(-10.0)+72.0)/17.0));
		tau_oa[i] = (1/(alpha+beta)/K_Q10);
		oa_infinity[i] = 1.0/(1.0+exp((vv-(-10.0)+10.47)/-17.54));

		// oi
		alpha = 1.0/(18.53+1.0*exp((vv-(-10.0)+103.7)/10.95));
		beta = 1.0/(35.56+1.0*exp((vv-(-10.0)-8.74)/-7.44));
		tau_oi[i] = (1/(alpha+beta)/K_Q10);
		oi_infinity[i] = 1.0/(1.0+exp((vv-(-10.0)+33.1)/5.3));

		// ua
		alpha = 0.65*1.0/(exp((vv-(-10.0))/-8.5)+exp((vv-(-10.0)-40.0)/-59.0));
		beta = 0.65*1.0/(2.5+exp((vv-(-10.0)+72.0)/17.0));
		tau_ua[i] = (1/(alpha+beta)/K_Q10);
		ua_infinity[i] = 1.0/(1.0+exp((vv-(-10.0)+20.3)/-9.6));

		// ui
		alpha = 1.0/(21.0+1.0*exp((vv-(-10.0)-195.0)/-28.0));
		beta = 1.0/exp((vv-(-10.0)-168.0)/-16.0);
		tau_ui[i] = (1.0/(alpha+beta)/K_Q10);
		ui_infinity[i] = 1.0/(1.0+exp((vv-(-10.0)-109.45)/27.48));

		// f
		f_infinity[i] = exp(-(vv+28.0)/6.9)/(1.0+exp(-(vv+28.0)/6.9));
		tau_f[i] = (9.0*1.0/(0.0197*exp(-0.0337*0.0337*(vv+10.0)*(vv+10.0))+0.02));

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
	}
}


__global__ void CRNKernel(double* u, double* vc, double* w, double* d, double*  f_Ca, double* f, double* h, double* jj, double* m, double* Ca_i, double* Ca_rel, double* Ca_up, double* K_i, double* Na_i, double* v, double* xr, double* xs, double* oa, double* oi, double* ua, double* ui,
	double* volt_ , bool *stim, const int num, double *remodelingParameter,
	const double *tau_m, const double *m_inf, const double *tau_h, const double *h_inf, const double *tau_j, const double *j_inf, const double *tau_xr, const double *xr_infinity,
	const double *tau_w, const double *w_infinity, const double *tau_xs, const double *xs_infinity, const double *tau_d, const double *d_infinity, const double *tau_oa, const double *oa_infinity,
	const double *tau_oi, const double *oi_infinity, const double *tau_ua, const double *ua_infinity, const double *tau_ui, const double *ui_infinity, const double *tau_f, const double *f_infinity, const double *j_NaK,
	const double *j_NaCa1, const double *j_NaCa2, const double *g_Kur, const double *j_K1, const double *j_Kr)
{
	int id = blockIdx.x*blockDim.x + threadIdx.x;
	int parameterIndex = NUM_PARAMETER*blockIdx.x;
	
	int i, repeat;
	double dtime;

	double qv1, qv2, temp1, temp2;
	int vindex1, vindex2;

	double Vm, Nai, Cai, Ki;
	double i_up_leak, i_rel, i_Ca_L, i_NaCa, Fn, i_up, E_Ca, E_Na, i_B_Na, i_B_Ca, E_K, i_B_K, i_Na;
	double i_NaK, i_K1, i_to, i_Kur, i_Kr, i_Ks, i_CaP, i_tr;


	Vm = volt_[id];
	Nai = Na_i[id];
	Cai = Ca_i[id];
	Ki = K_i[id];


	if(fabs((v[id]-Vm)/dt)>25.0)
	{
		dtime = dt / 5.0;
		repeat = 5;
	}
	else
	{
		dtime = dt;
		repeat = 1;
	}

	v[id] = Vm;

	for(i=0; i<repeat; i++)
	{
		// Lookup table index
		temp1 = (Vm-vlo)/dvt;
		vindex1 = (int)temp1;
		vindex2 = vindex1 + 1;
		qv1 = temp1 - (double)vindex1;
		qv2 = (double)vindex2 - temp1;

		i_up_leak = I_up_max*Ca_up[id]/Ca_up_max;

		i_rel = K_rel*u[id]*u[id]*vc[id]*w[id]*(Ca_rel[id]-Cai);

		i_Ca_L = remodelingParameter[parameterIndex+4]*Cm*g_Ca_L*d[id]*f[id]*f_Ca[id]*(Vm-65.0);

		temp1 = j_NaCa1[vindex1]*qv2 + j_NaCa1[vindex2]*qv1;
		temp2 = j_NaCa2[vindex1]*qv2 + j_NaCa2[vindex2]*qv1;
		i_NaCa = temp1*Nai*Nai*Nai - temp2*Cai;

		Fn = 1.0e-12*V_rel*i_rel-5.0e-13/F*(0.5*i_Ca_L-0.2*i_NaCa);

		temp1 = 1.0/(1.0+__expf(-(Fn-3.4175e-13)/13.67e-16));
		if(repeat == 1) u[id] = temp1-(temp1-u[id])*expdt_tau_u_1;
		else u[id] = temp1-(temp1-u[id])*expdt_tau_u_10;

		temp1 = 1.0-1.0/(1.0+__expf(-(Fn-6.835e-14)/13.67e-16));
		temp2 = 1.91+2.09*1.0/(1.0+__expf(-(Fn-3.4175e-13)/13.67e-16));
		vc[id] = temp1-(temp1-vc[id])*__expf(-dtime/temp2);

		temp1 = w_infinity[vindex1]*qv2 + w_infinity[vindex2]*qv1;
		temp2 = tau_w[vindex1]*qv2 + tau_w[vindex2]*qv1;
		w[id] = temp1-(temp1-w[id])*__expf(-dtime/temp2);

		i_up = I_up_max/(1.0+K_up/Cai);


		temp1 = d_infinity[vindex1]*qv2 + d_infinity[vindex2]*qv1;
		temp2 = tau_d[vindex1]*qv2 + tau_d[vindex2]*qv1;
		d[id] = temp1-(temp1-d[id])*__expf(-dtime/temp2);

		temp1 = 1.0/(1.0+Cai/0.00035);
		if(repeat == 1) f_Ca[id] = temp1-(temp1-f_Ca[id])*expdt_f_Ca_1;
		else f_Ca[id] = temp1-(temp1-f_Ca[id])*expdt_f_Ca_10;

		temp1 = f_infinity[vindex1]*qv2 + f_infinity[vindex2]*qv1;
		temp2 = tau_f[vindex1]*qv2 + tau_f[vindex2]*qv1;
		f[id] = temp1-(temp1-f[id])*__expf(-dtime/temp2);

		E_Ca = RT_F/2.0*__logf(Ca_o/Cai);
		E_Na = RT_F*__logf(Na_o/Nai);
		i_B_Na = Cm*g_B_Na*(Vm-E_Na);
		i_B_Ca = Cm*g_B_Ca*(Vm-E_Ca);
		E_K = RT_F*__logf(K_o/Ki);
		i_B_K = Cm*g_B_K*(Vm-E_K);

		i_Na = remodelingParameter[parameterIndex]*Cm*g_Na*m[id]*m[id]*m[id]*h[id]*jj[id]*(Vm-E_Na);


		temp1 = h_inf[vindex1]*qv2 + h_inf[vindex2]*qv1;
		temp2 = tau_h[vindex1]*qv2 + tau_h[vindex2]*qv1;
		h[id] = temp1-(temp1-h[id])*__expf(-dtime/temp2);

		temp1 = j_inf[vindex1]*qv2 + j_inf[vindex2]*qv1;
		temp2 = tau_j[vindex1]*qv2 + tau_j[vindex2]*qv1;
		jj[id] = temp1-(temp1-jj[id])*__expf(-dtime/temp2);

		temp1 = m_inf[vindex1]*qv2 + m_inf[vindex2]*qv1;
		temp2 = tau_m[vindex1]*qv2 + tau_m[vindex2]*qv1;
		m[id] = temp1-(temp1-m[id])*__expf(-dtime/temp2);

		temp1 = j_NaK[vindex1]*qv2 + j_NaK[vindex2]*qv1;
		temp2 = Km_Na_i/Nai;
		i_NaK = temp1 / (1.0+temp2*sqrt(temp2));

		temp1 = j_K1[vindex1]*qv2 + j_K1[vindex2]*qv1;
		i_K1 = remodelingParameter[parameterIndex+1]*temp1*(Vm-E_K);

		i_to = remodelingParameter[parameterIndex+2]*Cm*g_to*oa[id]*oa[id]*oa[id]*oi[id]*(Vm-E_K);

		temp1 = g_Kur[vindex1]*qv2 + g_Kur[vindex2]*qv1;
		i_Kur = remodelingParameter[parameterIndex+5]*temp1*ua[id]*ua[id]*ua[id]*ui[id]*(Vm-E_K);

		temp1 = j_Kr[vindex1]*qv2 + j_Kr[vindex2]*qv1;
		i_Kr = remodelingParameter[parameterIndex+3]*temp1*xr[id]*(Vm-E_K);

		i_Ks = Cm*g_Ks*xs[id]*xs[id]*(Vm-E_K);
		i_CaP = Cm*i_CaP_max*Cai/(0.0005+Cai);
		i_tr = (Ca_up[id]-Ca_rel[id])/tau_tr;


		temp1 = xr_infinity[vindex1]*qv2 + xr_infinity[vindex2]*qv1;
		temp2 = tau_xr[vindex1]*qv2 + tau_xr[vindex2]*qv1;
		xr[id] = temp1-(temp1-xr[id])*__expf(-dtime/temp2);

		temp1 = xs_infinity[vindex1]*qv2 + xs_infinity[vindex2]*qv1;
		temp2 = tau_xs[vindex1]*qv2 + tau_xs[vindex2]*qv1;
		xs[id] = temp1-(temp1-xs[id])*__expf(-dtime/temp2);

		temp1 = oa_infinity[vindex1]*qv2 + oa_infinity[vindex2]*qv1;
		temp2 = tau_oa[vindex1]*qv2 + tau_oa[vindex2]*qv1;
		oa[id] = temp1-(temp1-oa[id])*__expf(-dtime/temp2);

		temp1 = oi_infinity[vindex1]*qv2 + oi_infinity[vindex2]*qv1;
		temp2 = tau_oi[vindex1]*qv2 + tau_oi[vindex2]*qv1;
		oi[id] = temp1-(temp1-oi[id])*__expf(-dtime/temp2);

		temp1 = ua_infinity[vindex1]*qv2 + ua_infinity[vindex2]*qv1;
		temp2 = tau_ua[vindex1]*qv2 + tau_ua[vindex2]*qv1;
		ua[id] = temp1-(temp1-ua[id])*__expf(-dtime/temp2);

		temp1 = ui_infinity[vindex1]*qv2 + ui_infinity[vindex2]*qv1;
		temp2 = tau_ui[vindex1]*qv2 + tau_ui[vindex2]*qv1;
		ui[id] = temp1-(temp1-ui[id])*__expf(-dtime/temp2);


		temp1 = (2.0*i_NaCa-(i_CaP+i_Ca_L+i_B_Ca))/(2.0*V_i*F)+(V_up*(i_up_leak-i_up)+i_rel*V_rel)/V_i; // B1
		temp2 = 1.0+TRPN_max*Km_TRPN/((Cai+Km_TRPN)*(Cai+Km_TRPN))+CMDN_max*Km_CMDN/((Cai+Km_CMDN)*(Cai+Km_CMDN)); // B2


		if (stim[id] && Vm < 0.0)
		{
			Vm -= (stim_amplitude + i_Na+i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_Na+i_B_Ca+i_NaK+i_CaP+i_NaCa+i_Ca_L)/Cm*dtime;

			Cai += (temp1/temp2)*dtime;

			Ca_up[id] += (i_up-(i_up_leak+i_tr*V_rel/V_up))*dtime;

			Ca_rel[id] += ((i_tr-i_rel)*1.0/(1.0+CSQN_max*Km_CSQN/((Ca_rel[id]+Km_CSQN)*(Ca_rel[id]+Km_CSQN))))*dtime;

			Nai += (-3.0*i_NaK-(3.0*i_NaCa+i_B_Na+i_Na))/(V_i*F)*dtime;

			Ki += (2.0*i_NaK-(i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_K))/(V_i*F)*dtime;

		}
		else
		{
			Vm -= (i_Na+i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_Na+i_B_Ca+i_NaK+i_CaP+i_NaCa+i_Ca_L)/Cm*dtime;

			Cai += (temp1/temp2)*dtime;

			Ca_up[id] += (i_up-(i_up_leak+i_tr*V_rel/V_up))*dtime;

			Ca_rel[id] += ((i_tr-i_rel)*1.0/(1.0+CSQN_max*Km_CSQN/((Ca_rel[id]+Km_CSQN)*(Ca_rel[id]+Km_CSQN))))*dtime;

			Nai += (-3.0*i_NaK-(3.0*i_NaCa+i_B_Na+i_Na))/(V_i*F)*dtime;

			Ki += (2.0*i_NaK-(i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_K))/(V_i*F)*dtime;
		}
	}

	volt_[id] = Vm;
	Na_i[id] = Nai;
	Ca_i[id] = Cai;
	K_i[id] = Ki;
}


__global__ void diffKernel(double *volt_after, double *volt_before, const int num)
{
	int id = blockIdx.x*blockDim.x + threadIdx.x;

	double Dfudtdx2 = Dfu*dt/(dx_*dx_)/2.0;

	if(id%N_ == 0) volt_after[id] = volt_before[id] + (volt_before[id+1]-volt_before[id])*Dfudtdx2;
	else if(id%N_ == N_-1) volt_after[id] = volt_before[id] + (volt_before[id-1]-volt_before[id])*Dfudtdx2;
	else volt_after[id] = volt_before[id] + (volt_before[id-1]+volt_before[id+1]-2.0*volt_before[id])*Dfudtdx2;
}


__global__ void setAPDCalculator(bool *detectStart, bool *detectPeak, bool *detectEnd, int *apdCnt, int *apd, double *preVolt)
{
	int id = threadIdx.x;
	int i;

	detectStart[id] = true;
	detectPeak[id] = false;
	detectEnd[id] = false;
	apdCnt[id] = 0;

	preVolt[id] = 0.0;

	for(i=0; i<5; i++) apd[id*5 + i] = 0;
}


__global__ void apd90Kernel(bool *detectStart, bool *detectPeak, bool *detectEnd, int *apdCnt, int *apd, double *preVolt, double *startTime, double *restPotential, double *v90, double *volt_, double time_)
{
	int id = threadIdx.x;
	int voltage = volt_[id*N_ + N_/2];

	if(detectStart[id])
	{
		if((voltage-preVolt[id])/dt > 0.2)
		{
			detectStart[id] = false;
			detectPeak[id] = true;
			startTime[id] = time_;
			restPotential[id] = preVolt[id];
		}
	}
	
	if(detectPeak[id])
	{
		if(voltage > -10.0 && voltage <= preVolt[id])
		{
			detectPeak[id] = false;
			detectEnd[id] = true;
			v90[id] = voltage - (voltage - restPotential[id])*0.9;
		}
	}

	if(detectEnd[id])
	{
		if(voltage <= v90[id])
		{
			detectEnd[id] = false;
			detectStart[id] = true;
			apdCnt[id]++;

			if(apdCnt[id]%8 == 0)
			{
				apd[id*5 + (int)((apdCnt[id]+0.1)/8) - 1] = (int)(time_ - startTime[id] + 0.5);
			}
		}
	}

	preVolt[id] = voltage;
}


__global__ void apd74Kernel(bool *detectStart, bool *detectPeak, bool *detectEnd, int *apdCnt, int *apd, double *preVolt, double *startTime, double *volt_, double time_)
{
	int id = threadIdx.x;
	int voltage = volt_[id*N_ + N_/2];

	if(detectStart[id])
	{
		if((voltage-preVolt[id])/dt > 0.2)
		{
			detectStart[id] = false;
			detectPeak[id] = true;
			startTime[id] = time_;
		}
	}

	if(detectPeak[id])
	{
		if(voltage > -30.0)
		{
			detectPeak[id] = false;
			detectEnd[id] = true;
		}
	}

	if(detectEnd[id])
	{
		if(voltage <= -74.0)
		{
			detectEnd[id] = false;
			detectStart[id] = true;
			apdCnt[id]++;

			if(apdCnt[id]%8 == 0)
			{
				apd[id*5 + (int)((apdCnt[id]+0.1)/8) - 1] = (int)(time_ - startTime[id] + 0.5);
			}
		}
	}

	preVolt[id] = voltage;
}



int main()
{
	int i, j, iter, i_rept, pacingCnt = 0, apd_host[NUM_PARALLEL*5];
	bool pacingFinish = false, stimEnd = true;
	double tttt;

	printf("Input number of repetition: ");
	scanf("%d", &NUM_REPEAT);

	databaseFile = fopen("Database.txt", "a"); // output file

	pacingSetting();

	///////////////// CUDA Threads, Blocks Setting Start ///////////////
	dim3 grid(NUM_PARALLEL, 1, 1);
	dim3 threads(N_, 1, 1);
	///////////////// CUDA Threads, Blocks Setting End /////////////////

	preprocess();

	srand((unsigned)time(NULL));


	// ------------------------ Making CRN Lookup Table  Start ------------------------
	double *tau_m, *m_inf, *tau_h, *h_inf, *tau_j, *j_inf, *tau_xr, *xr_infinity, *tau_w, *w_infinity, *tau_xs, *xs_infinity, *tau_d, *d_infinity, *tau_oa, *oa_infinity;
	double *tau_oi, *oi_infinity, *tau_ua, *ua_infinity, *tau_ui, *ui_infinity, *tau_f, *f_infinity, *j_NaK, *j_NaCa1, *j_NaCa2, *g_Kur, *j_K1, *j_Kr;
	cudaMalloc((void**)&tau_m, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&m_inf, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&tau_h, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&h_inf, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&tau_j, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&j_inf, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&tau_xr, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&xr_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&tau_w, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&w_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&tau_xs, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&xs_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&tau_d, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&d_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&tau_oa, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&oa_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&tau_oi, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&oi_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&tau_ua, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&ua_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&tau_ui, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&ui_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&tau_f, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&f_infinity, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&j_NaK, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&j_NaCa1, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&j_NaCa2, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&g_Kur, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&j_K1, sizeof(double)*CRNTableVolt);
	cudaMalloc((void**)&j_Kr, sizeof(double)*CRNTableVolt);

	setCRNTableVolt<<<grid, threads>>>(tau_m, m_inf, tau_h, h_inf, tau_j, j_inf, tau_xr, xr_infinity, tau_w, w_infinity, tau_xs, xs_infinity, tau_d, d_infinity, tau_oa, oa_infinity,
		tau_oi, oi_infinity, tau_ua, ua_infinity, tau_ui, ui_infinity, tau_f, f_infinity, j_NaK, j_NaCa1, j_NaCa2, g_Kur, j_K1, j_Kr);
	cudaThreadSynchronize();
	// ------------------------ Making CRN Lookup Table  End --------------------------



	// ============================ Device Variables ============================
	double *parameter_d;
	double *volt_d, *volt_d_;
	bool *stimulation;
	double *yglob0, *yglob1, *yglob2, *yglob3, *yglob4, *yglob5, *yglob6, *yglob7, *yglob8, *yglob9, *yglob10, *yglob11, *yglob12, *yglob13, *yglob14, *yglob15, *yglob16, *yglob17, *yglob18, *yglob19, *yglob20;

	bool *detectStart, *detectPeak, *detectEnd;
	int *apdCnt, *apd;
	double *preVolt, *startTime, *restPotential, *v90;

	cudaMalloc((void**)&detectStart, sizeof(bool)*NUM_PARALLEL);
	cudaMalloc((void**)&detectPeak, sizeof(bool)*NUM_PARALLEL);
	cudaMalloc((void**)&detectEnd, sizeof(bool)*NUM_PARALLEL);
	cudaMalloc((void**)&apdCnt, sizeof(int)*NUM_PARALLEL);
	cudaMalloc((void**)&apd, sizeof(int)*NUM_PARALLEL*5);
	cudaMalloc((void**)&preVolt, sizeof(double)*NUM_PARALLEL);
	cudaMalloc((void**)&startTime, sizeof(double)*NUM_PARALLEL);
	cudaMalloc((void**)&restPotential, sizeof(double)*NUM_PARALLEL);
	cudaMalloc((void**)&v90, sizeof(double)*NUM_PARALLEL);

	cudaMalloc((void**)&parameter_d, sizeof(double)*NUM_PARALLEL*NUM_PARAMETER);

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

	cudaMalloc((void**)&volt_d, sizeof(double)*nde);
	cudaMalloc((void**)&volt_d_, sizeof(double)*nde);
	cudaMalloc((void**)&stimulation, sizeof(bool)*nde);


	// Initialize cells
	CRN_initialValue();


	printf("\nStart!\n\n");
	for(i_rept = 0; i_rept < NUM_REPEAT; i_rept++)
	{
		parameterSetting();
		cudaMemcpy(parameter_d, parameter, sizeof(double)*NUM_PARALLEL*NUM_PARAMETER, cudaMemcpyHostToDevice);

		// Initialize variables
		pacingCnt = 0;
		pacingFinish = false;
		stimEnd = true;
		
		setAPDCalculator<<<1, NUM_PARALLEL>>>(detectStart, detectPeak, detectEnd, apdCnt, apd, preVolt);


		// CUDA memcpy
		cudaMemcpy(volt_d, volt, sizeof(double)*nde, cudaMemcpyHostToDevice);
		cudaMemcpy(volt_d_, volt, sizeof(double)*nde, cudaMemcpyHostToDevice);

		for(i=0; i<nde; i++) Istim[i] = false;
		cudaMemcpy(stimulation, Istim, sizeof(bool)*nde, cudaMemcpyHostToDevice);

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


		for(iter=1; iter<itmax; iter++)
		{
			ttt = dt*iter;

			if(iter%5000 == 0) printf("iter: %7d      time (ms): %7.1lf\n", iter, ttt);

			if(!pacingFinish)
			{
				tttt = ttt - pacingStartTime[pacingCnt];

				if(stimEnd && tttt >= 0.0 && tttt <= 4.0)
				{
					for(i=0; i<nstim; i++) Istim[init_stim[i]] = true;
					stimEnd = false;
					cudaMemcpy(stimulation, Istim, sizeof(bool)*nde, cudaMemcpyHostToDevice);
				}

				if(!stimEnd && tttt > 4.0)
				{
					for(i=0; i<nde; i++) Istim[i] = false;
					stimEnd = true;
					pacingCnt++;
					cudaMemcpy(stimulation, Istim, sizeof(bool)*nde, cudaMemcpyHostToDevice);
				}

				if(pacingCnt == PACING_NUM)
				{
					pacingFinish = true;
				}
			}


			CRNKernel<<<grid, threads>>>(yglob0, yglob1, yglob2, yglob3, yglob4, yglob5, yglob6, yglob7, yglob8, yglob9, yglob10, yglob11, yglob12, yglob13, yglob14, yglob15, yglob16, yglob17, yglob18, yglob19, yglob20, volt_d, stimulation, nde, parameter_d,
				tau_m, m_inf, tau_h, h_inf, tau_j, j_inf, tau_xr, xr_infinity, tau_w, w_infinity, tau_xs, xs_infinity, tau_d, d_infinity, tau_oa, oa_infinity,
				tau_oi, oi_infinity, tau_ua, ua_infinity, tau_ui, ui_infinity, tau_f, f_infinity, j_NaK, j_NaCa1, j_NaCa2, g_Kur, j_K1, j_Kr);

			diffKernel<<<grid, threads>>>(volt_d_, volt_d, nde);

			diffKernel<<<grid, threads>>>(volt_d, volt_d_, nde);


			if(iter%10 == 0)
			{
				//apd90Kernel<<<1, NUM_PARALLEL>>>(detectStart, detectPeak, detectEnd, apdCnt, apd, preVolt, startTime, restPotential, v90, volt_d, ttt);
				apd74Kernel<<<1, NUM_PARALLEL>>>(detectStart, detectPeak, detectEnd, apdCnt, apd, preVolt, startTime, volt_d, ttt);
			}
		}


		// Print parameters and APDs.
		cudaMemcpy(apd_host, apd, sizeof(int)*NUM_PARALLEL*5, cudaMemcpyDeviceToHost);

		for(i=0; i<NUM_PARALLEL; i++)
		{
			for(j=0; j<NUM_PARAMETER; j++)
			{
				fprintf(databaseFile, "%.2lf  ", parameter[i*NUM_PARAMETER + j]);
				printf("%.2lf  ", parameter[i*NUM_PARAMETER + j]);
			}

			for(j=0; j<5; j++)
			{
				fprintf(databaseFile, "%d  ", apd_host[i*5 + j]);
				printf("%d  ", apd_host[i*5 + j]);
			}

			fprintf(databaseFile, "\n");
			printf("\n");
		}
	}

	fclose(databaseFile);

	return 0;
}