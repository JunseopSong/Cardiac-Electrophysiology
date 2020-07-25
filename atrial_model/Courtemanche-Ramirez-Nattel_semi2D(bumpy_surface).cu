/*
Atrial fibrillation simulation on semi-2D surface (bumpy geometry)
Reference:
Jun-Seop Song et al. (2018) Pro-Arrhythmogenic Effects of Heterogeneous Tissue Curvature - Suggestion for Role of Left Atrial Appendage in Atrial Fibrillation. Circulation Journal, 83(1), 32-40.
https://doi.org/10.1253/circj.CJ-18-0615
*/


// ------------------------------ PROGRAM SETTING START ------------------------------
#define g_Na  7.8   // nanoS_per_picoF (in fast_sodium_current)
#define g_K1  0.09*1.01   // nanoS_per_picoF (in time_independent_potassium_current)
#define g_to  0.1652   // nanoS_per_picoF (in transient_outward_K_current)
#define g_Kr  0.029411765*1.02   // nanoS_per_picoF (in rapid_delayed_rectifier_K_current)
#define g_Ca_L  0.12375*0.97   // nanoS_per_picoF (in L_type_Ca_channel)
#define GKur 1 // dimensionless
#define I_NaCa_max  1600.0   // picoA_per_picoF (in Na_Ca_exchanger_current)
#define I_up_max  0.005   // millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)

#define S1S2Protocol
#define S1S2_TARGET_VOLT -60.0

#define M 512
#define DiffCoeff 0.00097 // cm^2/ms
#define dx_ 0.025 // cm

// z = A*sin(2*PI/B*x)*cos(2*PI/B*y)
#define TWOPI 6.283185307179586
#define ALPHA 3.0 // electrical & bumpiness gradient
#define GEOMETRY_A 24.0
#define GEOMETRY_B (M/4.0)
// ------------------------------ PROGRAM SETTING END --------------------------------


#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>
#include <process.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


#define index(i,j) i*M+j

#define dt 0.02
#define itmax 252100 // time limit
#define tol 0.000001

double pacingStartTime[] = { 0 };
int PACING_NUM = sizeof(pacingStartTime) / sizeof(pacingStartTime[0]);
int nde = M*M;


////////////////////// CRN Constants Start //////////////////////
#define NUM_CRN_VAR 21

#define CMDN_max  0.05   // millimolar (in Ca_buffers)
#define CSQN_max  10.0   // millimolar (in Ca_buffers)
#define Km_CMDN  0.00238   // millimolar (in Ca_buffers)
#define Km_CSQN  0.8   // millimolar (in Ca_buffers)
#define Km_TRPN  0.0005   // millimolar (in Ca_buffers)
#define TRPN_max  0.07   // millimolar (in Ca_buffers)
#define Ca_up_max  15.0   // millimolar (in Ca_leak_current_by_the_NSR)
#define K_rel  30.0   // per_millisecond (in Ca_release_current_from_JSR)
#define K_up  0.00092   // millimolar (in Ca_uptake_current_by_the_NSR)
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
#define CRNTableVolt 24576
#define vlo (-90.0)
#define dvt 0.005

#define expdt_tau_u 0.9975031224
#define expdt_tau_f_Ca 0.9900498337
// --------------------------------------------------

int nstim, *init_stim;
bool *Istim, *ifscar;
double *yglob[NUM_CRN_VAR], yinitial[NUM_CRN_VAR], *volt, *volt_print, ttt, ttt_print;



void preprocess()
{
	int i;

	for (i = 0; i<NUM_CRN_VAR; i++) yglob[i] = (double*)malloc(sizeof(double)*nde);
	volt = (double*)malloc(sizeof(double)*nde);
	volt_print = (double*)malloc(sizeof(double)*nde);
	Istim = (bool*)malloc(sizeof(bool)*nde);
	ifscar = (bool*)malloc(sizeof(bool)*nde);

	for (i = 0; i<nde; i++)
	{
		Istim[i] = false;
		ifscar[i] = false;
	}

	nstim = M;
	init_stim = (int*)malloc(sizeof(int)*nstim);
	for (i = 0; i<M; i++)
	{
		init_stim[i] = i;
	}
}


void CRN_initialValue()
{
	int i, j;

	yinitial[0] = 0.0;  // u (dimensionless) (in Ca_release_current_from_JSR_u_gate)
	yinitial[1] = 1.0;        // v (dimensionless) (in Ca_release_current_from_JSR_v_gate)
	yinitial[2] = 0.9991;     // w (dimensionless) (in Ca_release_current_from_JSR_w_gate)
	yinitial[3] = 0.0002;   // d (dimensionless) (in L_type_Ca_channel_d_gate)
	yinitial[4] = 0.6687;   // f_Ca (dimensionless) (in L_type_Ca_channel_f_Ca_gate)
	yinitial[5] = 0.7741;   // f (dimensionless) (in L_type_Ca_channel_f_gate)
	yinitial[6] = 0.9461;   // h (dimensionless) (in fast_sodium_current_h_gate)
	yinitial[7] = 0.9589;   // j (dimensionless) (in fast_sodium_current_j_gate)
	yinitial[8] = 0.0040;   // m (dimensionless) (in fast_sodium_current_m_gate)
	yinitial[9] = 0.0002;   // Ca_i (millimolar) (in intracellular_ion_concentrations)
	yinitial[10] = 0.6842;     // Ca_rel (millimolar) (in intracellular_ion_concentrations)
	yinitial[11] = 1.7812;     // Ca_up (millimolar) (in intracellular_ion_concentrations)
	yinitial[12] = 138.6555;    // K_i (millimolar) (in intracellular_ion_concentrations)
	yinitial[13] = 11.3828;   // Na_i (millimolar) (in intracellular_ion_concentrations)
	yinitial[14] = -79.2567;    // V (millivolt) (in membrane)
	yinitial[15] = 0.0257;  // xr (dimensionless) (in rapid_delayed_rectifier_K_current_xr_gate)
	yinitial[16] = 0.0350;  // xs (dimensionless) (in slow_delayed_rectifier_K_current_xs_gate)
	yinitial[17] = 0.0339;  // oa (dimensionless) (in transient_outward_K_current_oa_gate)
	yinitial[18] = 0.9989;  // oi (dimensionless) (in transient_outward_K_current_oi_gate)
	yinitial[19] = 0.0061;  // ua (dimensionless) (in ultrarapid_delayed_rectifier_K_current_ua_gate)
	yinitial[20] = 0.9882;  // ui (dimensionless) (in ultrarapid_delayed_rectifier_K_current_ui_gate)

							/*
							for(j=0; j<NUM_CRN_VAR; j++) yglob[j][0] = yinitial[j];
							volt[0] = yinitial[14];

							// Make resting state
							for(i=0; i<500000; i++) CRNcompute(yglob[0], yglob[1], yglob[2], yglob[3], yglob[4], yglob[5], yglob[6], yglob[7], yglob[8], yglob[9], yglob[10], yglob[11], yglob[12], yglob[13], yglob[14], yglob[15], yglob[16], yglob[17], yglob[18], yglob[19], yglob[20], volt, Istim, 0);

							for(j=0; j<NUM_CRN_VAR; j++) yinitial[j] = yglob[j][0];
							*/

	for (i = 0; i<nde; i++)
	{
		for (j = 0; j<NUM_CRN_VAR; j++) yglob[j][i] = yinitial[j];
		volt[i] = yinitial[14];
	}
}


__global__ void setCRNTableVolt(double *expdt_tau_m, double *m_inf, double *expdt_tau_h, double *h_inf, double *expdt_tau_j, double *j_inf, double *expdt_tau_xr, double *xr_infinity,
	double *expdt_tau_w, double *w_infinity, double *expdt_tau_xs, double *xs_infinity, double *expdt_tau_d, double *d_infinity, double *expdt_tau_oa, double *oa_infinity,
	double *expdt_tau_oi, double *oi_infinity, double *expdt_tau_ua, double *ua_infinity, double *expdt_tau_ui, double *ui_infinity, double *expdt_tau_f, double *f_infinity, double *j_NaK,
	double *j_NaCa1, double *j_NaCa2, double *g_Kur, double *j_K1, double *j_Kr)
{
	double vv;

	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int bdim = blockDim.x;
	const unsigned int gdim = gridDim.x;
	int step = bdim*gdim;

	double alpha, beta;

	for (int i = bid*bdim + tid; i<CRNTableVolt; i += step)
	{
		vv = vlo + (double)i * dvt;

		// m
		if (fabs(vv + 47.13) < tol) alpha = 3.2;
		else alpha = 0.32*(vv + 47.13) / (1.0 - exp(-0.1*(vv + 47.13)));
		beta = 0.08*exp(-vv / 11.0);
		m_inf[i] = alpha / (alpha + beta);
		expdt_tau_m[i] = exp(-dt / (1.0 / (alpha + beta)));

		// h
		if (vv < -40.0) alpha = 0.135*exp((vv + 80.0) / -6.8);
		else alpha = 0.0;
		if (vv < -40.0) beta = 3.56*exp(0.079*vv) + 3.1e5*exp(0.35*vv);
		else beta = 1.0 / (0.13*(1.0 + exp((vv + 10.66) / -11.1)));
		h_inf[i] = alpha / (alpha + beta);
		expdt_tau_h[i] = exp(-dt / (1.0 / (alpha + beta)));

		// j
		if (vv < -40.0) alpha = (-1.2714e5*exp(0.2444*vv) - 3.474e-5*exp(-0.04391*vv))*(vv + 37.78) / (1.0 + exp(0.311*(vv + 79.23)));
		else alpha = 0.0;
		if (vv < -40.0) beta = 0.1212*exp(-0.01052*vv) / (1.0 + exp(-0.1378*(vv + 40.14)));
		else beta = 0.3*exp(-2.535e-7*vv) / (1.0 + exp(-0.1*(vv + 32.0)));
		j_inf[i] = alpha / (alpha + beta);
		expdt_tau_j[i] = exp(-dt / (1.0 / (alpha + beta)));

		// xr
		if (fabs(vv + 14.1) < tol) alpha = 0.0015;
		else alpha = 0.0003*(vv + 14.1) / (1.0 - exp((vv + 14.1) / -5.0));
		if (fabs(vv - 3.3328) < tol) beta = 3.7836118e-4;
		else beta = 0.000073898*(vv - 3.3328) / (exp((vv - 3.3328) / 5.1237) - 1.0);
		expdt_tau_xr[i] = exp(-dt / (1 / (alpha + beta)));
		xr_infinity[i] = 1 / (1.0 + exp((vv + 14.1) / -6.5));

		// w
		if (fabs(vv - 7.9) < tol) expdt_tau_w[i] = exp(-dt / (6.0*0.2 / 1.3));
		else expdt_tau_w[i] = exp(-dt / (6.0*(1.0 - exp(-(vv - 7.9) / 5.0)) / ((1.0 + 0.3*exp(-(vv - 7.9) / 5.0))*1.0*(vv - 7.9))));
		w_infinity[i] = 1.0 - 1.0 / (1.0 + exp(-(vv - 40.0) / 17.0));

		// xs
		if (fabs(vv - 19.9) < tol) alpha = 0.00068;
		else alpha = 0.00004*(vv - 19.9) / (1.0 - exp((vv - 19.9) / -17.0));
		if (fabs(vv - 19.9) < tol) beta = 0.000315;
		else beta = 0.000035*(vv - 19.9) / (exp((vv - 19.9) / 9.0) - 1.0);
		expdt_tau_xs[i] = exp(-dt / (0.5 * 1 / (alpha + beta)));
		xs_infinity[i] = pow(1.0 + exp((vv - 19.9) / -12.7), -0.5);

		// d
		d_infinity[i] = 1.0 / (1.0 + exp((vv + 10.0) / -8.0));
		if (fabs(vv + 10.0) < tol) expdt_tau_d[i] = exp(-dt / (4.579 / (1.0 + exp((vv + 10.0) / -6.24))));
		else expdt_tau_d[i] = exp(-dt / ((1.0 - exp((vv + 10.0) / -6.24)) / (0.035*(vv + 10.0)*(1.0 + exp((vv + 10.0) / -6.24)))));

		// oa
		alpha = 0.65*1.0 / (exp((vv - (-10.0)) / -8.5) + exp((vv - (-10.0) - 40.0) / -59.0));
		beta = 0.65*1.0 / (2.5 + exp((vv - (-10.0) + 72.0) / 17.0));
		expdt_tau_oa[i] = exp(-dt / (1 / (alpha + beta) / K_Q10));
		oa_infinity[i] = 1.0 / (1.0 + exp((vv - (-10.0) + 10.47) / -17.54));

		// oi
		alpha = 1.0 / (18.53 + 1.0*exp((vv - (-10.0) + 103.7) / 10.95));
		beta = 1.0 / (35.56 + 1.0*exp((vv - (-10.0) - 8.74) / -7.44));
		expdt_tau_oi[i] = exp(-dt / (1 / (alpha + beta) / K_Q10));
		oi_infinity[i] = 1.0 / (1.0 + exp((vv - (-10.0) + 33.1) / 5.3));

		// ua
		alpha = 0.65*1.0 / (exp((vv - (-10.0)) / -8.5) + exp((vv - (-10.0) - 40.0) / -59.0));
		beta = 0.65*1.0 / (2.5 + exp((vv - (-10.0) + 72.0) / 17.0));
		expdt_tau_ua[i] = exp(-dt / (1 / (alpha + beta) / K_Q10));
		ua_infinity[i] = 1.0 / (1.0 + exp((vv - (-10.0) + 20.3) / -9.6));

		// ui
		alpha = 1.0 / (21.0 + 1.0*exp((vv - (-10.0) - 195.0) / -28.0));
		beta = 1.0 / exp((vv - (-10.0) - 168.0) / -16.0);
		expdt_tau_ui[i] = exp(-dt / (1.0 / (alpha + beta) / K_Q10));
		ui_infinity[i] = 1.0 / (1.0 + exp((vv - (-10.0) - 109.45) / 27.48));

		// f
		f_infinity[i] = exp(-(vv + 28.0) / 6.9) / (1.0 + exp(-(vv + 28.0) / 6.9));
		expdt_tau_f[i] = exp(-dt / (9.0*1.0 / (0.0197*exp(-0.0337*0.0337*(vv + 10.0)*(vv + 10.0)) + 0.02)));

		// j_NaK
		j_NaK[i] = Cm*i_NaK_max / (1.0 + 0.1245*exp(-0.1*F*vv / (R*T)) + 0.0365*sigma*exp(-F*vv / (R*T)))*K_o / (K_o + Km_K_o);

		// j_NaCa
		j_NaCa1[i] = Cm*I_NaCa_max*exp(gamma*F*vv / (R*T))*Ca_o / ((K_mNa*K_mNa*K_mNa + Na_o*Na_o*Na_o)*(K_mCa + Ca_o)*(1.0 + K_sat*exp((gamma - 1.0)*vv*F / (R*T))));
		j_NaCa2[i] = Cm*I_NaCa_max*exp((gamma - 1.0)*F*vv / (R*T))*Na_o*Na_o*Na_o / ((K_mNa*K_mNa*K_mNa + Na_o*Na_o*Na_o)*(K_mCa + Ca_o)*(1.0 + K_sat*exp((gamma - 1.0)*vv*F / (R*T))));

		// g_Kur
		g_Kur[i] = Cm*GKur*(0.005 + 0.05 / (1.0 + exp((vv - 15.0) / -13.0)));

		// j_K1
		j_K1[i] = Cm*g_K1 / (1.0 + exp(0.07*(vv + 80.0)));

		// j_Kr
		j_Kr[i] = Cm*g_Kr / (1.0 + exp((vv + 15.0) / 22.4));
	}
}


__global__ void CRNKernel(double* u, double* vc, double* w, double* d, double*  f_Ca, double* f, double* h, double* jj, double* m, double* Ca_i, double* Ca_rel, double* Ca_up, double* K_i, double* Na_i, double* v, double* xr, double* xs, double* oa, double* oi, double* ua, double* ui,
	double* volt_, bool *stim, const int num,
	const double *expdt_tau_m, const double *m_inf, const double *expdt_tau_h, const double *h_inf, const double *expdt_tau_j, const double *j_inf, const double *expdt_tau_xr, const double *xr_infinity,
	const double *expdt_tau_w, const double *w_infinity, const double *expdt_tau_xs, const double *xs_infinity, const double *expdt_tau_d, const double *d_infinity, const double *expdt_tau_oa, const double *oa_infinity,
	const double *expdt_tau_oi, const double *oi_infinity, const double *expdt_tau_ua, const double *ua_infinity, const double *expdt_tau_ui, const double *ui_infinity, const double *expdt_tau_f, const double *f_infinity, const double *j_NaK,
	const double *j_NaCa1, const double *j_NaCa2, const double *g_Kur, const double *j_K1, const double *j_Kr, const double *aa)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int bdim = blockDim.x;
	const unsigned int gdim = gridDim.x;

	int step = bdim*gdim;

	double qv1, qv2, temp1, temp2;
	int vindex1, vindex2;

	double Vm, Nai, Cai, Ki;
	double i_up_leak, i_rel, i_Ca_L, i_NaCa, Fn, i_up, E_Ca, E_Na, i_B_Na, i_B_Ca, E_K, i_B_K, i_Na;
	double i_NaK, i_K1, i_to, i_Kur, i_Kr, i_Ks, i_CaP, i_tr;

	for (int id = bid*bdim + tid; id<num; id += step)
	{
		Vm = volt_[id];
		Nai = Na_i[id];
		Cai = Ca_i[id];
		Ki = K_i[id];

		v[id] = Vm;

		// Lookup table index
		temp1 = (Vm - vlo) / dvt;
		vindex1 = (int)temp1;
		vindex2 = vindex1 + 1;
		qv1 = temp1 - (double)vindex1;
		qv2 = (double)vindex2 - temp1;

		i_up_leak = I_up_max*Ca_up[id] / Ca_up_max;

		i_rel = K_rel*u[id] * u[id] * vc[id] * w[id] * (Ca_rel[id] - Cai);

		i_Ca_L = Cm*g_Ca_L*(1.0 + aa[id]*0.01)*d[id] * f[id] * f_Ca[id] * (Vm - 65.0);

		temp1 = j_NaCa1[vindex1] * qv2 + j_NaCa1[vindex2] * qv1;
		temp2 = j_NaCa2[vindex1] * qv2 + j_NaCa2[vindex2] * qv1;
		i_NaCa = temp1*Nai*Nai*Nai - temp2*Cai;

		Fn = 1.0e-12*V_rel*i_rel - 5.0e-13 / F*(0.5*i_Ca_L - 0.2*i_NaCa);

		temp1 = 1.0 / (1.0 + exp(-(Fn - 3.4175e-13) / 13.67e-16));
		u[id] = temp1 - (temp1 - u[id])*expdt_tau_u;

		temp1 = 1.0 - 1.0 / (1.0 + exp(-(Fn - 6.835e-14) / 13.67e-16));
		temp2 = 1.91 + 2.09*1.0 / (1.0 + exp(-(Fn - 3.4175e-13) / 13.67e-16));
		vc[id] = temp1 - (temp1 - vc[id])*exp(-dt / temp2);

		temp1 = w_infinity[vindex1] * qv2 + w_infinity[vindex2] * qv1;
		temp2 = expdt_tau_w[vindex1] * qv2 + expdt_tau_w[vindex2] * qv1;
		w[id] = temp1 - (temp1 - w[id])*temp2;

		i_up = I_up_max / (1.0 + K_up / Cai);


		temp1 = d_infinity[vindex1] * qv2 + d_infinity[vindex2] * qv1;
		temp2 = expdt_tau_d[vindex1] * qv2 + expdt_tau_d[vindex2] * qv1;
		d[id] = temp1 - (temp1 - d[id])*temp2;

		temp1 = 1.0 / (1.0 + Cai / 0.00035);
		f_Ca[id] = temp1 - (temp1 - f_Ca[id])*expdt_tau_f_Ca;

		temp1 = f_infinity[vindex1] * qv2 + f_infinity[vindex2] * qv1;
		temp2 = expdt_tau_f[vindex1] * qv2 + expdt_tau_f[vindex2] * qv1;
		f[id] = temp1 - (temp1 - f[id])*temp2;

		E_Ca = RT_F / 2.0*log(Ca_o / Cai);
		E_Na = RT_F*log(Na_o / Nai);
		i_B_Na = Cm*g_B_Na*(Vm - E_Na);
		i_B_Ca = Cm*g_B_Ca*(Vm - E_Ca);
		E_K = RT_F*log(K_o / Ki);
		i_B_K = Cm*g_B_K*(Vm - E_K);

		i_Na = Cm*g_Na*m[id] * m[id] * m[id] * h[id] * jj[id] * (Vm - E_Na);


		temp1 = h_inf[vindex1] * qv2 + h_inf[vindex2] * qv1;
		temp2 = expdt_tau_h[vindex1] * qv2 + expdt_tau_h[vindex2] * qv1;
		h[id] = temp1 - (temp1 - h[id])*temp2;

		temp1 = j_inf[vindex1] * qv2 + j_inf[vindex2] * qv1;
		temp2 = expdt_tau_j[vindex1] * qv2 + expdt_tau_j[vindex2] * qv1;
		jj[id] = temp1 - (temp1 - jj[id])*temp2;

		temp1 = m_inf[vindex1] * qv2 + m_inf[vindex2] * qv1;
		temp2 = expdt_tau_m[vindex1] * qv2 + expdt_tau_m[vindex2] * qv1;
		m[id] = temp1 - (temp1 - m[id])*temp2;

		temp1 = j_NaK[vindex1] * qv2 + j_NaK[vindex2] * qv1;
		temp2 = Km_Na_i / Nai;
		i_NaK = temp1 / (1.0 + temp2*sqrt(temp2));

		temp1 = j_K1[vindex1] * qv2 + j_K1[vindex2] * qv1;
		i_K1 = (1.0 + aa[id]*0.05)*temp1*(Vm - E_K);

		i_to = Cm*g_to*(1.0 - aa[id]*0.32)*oa[id] * oa[id] * oa[id] * oi[id] * (Vm - E_K);

		temp1 = g_Kur[vindex1] * qv2 + g_Kur[vindex2] * qv1;
		i_Kur = temp1*ua[id] * ua[id] * ua[id] * ui[id] * (Vm - E_K);

		temp1 = j_Kr[vindex1] * qv2 + j_Kr[vindex2] * qv1;
		i_Kr = (1.0 - aa[id]*0.30) * temp1*xr[id] * (Vm - E_K);

		i_Ks = Cm*g_Ks*xs[id] * xs[id] * (Vm - E_K);
		i_CaP = Cm*i_CaP_max*Cai / (0.0005 + Cai);
		i_tr = (Ca_up[id] - Ca_rel[id]) / tau_tr;


		temp1 = xr_infinity[vindex1] * qv2 + xr_infinity[vindex2] * qv1;
		temp2 = expdt_tau_xr[vindex1] * qv2 + expdt_tau_xr[vindex2] * qv1;
		xr[id] = temp1 - (temp1 - xr[id])*temp2;

		temp1 = xs_infinity[vindex1] * qv2 + xs_infinity[vindex2] * qv1;
		temp2 = expdt_tau_xs[vindex1] * qv2 + expdt_tau_xs[vindex2] * qv1;
		xs[id] = temp1 - (temp1 - xs[id])*temp2;

		temp1 = oa_infinity[vindex1] * qv2 + oa_infinity[vindex2] * qv1;
		temp2 = expdt_tau_oa[vindex1] * qv2 + expdt_tau_oa[vindex2] * qv1;
		oa[id] = temp1 - (temp1 - oa[id])*temp2;

		temp1 = oi_infinity[vindex1] * qv2 + oi_infinity[vindex2] * qv1;
		temp2 = expdt_tau_oi[vindex1] * qv2 + expdt_tau_oi[vindex2] * qv1;
		oi[id] = temp1 - (temp1 - oi[id])*temp2;

		temp1 = ua_infinity[vindex1] * qv2 + ua_infinity[vindex2] * qv1;
		temp2 = expdt_tau_ua[vindex1] * qv2 + expdt_tau_ua[vindex2] * qv1;
		ua[id] = temp1 - (temp1 - ua[id])*temp2;

		temp1 = ui_infinity[vindex1] * qv2 + ui_infinity[vindex2] * qv1;
		temp2 = expdt_tau_ui[vindex1] * qv2 + expdt_tau_ui[vindex2] * qv1;
		ui[id] = temp1 - (temp1 - ui[id])*temp2;


		temp1 = (2.0*i_NaCa - (i_CaP + i_Ca_L + i_B_Ca)) / (2.0*V_i*F) + (V_up*(i_up_leak - i_up) + i_rel*V_rel) / V_i; // B1
		temp2 = 1.0 + TRPN_max*Km_TRPN / ((Cai + Km_TRPN)*(Cai + Km_TRPN)) + CMDN_max*Km_CMDN / ((Cai + Km_CMDN)*(Cai + Km_CMDN)); // B2


		if (stim[id] && Vm<0.0)
		{
			Vm -= (stim_amplitude + i_Na + i_K1 + i_to + i_Kur + i_Kr + i_Ks + i_B_Na + i_B_Ca + i_NaK + i_CaP + i_NaCa + i_Ca_L) / Cm*dt;

			Cai += (temp1 / temp2)*dt;

			Ca_up[id] += (i_up - (i_up_leak + i_tr*V_rel / V_up))*dt;

			Ca_rel[id] += ((i_tr - i_rel)*1.0 / (1.0 + CSQN_max*Km_CSQN / ((Ca_rel[id] + Km_CSQN)*(Ca_rel[id] + Km_CSQN))))*dt;

			Nai += (-3.0*i_NaK - (3.0*i_NaCa + i_B_Na + i_Na)) / (V_i*F)*dt;

			Ki += (2.0*i_NaK - (i_K1 + i_to + i_Kur + i_Kr + i_Ks + i_B_K)) / (V_i*F)*dt;

		}
		else
		{
			Vm -= (i_Na + i_K1 + i_to + i_Kur + i_Kr + i_Ks + i_B_Na + i_B_Ca + i_NaK + i_CaP + i_NaCa + i_Ca_L) / Cm*dt;

			Cai += (temp1 / temp2)*dt;

			Ca_up[id] += (i_up - (i_up_leak + i_tr*V_rel / V_up))*dt;

			Ca_rel[id] += ((i_tr - i_rel)*1.0 / (1.0 + CSQN_max*Km_CSQN / ((Ca_rel[id] + Km_CSQN)*(Ca_rel[id] + Km_CSQN))))*dt;

			Nai += (-3.0*i_NaK - (3.0*i_NaCa + i_B_Na + i_Na)) / (V_i*F)*dt;

			Ki += (2.0*i_NaK - (i_K1 + i_to + i_Kur + i_Kr + i_Ks + i_B_K)) / (V_i*F)*dt;
		}

		volt_[id] = Vm;
		Na_i[id] = Nai;
		Ca_i[id] = Cai;
		K_i[id] = Ki;
	}
}


__global__ void calcDfudtdx2(double *Dfudtdx2_x, double *Dfudtdx2_y, double *aa_, const int num)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int bdim = blockDim.x;
	const unsigned int gdim = gridDim.x;
	int step = bdim*gdim;

	double g_xx, g_yy, aa;
	int x, y;

	for (int id = bid * bdim + tid; id<num; id += step)
	{
		x = (int)(id / M);
		y = id%M;

		aa = GEOMETRY_A * exp(ALPHA * x / (M - 1)) / exp(ALPHA);
		g_xx = 1.0 + pow(aa * ALPHA / (M - 1) * sin(TWOPI / GEOMETRY_B*x)*cos(TWOPI / GEOMETRY_B*y) + aa * TWOPI / GEOMETRY_B*cos(TWOPI / GEOMETRY_B*x)*cos(TWOPI / GEOMETRY_B*y), 2);
		g_yy = 1.0 + pow(aa * TWOPI / GEOMETRY_B*sin(TWOPI / GEOMETRY_B*x)*sin(TWOPI / GEOMETRY_B*y), 2);

		aa_[id] = exp(ALPHA * x / (M - 1)) / exp(ALPHA);
		Dfudtdx2_x[id] = DiffCoeff*dt / (dx_*dx_) / 2.0 / g_xx;
		Dfudtdx2_y[id] = DiffCoeff*dt / (dx_*dx_) / 2.0 / g_yy;
	}
}


__global__ void diffKernel(double *volt_after, double *volt_before, const bool *nonExcitable, const double *Dfudtdx2_x, const double *Dfudtdx2_y, const int num)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int bdim = blockDim.x;
	const unsigned int gdim = gridDim.x;
	int step = bdim*gdim;

	double volt;
	int x, y;

	for (int id = bid * bdim + tid; id<num; id += step)
	{
		if (nonExcitable[id]) continue;

		volt = volt_before[id];
		volt_after[id] = volt;

		x = (int)(id / M);
		y = id%M;

		if (x - 1 >= 0)
		{
			volt_after[id] += (Dfudtdx2_x[id] + Dfudtdx2_x[id - M]) / 2.0*(volt_before[id - M] - volt);
		}
		if (x + 1 <M)
		{
			volt_after[id] += (Dfudtdx2_x[id] + Dfudtdx2_x[id + M]) / 2.0*(volt_before[id + M] - volt);
		}
		if (y - 1 >= 0)
		{
			volt_after[id] += (Dfudtdx2_y[id] + Dfudtdx2_y[id - 1]) / 2.0*(volt_before[id - 1] - volt);
		}
		if (y + 1 <M)
		{
			volt_after[id] += (Dfudtdx2_y[id] + Dfudtdx2_y[id + 1]) / 2.0*(volt_before[id + 1] - volt);
		}
	}
}


unsigned int __stdcall printVoltage(void* arg)
{
	FILE *fout;
	char stringarr[30];

	sprintf(stringarr, "vm%d.txt", (int)(ttt_print + 0.5));
	fout = fopen(stringarr, "wb");

	fwrite(volt_print, sizeof(double), nde, fout);

	fclose(fout);

	return 0;
}



int main()
{
	int i, j, iter, pacingCnt = 0;
	bool pacingFinish = false, stimEnd = true, S2Finish = false, S2Check = false, S2Stim = false;
	double S1S2targetvolt, S2_stim_start_ttt;

	///////////////// CUDA Threads, Blocks Setting Start ///////////////
	dim3 grid(256, 1, 1);
	dim3 threads(128, 1, 1);
	///////////////// CUDA Threads, Blocks Setting End /////////////////

	printf("expdt_tau_u: %.10lf\n", exp(-dt / tau_u));
	printf("expdt_tau_f_Ca: %.10lf\n\n", exp(-dt / tau_f_Ca));

	preprocess();

	unsigned int printVoltageId;
	HANDLE printVoltageHandle;


	// ============================ Device Variables ============================
	double *volt_d, *volt_d_;
	bool *stimulation, *ifscar_d;
	double *yglob0, *yglob1, *yglob2, *yglob3, *yglob4, *yglob5, *yglob6, *yglob7, *yglob8, *yglob9, *yglob10, *yglob11, *yglob12, *yglob13, *yglob14, *yglob15, *yglob16, *yglob17, *yglob18, *yglob19, *yglob20;
	double *Dfudtdx2_x, *Dfudtdx2_y, *aa;

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
	cudaMalloc((void**)&ifscar_d, sizeof(bool)*nde);

	cudaMalloc((void**)&Dfudtdx2_x, sizeof(double)*nde);
	cudaMalloc((void**)&Dfudtdx2_y, sizeof(double)*nde);
	cudaMalloc((void**)&aa, sizeof(double)*nde);


	// Initialize cells
	CRN_initialValue();


	// Calculate metric coefficient
	calcDfudtdx2 << <grid, threads >> >(Dfudtdx2_x, Dfudtdx2_y, aa, nde);


	// CUDA memcpy
	cudaMemcpy(volt_d, volt, sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(volt_d_, volt, sizeof(double)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(stimulation, Istim, sizeof(bool)*nde, cudaMemcpyHostToDevice);
	cudaMemcpy(ifscar_d, ifscar, sizeof(bool)*nde, cudaMemcpyHostToDevice);

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


	// ------------------------ Making CRN Lookup Table  Start ------------------------
	double *expdt_tau_m, *m_inf, *expdt_tau_h, *h_inf, *expdt_tau_j, *j_inf, *expdt_tau_xr, *xr_infinity, *expdt_tau_w, *w_infinity, *expdt_tau_xs, *xs_infinity, *expdt_tau_d, *d_infinity, *expdt_tau_oa, *oa_infinity;
	double *expdt_tau_oi, *oi_infinity, *expdt_tau_ua, *ua_infinity, *expdt_tau_ui, *ui_infinity, *expdt_tau_f, *f_infinity, *j_NaK, *j_NaCa1, *j_NaCa2, *g_Kur, *j_K1, *j_Kr;
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

	setCRNTableVolt << <grid, threads >> >(expdt_tau_m, m_inf, expdt_tau_h, h_inf, expdt_tau_j, j_inf, expdt_tau_xr, xr_infinity, expdt_tau_w, w_infinity, expdt_tau_xs, xs_infinity, expdt_tau_d, d_infinity, expdt_tau_oa, oa_infinity,
		expdt_tau_oi, oi_infinity, expdt_tau_ua, ua_infinity, expdt_tau_ui, ui_infinity, expdt_tau_f, f_infinity, j_NaK, j_NaCa1, j_NaCa2, g_Kur, j_K1, j_Kr);
	cudaThreadSynchronize();
	// ------------------------ Making CRN Lookup Table  End --------------------------


	for (iter = 1; iter <= itmax; iter++)
	{
		ttt = dt*iter;

		if (!pacingFinish)
		{
			double tttt = ttt - pacingStartTime[pacingCnt];

			if (stimEnd && tttt >= 0.0 && tttt <= 4.0)
			{
				for (i = 0; i<nstim; i++) Istim[init_stim[i]] = true;
				stimEnd = false;
				cudaMemcpy(stimulation, Istim, sizeof(bool)*nde, cudaMemcpyHostToDevice);
			}

			if (!stimEnd && tttt > 4.0)
			{
				for (i = 0; i<nde; i++) Istim[i] = false;
				stimEnd = true;
				pacingCnt++;
				cudaMemcpy(stimulation, Istim, sizeof(bool)*nde, cudaMemcpyHostToDevice);
			}

			if (pacingCnt == PACING_NUM)
			{
				pacingFinish = true;
			}
		}


#ifdef S1S2Protocol
		if (!S2Finish)
		{
			cudaMemcpy(&S1S2targetvolt, &volt_d[index(M / 2, M / 2)], sizeof(double), cudaMemcpyDeviceToHost);

			if (pacingFinish && !S2Check && S1S2targetvolt > 0.0) S2Check = true;

			if (S2Check && S1S2targetvolt < S1S2_TARGET_VOLT)
			{
				S2Finish = true;

				cudaMemcpy(volt, volt_d, sizeof(double)*nde, cudaMemcpyDeviceToHost);

				for (i = 0; i<M; i++)
				{
					for (j = 0; j<M / 2; j++)
					{
						volt[index(i, j)] = 0.0;
					}
				}

				cudaMemcpy(volt_d, volt, sizeof(double)*nde, cudaMemcpyHostToDevice);
			}
		}
#endif


		CRNKernel << <grid, threads >> >(yglob0, yglob1, yglob2, yglob3, yglob4, yglob5, yglob6, yglob7, yglob8, yglob9, yglob10, yglob11, yglob12, yglob13, yglob14, yglob15, yglob16, yglob17, yglob18, yglob19, yglob20, volt_d, stimulation, nde,
			expdt_tau_m, m_inf, expdt_tau_h, h_inf, expdt_tau_j, j_inf, expdt_tau_xr, xr_infinity, expdt_tau_w, w_infinity, expdt_tau_xs, xs_infinity, expdt_tau_d, d_infinity, expdt_tau_oa, oa_infinity,
			expdt_tau_oi, oi_infinity, expdt_tau_ua, ua_infinity, expdt_tau_ui, ui_infinity, expdt_tau_f, f_infinity, j_NaK, j_NaCa1, j_NaCa2, g_Kur, j_K1, j_Kr, aa);

		diffKernel << <grid, threads >> >(volt_d_, volt_d, ifscar_d, Dfudtdx2_x, Dfudtdx2_y, nde);

		diffKernel << <grid, threads >> >(volt_d, volt_d_, ifscar_d, Dfudtdx2_x, Dfudtdx2_y, nde);


		if (iter % 500 == 0)
		{
			printf("iter: %6d   time(ms): %6.1lf\n", iter, ttt);

			if (iter > 500)
			{
				WaitForSingleObject(printVoltageHandle, INFINITE);
				Sleep(5);
				CloseHandle(printVoltageHandle);
			}
			cudaMemcpy(volt_print, volt_d, sizeof(double)*nde, cudaMemcpyDeviceToHost);
			ttt_print = ttt;
			printVoltageHandle = (HANDLE)_beginthreadex(NULL, 0, printVoltage, NULL, 0, (unsigned*)&printVoltageId);
		}
	}


	return 0;
}
