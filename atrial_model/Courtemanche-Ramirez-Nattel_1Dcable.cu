// ------------------------------ PROGRAM SETTING START ------------------------------
#define g_Na  7.8   // nanoS_per_picoF (in fast_sodium_current)
#define g_K1  0.09*1.5   // nanoS_per_picoF (in time_independent_potassium_current)
#define g_to  0.1652*0.2   // nanoS_per_picoF (in transient_outward_K_current)
#define g_Kr  0.029411765   // nanoS_per_picoF (in rapid_delayed_rectifier_K_current)
#define g_Ca_L  0.12375*0.6   // nanoS_per_picoF (in L_type_Ca_channel)
#define GKur 0.5 // dimensionless

#define N_ 500
#define Dfu 0.00154*0.3 // cm^2/ms
#define dx_ 0.0125 // cm
// ------------------------------ PROGRAM SETTING END --------------------------------


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>



#define dt 0.05 // Depend on host program(Fortran) setting
#define dtinv 20.0 // dtinv = 1.0 / dt
#define itmax 1000000 // time limit
#define tol 0.000001

double pacingStartTime[200];
int PACING_NUM;



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
#define CRNTableVolt 14400
#define vlo (-90.0)
#define dvt 0.01

#define expdt_tau_u_1 0.993769491
#define expdt_tau_u_10 0.998750781
#define expdt_f_Ca_1 0.975309912
#define expdt_f_Ca_10 0.995012479
// --------------------------------------------------

int nstim, *init_stim;
bool *Istim, *ifscar;
double *yglob[NUM_CRN_VAR], yinitial[NUM_CRN_VAR], *volt, volt_print, ttt;



void setting()
{
	int i;

	for(i=0; i<8; i++) pacingStartTime[i] = 600.0*(double)i;
	for(i=8*1; i<8*2; i++) pacingStartTime[i] = pacingStartTime[i-1] + 400.0;
	for(i=8*2; i<8*3; i++) pacingStartTime[i] = pacingStartTime[i-1] + 300.0;
	for(i=8*3; i<8*4; i++) pacingStartTime[i] = pacingStartTime[i-1] + 250.0;
	for(i=8*4; i<8*5; i++) pacingStartTime[i] = pacingStartTime[i-1] + 220.0;
	for(i=8*5; i<8*6; i++) pacingStartTime[i] = pacingStartTime[i-1] + 200.0;
	for(i=8*6; i<8*7; i++) pacingStartTime[i] = pacingStartTime[i-1] + 190.0;
	for(i=8*7; i<8*8; i++) pacingStartTime[i] = pacingStartTime[i-1] + 180.0;
	for(i=8*8; i<8*9; i++) pacingStartTime[i] = pacingStartTime[i-1] + 170.0;
	for(i=8*9; i<8*10; i++) pacingStartTime[i] = pacingStartTime[i-1] + 160.0;
	for(i=8*10; i<8*11; i++) pacingStartTime[i] = pacingStartTime[i-1] + 150.0;
	for(i=8*11; i<8*12; i++) pacingStartTime[i] = pacingStartTime[i-1] + 120.0;
	PACING_NUM = 8*12;
}


void preprocess()
{
	int i;

	for(i=0; i<NUM_CRN_VAR; i++) yglob[i] = (double*)malloc(sizeof(double)*N_);
	volt = (double*)malloc(sizeof(double)*N_);
	Istim = (bool*)malloc(sizeof(bool)*N_);
	ifscar = (bool*)malloc(sizeof(bool)*N_);

	for(i=0; i<N_; i++)
	{
		Istim[i] = false;
		ifscar[i] = false;
	}

	nstim = 10;
	init_stim = (int*)malloc(sizeof(int)*nstim);
	for(i=0; i<nstim; i++) init_stim[i] = i;
}


void CRNcompute(double *yglob_0, double *yglob_1, double *yglob_2, double *yglob_3, double *yglob_4, double *yglob_5, double *yglob_6, double *yglob_7, double *yglob_8, double *yglob_9, double *yglob_10, double *yglob_11, double *yglob_12, double *yglob_13, double *yglob_14, double *yglob_15, double *yglob_16, double *yglob_17, double *yglob_18, double *yglob_19, double *yglob_20,
	double *volt_, const bool *stimulation, int cellNum)
{
	//double Ca_CMDN;   // millimolar (in Ca_buffers)
	//double Ca_CSQN;   // millimolar (in Ca_buffers)
	//double Ca_TRPN;   // millimolar (in Ca_buffers)
	double i_up_leak;   // millimolar_per_millisecond (in Ca_leak_current_by_the_NSR)
	double u_infinity;   // dimensionless (in Ca_release_current_from_JSR_u_gate)
	double tau_v;   // millisecond (in Ca_release_current_from_JSR_v_gate)
	double v_infinity;   // dimensionless (in Ca_release_current_from_JSR_v_gate)
	double tau_w;   // millisecond (in Ca_release_current_from_JSR_w_gate)
	double w_infinity;   // dimensionless (in Ca_release_current_from_JSR_w_gate)
	double Fn;   // dimensionless (in Ca_release_current_from_JSR)
	double i_rel;   // millimolar_per_millisecond (in Ca_release_current_from_JSR)
	double i_up;   // millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)
	double d_infinity;   // dimensionless (in L_type_Ca_channel_d_gate)
	double tau_d;   // millisecond (in L_type_Ca_channel_d_gate)
	double f_Ca_infinity;   // dimensionless (in L_type_Ca_channel_f_Ca_gate)
	double f_infinity;   // dimensionless (in L_type_Ca_channel_f_gate)
	double tau_f;   // millisecond (in L_type_Ca_channel_f_gate)
	double i_Ca_L;   // picoA (in L_type_Ca_channel)
	double i_NaCa;   // picoA (in Na_Ca_exchanger_current)
	double E_Ca;   // millivolt (in background_currents)
	double i_B_Ca;   // picoA (in background_currents)
	double i_B_K;   // picoA (in background_currents)
	double i_B_Na;   // picoA (in background_currents)
	double alpha_h;   // per_millisecond (in fast_sodium_current_h_gate)
	double beta_h;   // per_millisecond (in fast_sodium_current_h_gate)
	double h_inf;   // dimensionless (in fast_sodium_current_h_gate)
	double tau_h;   // millisecond (in fast_sodium_current_h_gate)
	double alpha_j;   // per_millisecond (in fast_sodium_current_j_gate)
	double beta_j;   // per_millisecond (in fast_sodium_current_j_gate)
	double j_inf;   // dimensionless (in fast_sodium_current_j_gate)
	double tau_j;   // millisecond (in fast_sodium_current_j_gate)
	double alpha_m;   // per_millisecond (in fast_sodium_current_m_gate)
	double beta_m;   // per_millisecond (in fast_sodium_current_m_gate)
	double m_inf;   // dimensionless (in fast_sodium_current_m_gate)
	double tau_m;   // millisecond (in fast_sodium_current_m_gate)
	double E_Na;   // millivolt (in fast_sodium_current)
	double i_Na;   // picoA (in fast_sodium_current)
	double B1;   // millimolar_per_millisecond (in intracellular_ion_concentrations)
	double B2;   // dimensionless (in intracellular_ion_concentrations)
	double i_st;   // picoA (in membrane)
	double alpha_xr;   // per_millisecond (in rapid_delayed_rectifier_K_current_xr_gate)
	double beta_xr;   // per_millisecond (in rapid_delayed_rectifier_K_current_xr_gate)
	double tau_xr;   // millisecond (in rapid_delayed_rectifier_K_current_xr_gate)
	double xr_infinity;   // dimensionless (in rapid_delayed_rectifier_K_current_xr_gate)
	double i_Kr;   // picoA (in rapid_delayed_rectifier_K_current)
	double i_CaP;   // picoA (in sarcolemmal_calcium_pump_current)
	double alpha_xs;   // per_millisecond (in slow_delayed_rectifier_K_current_xs_gate)
	double beta_xs;   // per_millisecond (in slow_delayed_rectifier_K_current_xs_gate)
	double tau_xs;   // millisecond (in slow_delayed_rectifier_K_current_xs_gate)
	double xs_infinity;   // dimensionless (in slow_delayed_rectifier_K_current_xs_gate)
	double i_Ks;   // picoA (in slow_delayed_rectifier_K_current)
	double f_NaK;   // dimensionless (in sodium_potassium_pump)
	double i_NaK;   // picoA (in sodium_potassium_pump)
	double E_K;   // millivolt (in time_independent_potassium_current)
	double i_K1;   // picoA (in time_independent_potassium_current)
	double i_tr;   // millimolar_per_millisecond (in transfer_current_from_NSR_to_JSR)
	double alpha_oa;   // per_millisecond (in transient_outward_K_current_oa_gate)
	double beta_oa;   // per_millisecond (in transient_outward_K_current_oa_gate)
	double oa_infinity;   // dimensionless (in transient_outward_K_current_oa_gate)
	double tau_oa;   // millisecond (in transient_outward_K_current_oa_gate)
	double alpha_oi;   // per_millisecond (in transient_outward_K_current_oi_gate)
	double beta_oi;   // per_millisecond (in transient_outward_K_current_oi_gate)
	double oi_infinity;   // dimensionless (in transient_outward_K_current_oi_gate)
	double tau_oi;   // millisecond (in transient_outward_K_current_oi_gate)
	double i_to;   // picoA (in transient_outward_K_current)
	double alpha_ua;   // per_millisecond (in ultrarapid_delayed_rectifier_K_current_ua_gate)
	double beta_ua;   // per_millisecond (in ultrarapid_delayed_rectifier_K_current_ua_gate)
	double tau_ua;   // millisecond (in ultrarapid_delayed_rectifier_K_current_ua_gate)
	double ua_infinity;   // dimensionless (in ultrarapid_delayed_rectifier_K_current_ua_gate)
	double alpha_ui;   // per_millisecond (in ultrarapid_delayed_rectifier_K_current_ui_gate)
	double beta_ui;   // per_millisecond (in ultrarapid_delayed_rectifier_K_current_ui_gate)
	double tau_ui;   // millisecond (in ultrarapid_delayed_rectifier_K_current_ui_gate)
	double ui_infinity;   // dimensionless (in ultrarapid_delayed_rectifier_K_current_ui_gate)
	double g_Kur;   // nanoS_per_picoF (in ultrarapid_delayed_rectifier_K_current)
	double i_Kur;   // picoA (in ultrarapid_delayed_rectifier_K_current)

	yglob_14[cellNum] = volt_[cellNum];

	//Ca_CMDN = CMDN_max*yglob_9[cellNum]/(yglob_9[cellNum]+Km_CMDN);
	//Ca_TRPN = TRPN_max*yglob_9[cellNum]/(yglob_9[cellNum]+Km_TRPN);
	//Ca_CSQN = CSQN_max*yglob_10[cellNum]/(yglob_10[cellNum]+Km_CSQN);
	i_up_leak = I_up_max*yglob_11[cellNum]/Ca_up_max;
	i_rel = K_rel*yglob_0[cellNum]*yglob_0[cellNum]*yglob_1[cellNum]*yglob_2[cellNum]*(yglob_10[cellNum]-yglob_9[cellNum]);
	i_Ca_L = Cm*g_Ca_L*yglob_3[cellNum]*yglob_5[cellNum]*yglob_4[cellNum]*(yglob_14[cellNum]-65.0);
	i_NaCa = Cm*I_NaCa_max*(exp(gamma*F*yglob_14[cellNum]/(R*T))*yglob_13[cellNum]*yglob_13[cellNum]*yglob_13[cellNum]*Ca_o-exp((gamma-1.0)*F*yglob_14[cellNum]/(R*T))*pow(Na_o, 3.0)*yglob_9[cellNum])/((pow(K_mNa, 3.0)+pow(Na_o, 3.0))*(K_mCa+Ca_o)*(1.0+K_sat*exp((gamma-1.0)*yglob_14[cellNum]*F/(R*T))));
	Fn = 1.0e3*(1.0e-15*V_rel*i_rel-1.0e-15/(2.0*F)*(0.5*i_Ca_L-0.2*i_NaCa));
	u_infinity = pow(1.0+exp(-(Fn-3.4175e-13)/13.67e-16), -1.0);
	yglob_0[cellNum] = u_infinity + (yglob_0[cellNum]-u_infinity)*exp(-dt/tau_u);
	tau_v = 1.91+2.09*pow(1.0+exp(-(Fn-3.4175e-13)/13.67e-16), -1.0);
	v_infinity = 1.0-pow(1.0+exp(-(Fn-6.835e-14)/13.67e-16), -1.0);
	yglob_1[cellNum] = v_infinity + (yglob_1[cellNum]-v_infinity)*exp(-dt/tau_v);

	if (fabs(yglob_14[cellNum]-7.9) < 1.0e-10)
		tau_w = 6.0*0.2/1.3;
	else
		tau_w = 6.0*(1.0-exp(-(yglob_14[cellNum]-7.9)/5.0))/((1.0+0.3*exp(-(yglob_14[cellNum]-7.9)/5.0))*1.0*(yglob_14[cellNum]-7.9));

	w_infinity = 1.0-pow(1.0+exp(-(yglob_14[cellNum]-40.0)/17.0), -1.0);
	yglob_2[cellNum] = w_infinity + (yglob_2[cellNum]-w_infinity)*exp(-dt/tau_w);
	i_up = I_up_max/(1.0+K_up/yglob_9[cellNum]);
	d_infinity = pow(1.0+exp((yglob_14[cellNum]+10.0)/-8.0), -1.0);

	if (fabs(yglob_14[cellNum]+10.0) < 1.0e-10)
		tau_d = 4.579/(1.0+exp((yglob_14[cellNum]+10.0)/-6.24));
	else
		tau_d = (1.0-exp((yglob_14[cellNum]+10.0)/-6.24))/(0.035*(yglob_14[cellNum]+10.0)*(1.0+exp((yglob_14[cellNum]+10.0)/-6.24)));

	yglob_3[cellNum] = d_infinity + (yglob_3[cellNum]-d_infinity)*exp(-dt/tau_d);
	f_Ca_infinity = pow(1.0+yglob_9[cellNum]/0.00035, -1.0);
	yglob_4[cellNum] = f_Ca_infinity + (yglob_4[cellNum]-f_Ca_infinity)*exp(-dt/tau_f_Ca);
	f_infinity = exp(-(yglob_14[cellNum]+28.0)/6.9)/(1.0+exp(-(yglob_14[cellNum]+28.0)/6.9));
	tau_f = 9.0*pow(0.0197*exp(-pow(0.0337, 2.0)*pow(yglob_14[cellNum]+10.0, 2.0))+0.02, -1.0);
	yglob_5[cellNum] = f_infinity + (yglob_5[cellNum]-f_infinity)*exp(-dt/tau_f);
	E_Ca = R*T/(2.0*F)*log(Ca_o/yglob_9[cellNum]);
	E_Na = R*T/F*log(Na_o/yglob_13[cellNum]);
	i_B_Na = Cm*g_B_Na*(yglob_14[cellNum]-E_Na);
	i_B_Ca = Cm*g_B_Ca*(yglob_14[cellNum]-E_Ca);
	E_K = R*T/F*log(K_o/yglob_12[cellNum]);
	i_B_K = Cm*g_B_K*(yglob_14[cellNum]-E_K);
	i_Na = Cm*g_Na*yglob_8[cellNum]*yglob_8[cellNum]*yglob_8[cellNum]*yglob_6[cellNum]*yglob_7[cellNum]*(yglob_14[cellNum]-E_Na);

	if (yglob_14[cellNum] < -40.0)
		alpha_h = 0.135*exp((yglob_14[cellNum]+80.0)/-6.8);
	else
		alpha_h = 0.0;

	if (yglob_14[cellNum] < -40.0)
		beta_h = 3.56*exp(0.079*yglob_14[cellNum])+3.1e5*exp(0.35*yglob_14[cellNum]);
	else
		beta_h = 1.0/(0.13*(1.0+exp((yglob_14[cellNum]+10.66)/-11.1)));

	h_inf = alpha_h/(alpha_h+beta_h);
	tau_h = 1.0/(alpha_h+beta_h);
	yglob_6[cellNum] = h_inf + (yglob_6[cellNum]-h_inf)*exp(-dt/tau_h);

	if (yglob_14[cellNum] < -40.0)
		alpha_j = (-1.2714e5*exp(0.2444*yglob_14[cellNum])-3.474e-5*exp(-0.04391*yglob_14[cellNum]))*(yglob_14[cellNum]+37.78)/(1.0+exp(0.311*(yglob_14[cellNum]+79.23)));
	else
		alpha_j = 0.0;

	if (yglob_14[cellNum] < -40.0)
		beta_j = 0.1212*exp(-0.01052*yglob_14[cellNum])/(1.0+exp(-0.1378*(yglob_14[cellNum]+40.14)));
	else
		beta_j = 0.3*exp(-2.535e-7*yglob_14[cellNum])/(1.0+exp(-0.1*(yglob_14[cellNum]+32.0)));

	j_inf = alpha_j/(alpha_j+beta_j);
	tau_j = 1.0/(alpha_j+beta_j);
	yglob_7[cellNum] = j_inf + (yglob_7[cellNum]-j_inf)*exp(-dt/tau_j);

	if (yglob_14[cellNum] == -47.13)
		alpha_m = 3.2;
	else
		alpha_m = 0.32*(yglob_14[cellNum]+47.13)/(1.0-exp(-0.1*(yglob_14[cellNum]+47.13)));

	beta_m = 0.08*exp(-yglob_14[cellNum]/11.0);
	m_inf = alpha_m/(alpha_m+beta_m);
	tau_m = 1.0/(alpha_m+beta_m);
	yglob_8[cellNum] = m_inf + (yglob_8[cellNum]-m_inf)*exp(-dt/tau_m);
	f_NaK = pow(1.0+0.1245*exp(-0.1*F*yglob_14[cellNum]/(R*T))+0.0365*sigma*exp(-F*yglob_14[cellNum]/(R*T)), -1.0);
	i_NaK = Cm*i_NaK_max*f_NaK*1.0/(1.0+pow(Km_Na_i/yglob_13[cellNum], 1.5))*K_o/(K_o+Km_K_o);
	yglob_13[cellNum] += (-3.0*i_NaK-(3.0*i_NaCa+i_B_Na+i_Na))/(V_i*F)*dt;
	i_K1 = Cm*g_K1*(yglob_14[cellNum]-E_K)/(1.0+exp(0.07*(yglob_14[cellNum]+80.0)));
	i_to = Cm*g_to*yglob_17[cellNum]*yglob_17[cellNum]*yglob_17[cellNum]*yglob_18[cellNum]*(yglob_14[cellNum]-E_K);
	g_Kur = 0.005+0.05/(1.0+exp((yglob_14[cellNum]-15.0)/-13.0));
	i_Kur = Cm*g_Kur*yglob_19[cellNum]*yglob_19[cellNum]*yglob_19[cellNum]*yglob_20[cellNum]*(yglob_14[cellNum]-E_K);
	i_Kr = Cm*g_Kr*yglob_15[cellNum]*(yglob_14[cellNum]-E_K)/(1.0+exp((yglob_14[cellNum]+15.0)/22.4));
	i_Ks = Cm*g_Ks*yglob_16[cellNum]*yglob_16[cellNum]*(yglob_14[cellNum]-E_K);
	yglob_12[cellNum] += (2.0*i_NaK-(i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_K))/(V_i*F)*dt;
	i_CaP = Cm*i_CaP_max*yglob_9[cellNum]/(0.0005+yglob_9[cellNum]);
	B1 = (2.0*i_NaCa-(i_CaP+i_Ca_L+i_B_Ca))/(2.0*V_i*F)+(V_up*(i_up_leak-i_up)+i_rel*V_rel)/V_i;
	B2 = 1.0+TRPN_max*Km_TRPN/pow(yglob_9[cellNum]+Km_TRPN, 2.0)+CMDN_max*Km_CMDN/pow(yglob_9[cellNum]+Km_CMDN, 2.0);
	yglob_9[cellNum] += B1/B2*dt;
	i_tr = (yglob_11[cellNum]-yglob_10[cellNum])/tau_tr;
	yglob_11[cellNum] += (i_up-(i_up_leak+i_tr*V_rel/V_up))*dt;
	yglob_10[cellNum] += (i_tr-i_rel)*pow(1.0+CSQN_max*Km_CSQN/pow(yglob_10[cellNum]+Km_CSQN, 2.0), -1.0)*dt;

	if (stimulation[cellNum])
		i_st = stim_amplitude;
	else
		i_st = 0.0;

	yglob_14[cellNum] += -(i_Na+i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_Na+i_B_Ca+i_NaK+i_CaP+i_NaCa+i_Ca_L+i_st)/Cm*dt;

	if (fabs(yglob_14[cellNum]+14.1) < 1.0e-10)
		alpha_xr = 0.0015;
	else
		alpha_xr = 0.0003*(yglob_14[cellNum]+14.1)/(1.0-exp((yglob_14[cellNum]+14.1)/-5.0));

	if (fabs(yglob_14[cellNum]-3.3328) < 1.0e-10)
		beta_xr = 3.7836118e-4;
	else
		beta_xr = 0.000073898*(yglob_14[cellNum]-3.3328)/(exp((yglob_14[cellNum]-3.3328)/5.1237)-1.0);

	tau_xr = 1.0/(alpha_xr+beta_xr);
	xr_infinity = pow(1.0+exp((yglob_14[cellNum]+14.1)/-6.5), -1.0);
	yglob_15[cellNum] = xr_infinity + (yglob_15[cellNum]-xr_infinity)*exp(-dt/tau_xr);

	if (fabs(yglob_14[cellNum]-19.9) < 1.0e-10)
		alpha_xs = 0.00068;
	else
		alpha_xs = 0.00004*(yglob_14[cellNum]-19.9)/(1.0-exp((yglob_14[cellNum]-19.9)/-17.0));

	if (fabs(yglob_14[cellNum]-19.9) < 1.0e-10)
		beta_xs = 0.000315;
	else
		beta_xs = 0.000035*(yglob_14[cellNum]-19.9)/(exp((yglob_14[cellNum]-19.9)/9.0)-1.0);

	tau_xs = 0.5*1.0/(alpha_xs+beta_xs);
	xs_infinity = pow(1.0+exp((yglob_14[cellNum]-19.9)/-12.7), -0.5);
	yglob_16[cellNum] = xs_infinity + (yglob_16[cellNum]-xs_infinity)*exp(-dt/tau_xs);
	alpha_oa = 0.65*pow(exp((yglob_14[cellNum]-(-10.0))/-8.5)+exp((yglob_14[cellNum]-(-10.0)-40.0)/-59.0), -1.0);
	beta_oa = 0.65*pow(2.5+exp((yglob_14[cellNum]-(-10.0)+72.0)/17.0), -1.0);
	tau_oa = 1.0/(alpha_oa+beta_oa)/K_Q10;
	oa_infinity = pow(1.0+exp((yglob_14[cellNum]-(-10.0)+10.47)/-17.54), -1.0);
	yglob_17[cellNum] = oa_infinity + (yglob_17[cellNum]-oa_infinity)*exp(-dt/tau_oa);
	alpha_oi = pow(18.53+1.0*exp((yglob_14[cellNum]-(-10.0)+103.7)/10.95), -1.0);
	beta_oi = pow(35.56+1.0*exp((yglob_14[cellNum]-(-10.0)-8.74)/-7.44), -1.0);
	tau_oi = 1.0/(alpha_oi+beta_oi)/K_Q10;
	oi_infinity = pow(1.0+exp((yglob_14[cellNum]-(-10.0)+33.1)/5.3), -1.0);
	yglob_18[cellNum] = oi_infinity + (yglob_18[cellNum]-oi_infinity)*exp(-dt/tau_oi);
	alpha_ua = 0.65*pow(exp((yglob_14[cellNum]-(-10.0))/-8.5)+exp((yglob_14[cellNum]-(-10.0)-40.0)/-59.0), -1.0);
	beta_ua = 0.65*pow(2.5+exp((yglob_14[cellNum]-(-10.0)+72.0)/17.0), -1.0);
	tau_ua = 1.0/(alpha_ua+beta_ua)/K_Q10;
	ua_infinity = pow(1.0+exp((yglob_14[cellNum]-(-10.0)+20.3)/-9.6), -1.0);
	yglob_19[cellNum] = ua_infinity + (yglob_19[cellNum]-ua_infinity)*exp(-dt/tau_ua);
	alpha_ui = pow(21.0+1.0*exp((yglob_14[cellNum]-(-10.0)-195.0)/-28.0), -1.0);
	beta_ui = 1.0/exp((yglob_14[cellNum]-(-10.0)-168.0)/-16.0);
	tau_ui = 1.0/(alpha_ui+beta_ui)/K_Q10;
	ui_infinity = pow(1.0+exp((yglob_14[cellNum]-(-10.0)-109.45)/27.48), -1.0);
	yglob_20[cellNum] = ui_infinity + (yglob_20[cellNum]-ui_infinity)*exp(-dt/tau_ui);

	volt_[cellNum] = yglob_14[cellNum];
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

	for(j=0; j<NUM_CRN_VAR; j++) yglob[j][0] = yinitial[j];
	volt[0] = yinitial[14];

	// Make resting state
	for(i=0; i<500000; i++) CRNcompute(yglob[0], yglob[1], yglob[2], yglob[3], yglob[4], yglob[5], yglob[6], yglob[7], yglob[8], yglob[9], yglob[10], yglob[11], yglob[12], yglob[13], yglob[14], yglob[15], yglob[16], yglob[17], yglob[18], yglob[19], yglob[20], volt, Istim, 0);

	for(j=0; j<NUM_CRN_VAR; j++) yinitial[j] = yglob[j][0];

	for(i=0; i<N_; i++)
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
	double* volt_ , bool *stim, const bool *nonExcitable, const double restPotential, const int num,
	const double *tau_m, const double *m_inf, const double *tau_h, const double *h_inf, const double *tau_j, const double *j_inf, const double *tau_xr, const double *xr_infinity,
	const double *tau_w, const double *w_infinity, const double *tau_xs, const double *xs_infinity, const double *tau_d, const double *d_infinity, const double *tau_oa, const double *oa_infinity,
	const double *tau_oi, const double *oi_infinity, const double *tau_ua, const double *ua_infinity, const double *tau_ui, const double *ui_infinity, const double *tau_f, const double *f_infinity, const double *j_NaK,
	const double *j_NaCa1, const double *j_NaCa2, const double *g_Kur, const double *j_K1, const double *j_Kr)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int bdim = blockDim.x;
	const unsigned int gdim = gridDim.x;

	int step = bdim*gdim;
	int i, repeat;
	double dtime;

	double qv1, qv2, temp1, temp2;
	int vindex1, vindex2;

	double Vm, Nai, Cai, Ki;
	double i_up_leak, i_rel, i_Ca_L, i_NaCa, Fn, i_up, E_Ca, E_Na, i_B_Na, i_B_Ca, E_K, i_B_K, i_Na;
	double i_NaK, i_K1, i_to, i_Kur, i_Kr, i_Ks, i_CaP, i_tr;

	for (int id = bid*bdim + tid; id<num; id+=step)
	{
		if(nonExcitable[id])
		{
			volt_[id] = restPotential;
			continue;
		}

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

			i_Ca_L = Cm*g_Ca_L*d[id]*f[id]*f_Ca[id]*(Vm-65.0);

			temp1 = j_NaCa1[vindex1]*qv2 + j_NaCa1[vindex2]*qv1;
			temp2 = j_NaCa2[vindex1]*qv2 + j_NaCa2[vindex2]*qv1;
			i_NaCa = temp1*Nai*Nai*Nai - temp2*Cai;

			Fn = 1.0e-12*V_rel*i_rel-5.0e-13/F*(0.5*i_Ca_L-0.2*i_NaCa);

			temp1 = 1.0/(1.0+exp(-(Fn-3.4175e-13)/13.67e-16));
			if(repeat == 1) u[id] = temp1-(temp1-u[id])*expdt_tau_u_1;
			else u[id] = temp1-(temp1-u[id])*expdt_tau_u_10;

			temp1 = 1.0-1.0/(1.0+exp(-(Fn-6.835e-14)/13.67e-16));
			temp2 = 1.91+2.09*1.0/(1.0+exp(-(Fn-3.4175e-13)/13.67e-16));
			vc[id] = temp1-(temp1-vc[id])*exp(-dtime/temp2);

			temp1 = w_infinity[vindex1]*qv2 + w_infinity[vindex2]*qv1;
			temp2 = tau_w[vindex1]*qv2 + tau_w[vindex2]*qv1;
			w[id] = temp1-(temp1-w[id])*exp(-dtime/temp2);

			i_up = I_up_max/(1.0+K_up/Cai);


			temp1 = d_infinity[vindex1]*qv2 + d_infinity[vindex2]*qv1;
			temp2 = tau_d[vindex1]*qv2 + tau_d[vindex2]*qv1;
			d[id] = temp1-(temp1-d[id])*exp(-dtime/temp2);

			temp1 = 1.0/(1.0+Cai/0.00035);
			if(repeat == 1) f_Ca[id] = temp1-(temp1-f_Ca[id])*expdt_f_Ca_1;
			else f_Ca[id] = temp1-(temp1-f_Ca[id])*expdt_f_Ca_10;

			temp1 = f_infinity[vindex1]*qv2 + f_infinity[vindex2]*qv1;
			temp2 = tau_f[vindex1]*qv2 + tau_f[vindex2]*qv1;
			f[id] = temp1-(temp1-f[id])*exp(-dtime/temp2);

			E_Ca = RT_F/2.0*log(Ca_o/Cai);
			E_Na = RT_F*log(Na_o/Nai);
			i_B_Na = Cm*g_B_Na*(Vm-E_Na);
			i_B_Ca = Cm*g_B_Ca*(Vm-E_Ca);
			E_K = RT_F*log(K_o/Ki);
			i_B_K = Cm*g_B_K*(Vm-E_K);

			i_Na = Cm*g_Na*m[id]*m[id]*m[id]*h[id]*jj[id]*(Vm-E_Na);


			temp1 = h_inf[vindex1]*qv2 + h_inf[vindex2]*qv1;
			temp2 = tau_h[vindex1]*qv2 + tau_h[vindex2]*qv1;
			h[id] = temp1-(temp1-h[id])*exp(-dtime/temp2);

			temp1 = j_inf[vindex1]*qv2 + j_inf[vindex2]*qv1;
			temp2 = tau_j[vindex1]*qv2 + tau_j[vindex2]*qv1;
			jj[id] = temp1-(temp1-jj[id])*exp(-dtime/temp2);

			temp1 = m_inf[vindex1]*qv2 + m_inf[vindex2]*qv1;
			temp2 = tau_m[vindex1]*qv2 + tau_m[vindex2]*qv1;
			m[id] = temp1-(temp1-m[id])*exp(-dtime/temp2);

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
			temp2 = tau_xr[vindex1]*qv2 + tau_xr[vindex2]*qv1;
			xr[id] = temp1-(temp1-xr[id])*exp(-dtime/temp2);

			temp1 = xs_infinity[vindex1]*qv2 + xs_infinity[vindex2]*qv1;
			temp2 = tau_xs[vindex1]*qv2 + tau_xs[vindex2]*qv1;
			xs[id] = temp1-(temp1-xs[id])*exp(-dtime/temp2);

			temp1 = oa_infinity[vindex1]*qv2 + oa_infinity[vindex2]*qv1;
			temp2 = tau_oa[vindex1]*qv2 + tau_oa[vindex2]*qv1;
			oa[id] = temp1-(temp1-oa[id])*exp(-dtime/temp2);

			temp1 = oi_infinity[vindex1]*qv2 + oi_infinity[vindex2]*qv1;
			temp2 = tau_oi[vindex1]*qv2 + tau_oi[vindex2]*qv1;
			oi[id] = temp1-(temp1-oi[id])*exp(-dtime/temp2);

			temp1 = ua_infinity[vindex1]*qv2 + ua_infinity[vindex2]*qv1;
			temp2 = tau_ua[vindex1]*qv2 + tau_ua[vindex2]*qv1;
			ua[id] = temp1-(temp1-ua[id])*exp(-dtime/temp2);

			temp1 = ui_infinity[vindex1]*qv2 + ui_infinity[vindex2]*qv1;
			temp2 = tau_ui[vindex1]*qv2 + tau_ui[vindex2]*qv1;
			ui[id] = temp1-(temp1-ui[id])*exp(-dtime/temp2);


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
}


__global__ void diffKernel(double *volt_after, double *volt_before, const bool *nonExcitable, const int num)
{
    const unsigned int tid = threadIdx.x;
    const unsigned int bid = blockIdx.x;
    const unsigned int bdim = blockDim.x;
    const unsigned int gdim = gridDim.x;
    int step=bdim*gdim;

	double Dfudtdx2 = Dfu*dt/(dx_*dx_)/2.0;

	for (int id=bid * bdim + tid;id<num;id+=step)
	{
		if(id == 0) volt_after[id] = volt_before[id] + 2.0*(volt_before[id+1]-volt_before[id])*Dfudtdx2;
		else if (id == N_-1) volt_after[id] = volt_before[id] + 2.0*(volt_before[id-1]-volt_before[id])*Dfudtdx2;
		else volt_after[id] = volt_before[id] + (volt_before[id-1]+volt_before[id+1]-2.0*volt_before[id])*Dfudtdx2;
	}
}



int main()
{
	int i, iter, pacingCnt = 0;
	bool pacingFinish = false, stimEnd = true;

	// ------------------- APD Calculator Variables Start -------------------
	bool detectStart = true, detectPeak = false, detectEnd = false;
	double preVolt = 0.0, startTime, restPotential, v90;
	// ------------------- APD Calculator Variables End ---------------------

	FILE *vmrecord = fopen("VoltageRecord.txt", "wb");
	FILE *apdrecord = fopen("APDRecord.txt", "w");

	setting();

	///////////////// CUDA Threads, Blocks Setting Start ///////////////
	dim3 grid(4, 1, 1);
	dim3 threads(128, 1, 1);
	///////////////// CUDA Threads, Blocks Setting End /////////////////

	preprocess();


	// ============================ Device Variables ============================
	double *volt_d, *volt_d_;
	bool *stimulation, *ifscar_d;
	double *yglob0, *yglob1, *yglob2, *yglob3, *yglob4, *yglob5, *yglob6, *yglob7, *yglob8, *yglob9, *yglob10, *yglob11, *yglob12, *yglob13, *yglob14, *yglob15, *yglob16, *yglob17, *yglob18, *yglob19, *yglob20;
	
	cudaMalloc((void**)&yglob0, sizeof(double)*N_);
	cudaMalloc((void**)&yglob1, sizeof(double)*N_);
	cudaMalloc((void**)&yglob2, sizeof(double)*N_);
	cudaMalloc((void**)&yglob3, sizeof(double)*N_);
	cudaMalloc((void**)&yglob4, sizeof(double)*N_);
	cudaMalloc((void**)&yglob5, sizeof(double)*N_);
	cudaMalloc((void**)&yglob6, sizeof(double)*N_);
	cudaMalloc((void**)&yglob7, sizeof(double)*N_);
	cudaMalloc((void**)&yglob8, sizeof(double)*N_);
	cudaMalloc((void**)&yglob9, sizeof(double)*N_);
	cudaMalloc((void**)&yglob10, sizeof(double)*N_);
	cudaMalloc((void**)&yglob11, sizeof(double)*N_);
	cudaMalloc((void**)&yglob12, sizeof(double)*N_);
	cudaMalloc((void**)&yglob13, sizeof(double)*N_);
	cudaMalloc((void**)&yglob14, sizeof(double)*N_);
	cudaMalloc((void**)&yglob15, sizeof(double)*N_);
	cudaMalloc((void**)&yglob16, sizeof(double)*N_);
	cudaMalloc((void**)&yglob17, sizeof(double)*N_);
	cudaMalloc((void**)&yglob18, sizeof(double)*N_);
	cudaMalloc((void**)&yglob19, sizeof(double)*N_);
	cudaMalloc((void**)&yglob20, sizeof(double)*N_);

	cudaMalloc((void**)&volt_d, sizeof(double)*N_);
	cudaMalloc((void**)&volt_d_, sizeof(double)*N_);
	cudaMalloc((void**)&stimulation, sizeof(bool)*N_);
	cudaMalloc((void**)&ifscar_d, sizeof(bool)*N_);


	// Initialize cells
	CRN_initialValue();


	// CUDA memcpy
	cudaMemcpy(volt_d, volt, sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(volt_d_, volt, sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(stimulation, Istim, sizeof(bool)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(ifscar_d, ifscar, sizeof(bool)*N_, cudaMemcpyHostToDevice);

	cudaMemcpy(yglob0, yglob[0], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob1, yglob[1], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob2, yglob[2], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob3, yglob[3], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob4, yglob[4], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob5, yglob[5], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob6, yglob[6], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob7, yglob[7], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob8, yglob[8], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob9, yglob[9], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob10, yglob[10], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob11, yglob[11], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob12, yglob[12], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob13, yglob[13], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob14, yglob[14], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob15, yglob[15], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob16, yglob[16], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob17, yglob[17], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob18, yglob[18], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob19, yglob[19], sizeof(double)*N_, cudaMemcpyHostToDevice);
	cudaMemcpy(yglob20, yglob[20], sizeof(double)*N_, cudaMemcpyHostToDevice);


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

	
	for(iter=1; iter<itmax; iter++)
	{
		ttt = dt*iter;

		if(!pacingFinish)
		{
			double tttt = ttt - pacingStartTime[pacingCnt];

			if(stimEnd && tttt >= 0.0 && tttt <= 4.0)
			{
				for(i=0; i<nstim; i++) Istim[init_stim[i]] = true;
				stimEnd = false;
				cudaMemcpy(stimulation, Istim, sizeof(bool)*N_, cudaMemcpyHostToDevice);
			}

			if(!stimEnd && tttt > 4.0)
			{
				for(i=0; i<N_; i++) Istim[i] = false;
				stimEnd = true;
				pacingCnt++;
				cudaMemcpy(stimulation, Istim, sizeof(bool)*N_, cudaMemcpyHostToDevice);
			}

			if(pacingCnt == PACING_NUM)
			{
				pacingFinish = true;
			}
		}


		CRNKernel<<<grid, threads>>>(yglob0, yglob1, yglob2, yglob3, yglob4, yglob5, yglob6, yglob7, yglob8, yglob9, yglob10, yglob11, yglob12, yglob13, yglob14, yglob15, yglob16, yglob17, yglob18, yglob19, yglob20, volt_d, stimulation, ifscar_d, yinitial[14], N_,
			tau_m, m_inf, tau_h, h_inf, tau_j, j_inf, tau_xr, xr_infinity, tau_w, w_infinity, tau_xs, xs_infinity, tau_d, d_infinity, tau_oa, oa_infinity,
			tau_oi, oi_infinity, tau_ua, ua_infinity, tau_ui, ui_infinity, tau_f, f_infinity, j_NaK, j_NaCa1, j_NaCa2, g_Kur, j_K1, j_Kr);
		cudaThreadSynchronize();

		diffKernel<<<grid, threads>>>(volt_d_, volt_d, ifscar_d, N_);
		cudaThreadSynchronize();

		diffKernel<<<grid, threads>>>(volt_d, volt_d_, ifscar_d, N_);
		cudaThreadSynchronize();


		if(iter%10 == 0)
		{
			cudaMemcpy(&volt_print, &volt_d[N_/2], sizeof(double), cudaMemcpyDeviceToHost);
			fprintf(vmrecord, "%.1lf %.2lf\n", ttt, volt_print);

			if(detectStart)
			{
				if((volt_print-preVolt)/dt > 0.1)
				{
					detectStart = false;
					detectPeak = true;
					startTime = ttt;
					restPotential = preVolt;
				}
			}
			if(detectPeak)
			{
				if(volt_print > -10.0 && volt_print <= preVolt)
				{
					detectPeak = false;
					detectEnd = true;
					v90 = volt_print - (volt_print - restPotential)*0.9;
				}
			}
			if(detectEnd)
			{
				if(volt_print <= v90)
				{
					detectEnd = false;
					detectStart = true;
					fprintf(apdrecord, "%.1lf\n", ttt - startTime);
				}
			}

			preVolt = volt_print;
		}

		if(iter%400 == 0)
		{
			printf("iter: %6d     time(ms): %6.1lf     volt: %6.2lf\n", iter, ttt, volt_print);
			if(pacingFinish && volt_print < -80.0) break;
		}
	}
	
	fclose(vmrecord);
	fclose(apdrecord);

	return 0;
}