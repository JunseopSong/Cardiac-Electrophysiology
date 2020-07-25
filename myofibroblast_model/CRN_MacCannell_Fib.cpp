#include "CRN_MacCannell_Fib.h"
#include <math.h>
#include <stdio.h>

//------------------------------------------------------------------------------
// State variables
//------------------------------------------------------------------------------

double Y[_NB_OF_STATE_VARIABLES_];
double dY[_NB_OF_STATE_VARIABLES_];
// 0: u (dimensionless) (in Ca_release_current_from_JSR_u_gate)
// 1: v (dimensionless) (in Ca_release_current_from_JSR_v_gate)
// 2: w (dimensionless) (in Ca_release_current_from_JSR_w_gate)
// 3: d (dimensionless) (in L_type_Ca_channel_d_gate)
// 4: f_Ca (dimensionless) (in L_type_Ca_channel_f_Ca_gate)
// 5: f (dimensionless) (in L_type_Ca_channel_f_gate)
// 6: h (dimensionless) (in fast_sodium_current_h_gate)
// 7: j (dimensionless) (in fast_sodium_current_j_gate)
// 8: m (dimensionless) (in fast_sodium_current_m_gate)
// 9: Ca_i (millimolar) (in intracellular_ion_concentrations)
// 10: Ca_rel (millimolar) (in intracellular_ion_concentrations)
// 11: Ca_up (millimolar) (in intracellular_ion_concentrations)
// 12: K_i (millimolar) (in intracellular_ion_concentrations)
// 13: Na_i (millimolar) (in intracellular_ion_concentrations)
// 14: V (millivolt) (in membrane)
// 15: xr (dimensionless) (in rapid_delayed_rectifier_K_current_xr_gate)
// 16: xs (dimensionless) (in slow_delayed_rectifier_K_current_xs_gate)
// 17: oa (dimensionless) (in transient_outward_K_current_oa_gate)
// 18: oi (dimensionless) (in transient_outward_K_current_oi_gate)
// 19: ua (dimensionless) (in ultrarapid_delayed_rectifier_K_current_ua_gate)
// 20: ui (dimensionless) (in ultrarapid_delayed_rectifier_K_current_ui_gate)

double YY[_NB_OF_STATE_VARIABLES_Fibroblast_]; // Fibroblast

//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------

// Fibroblast
double K_o_;
double g_Kv_;
double P_NaK_;
double K_mK_;
double K_mNa_;
double V_rev_;
double G_Nab_; // nS/pF
double Na_o_;
double C_m_;
double beta_sf_; // 1.0 for fibroblast, 8.0 for myofibroblast
double G_gap_; // nS
double S_mf_;
double S_m_;
double V_C_;
double g_K1_;

double V_m_;
double I_f_to_m_;
double E_K_;
double K_i_; //129.4349
double I_Kv_;
double r_; //0
double s_; //1
double r_inf_;
double tau_r_;
double s_inf_;
double tau_s_;
double alpha_K1_;
double beta_K1_;
double I_K1_;
double I_NaK_;
double Na_i_; //8.5547
double E_Na_;
double I_Nab_;
double I_ion_;


double CMDN_max;   // millimolar (in Ca_buffers)
double CSQN_max;   // millimolar (in Ca_buffers)
double Km_CMDN;   // millimolar (in Ca_buffers)
double Km_CSQN;   // millimolar (in Ca_buffers)
double Km_TRPN;   // millimolar (in Ca_buffers)
double TRPN_max;   // millimolar (in Ca_buffers)
double Ca_up_max;   // millimolar (in Ca_leak_current_by_the_NSR)
double K_rel;   // per_millisecond (in Ca_release_current_from_JSR)
double I_up_max;   // millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)
double K_up;   // millimolar (in Ca_uptake_current_by_the_NSR)
double g_Ca_L;   // nanoS_per_picoF (in L_type_Ca_channel)
double I_NaCa_max;   // picoA_per_picoF (in Na_Ca_exchanger_current)
double K_mCa;   // millimolar (in Na_Ca_exchanger_current)
double K_mNa;   // millimolar (in Na_Ca_exchanger_current)
double K_sat;   // dimensionless (in Na_Ca_exchanger_current)
double _gamma;   // dimensionless (in Na_Ca_exchanger_current)
double g_B_Ca;   // nanoS_per_picoF (in background_currents)
double g_B_K;   // nanoS_per_picoF (in background_currents)
double g_B_Na;   // nanoS_per_picoF (in background_currents)
double g_Na;   // nanoS_per_picoF (in fast_sodium_current)
double V_cell;   // micrometre_3 (in intracellular_ion_concentrations)
double Cm;   // picoF (in membrane)
double F;   // coulomb_per_millimole (in membrane)
double R;   // joule_per_mole_kelvin (in membrane)
double T;   // kelvin (in membrane)
double g_Kr;   // nanoS_per_picoF (in rapid_delayed_rectifier_K_current)
double i_CaP_max;   // picoA_per_picoF (in sarcolemmal_calcium_pump_current)
double g_Ks;   // nanoS_per_picoF (in slow_delayed_rectifier_K_current)
double Km_K_o;   // millimolar (in sodium_potassium_pump)
double Km_Na_i;   // millimolar (in sodium_potassium_pump)
double i_NaK_max;   // picoA_per_picoF (in sodium_potassium_pump)
double Ca_o;   // millimolar (in standard_ionic_concentrations)
double K_o;   // millimolar (in standard_ionic_concentrations)
double Na_o;   // millimolar (in standard_ionic_concentrations)
double g_K1;   // nanoS_per_picoF (in time_independent_potassium_current)
double tau_tr;   // millisecond (in transfer_current_from_NSR_to_JSR)
double K_Q10;   // dimensionless (in transient_outward_K_current)
double g_to;   // nanoS_per_picoF (in transient_outward_K_current)
double stim_amplitude;   // picoA (in membrane)

//------------------------------------------------------------------------------
// Computed variables
//------------------------------------------------------------------------------

double Ca_CMDN;   // millimolar (in Ca_buffers)
double Ca_CSQN;   // millimolar (in Ca_buffers)
double Ca_TRPN;   // millimolar (in Ca_buffers)
double i_up_leak;   // millimolar_per_millisecond (in Ca_leak_current_by_the_NSR)
double tau_u;   // millisecond (in Ca_release_current_from_JSR_u_gate)
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
double tau_f_Ca;   // millisecond (in L_type_Ca_channel_f_Ca_gate)
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
double V_i;   // micrometre_3 (in intracellular_ion_concentrations)
double V_rel;   // micrometre_3 (in intracellular_ion_concentrations)
double V_up;   // micrometre_3 (in intracellular_ion_concentrations)
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
double sigma;   // dimensionless (in sodium_potassium_pump)
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

//------------------------------------------------------------------------------
// Initialisation
//------------------------------------------------------------------------------

void CRN_Fib_init()
{
	//---------------------------------------------------------------------------
	// State variables
	//---------------------------------------------------------------------------

	Y[0] = 2.35e-112;   // u (dimensionless) (in Ca_release_current_from_JSR_u_gate)
	Y[1] = 1.0;   // v (dimensionless) (in Ca_release_current_from_JSR_v_gate)
	Y[2] = 0.9992;   // w (dimensionless) (in Ca_release_current_from_JSR_w_gate)
	Y[3] = 1.367e-4;   // d (dimensionless) (in L_type_Ca_channel_d_gate)
	Y[4] = 7.755e-1;   // f_Ca (dimensionless) (in L_type_Ca_channel_f_Ca_gate)
	Y[5] = 9.996e-1;   // f (dimensionless) (in L_type_Ca_channel_f_gate)
	Y[6] = 9.649e-1;   // h (dimensionless) (in fast_sodium_current_h_gate)
	Y[7] = 9.775e-1;   // j (dimensionless) (in fast_sodium_current_j_gate)
	Y[8] = 2.908e-3;   // m (dimensionless) (in fast_sodium_current_m_gate)
	Y[9] = 1.013e-4;   // Ca_i (millimolar) (in intracellular_ion_concentrations)
	Y[10] = 1.488;   // Ca_rel (millimolar) (in intracellular_ion_concentrations)
	Y[11] = 1.488;   // Ca_up (millimolar) (in intracellular_ion_concentrations)
	Y[12] = 1.39e2;   // K_i (millimolar) (in intracellular_ion_concentrations)
	Y[13] = 1.117e1;   // Na_i (millimolar) (in intracellular_ion_concentrations)
	Y[14] = -81.18;   // V (millivolt) (in membrane)
	Y[15] = 3.296e-5;   // xr (dimensionless) (in rapid_delayed_rectifier_K_current_xr_gate)
	Y[16] = 1.869e-2;   // xs (dimensionless) (in slow_delayed_rectifier_K_current_xs_gate)
	Y[17] = 3.043e-2;   // oa (dimensionless) (in transient_outward_K_current_oa_gate)
	Y[18] = 9.992e-1;   // oi (dimensionless) (in transient_outward_K_current_oi_gate)
	Y[19] = 4.966e-3;   // ua (dimensionless) (in ultrarapid_delayed_rectifier_K_current_ua_gate)
	Y[20] = 9.986e-1;   // ui (dimensionless) (in ultrarapid_delayed_rectifier_K_current_ui_gate)

	//---------------------------------------------------------------------------
	// Constants
	//---------------------------------------------------------------------------

	// Fibroblast
	K_o_ = 5.4; // mmol/L
	g_Kv_ = 0.25; // nS/pF
	P_NaK_ = 2.002; // pA/pF
	K_mK_ = 1.0; // mmol/L
	K_mNa_ = 11.0; // mmol/L
	V_rev_ = -150.0; // mV
	G_Nab_ = 9.5e-3; // nS/pF
	Na_o_ = 140.0; // mmol/L
	C_m_ = 6.3; // pF
	beta_sf_ = 1.0; // 1.0 for fibroblast, 8.0 for myofibroblast
	G_gap_ = 3.0; // nS
	S_mf_ = 6.3e-6;
	S_m_ = 1.0e-4;
    g_K1_ = 0.4822; // nS/pF

	YY[0] = -49.6; // mV
	YY[3] = 129.4349; // mmol/L
	YY[5] = 0.0; // dimensionless
	YY[6] = 1.0; // dimensionless
	YY[16] = 8.5547; // mmol/L


	CMDN_max = 0.05;   // millimolar (in Ca_buffers)
	CSQN_max = 10.0;   // millimolar (in Ca_buffers)
	Km_CMDN = 0.00238;   // millimolar (in Ca_buffers)
	Km_CSQN = 0.8;   // millimolar (in Ca_buffers)
	Km_TRPN = 0.0005;   // millimolar (in Ca_buffers)
	TRPN_max = 0.07;   // millimolar (in Ca_buffers)
	Ca_up_max = 15.0;   // millimolar (in Ca_leak_current_by_the_NSR)
	K_rel = 30.0;   // per_millisecond (in Ca_release_current_from_JSR)
	I_up_max = 0.005;   // millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)
	K_up = 0.00092;   // millimolar (in Ca_uptake_current_by_the_NSR)
	g_Ca_L = 0.12375;   // nanoS_per_picoF (in L_type_Ca_channel)
	I_NaCa_max = 1600.0;   // picoA_per_picoF (in Na_Ca_exchanger_current)
	K_mCa = 1.38;   // millimolar (in Na_Ca_exchanger_current)
	K_mNa = 87.5;   // millimolar (in Na_Ca_exchanger_current)
	K_sat = 0.1;   // dimensionless (in Na_Ca_exchanger_current)
	_gamma = 0.35;   // dimensionless (in Na_Ca_exchanger_current)
	g_B_Ca = 0.001131;   // nanoS_per_picoF (in background_currents)
	g_B_K = 0.0;   // nanoS_per_picoF (in background_currents)
	g_B_Na = 0.0006744375;   // nanoS_per_picoF (in background_currents)
	g_Na = 7.8;   // nanoS_per_picoF (in fast_sodium_current)
	V_cell = 20100.0;   // micrometre_3 (in intracellular_ion_concentrations)
	Cm = 100.0;   // picoF (in membrane)
	F = 96.4867;   // coulomb_per_millimole (in membrane)
	R = 8.3143;   // joule_per_mole_kelvin (in membrane)
	T = 310.0;   // kelvin (in membrane)
	g_Kr = 0.029411765;   // nanoS_per_picoF (in rapid_delayed_rectifier_K_current)
	i_CaP_max = 0.275;   // picoA_per_picoF (in sarcolemmal_calcium_pump_current)
	g_Ks = 0.12941176;   // nanoS_per_picoF (in slow_delayed_rectifier_K_current)
	Km_K_o = 1.5;   // millimolar (in sodium_potassium_pump)
	Km_Na_i = 10.0;   // millimolar (in sodium_potassium_pump)
	i_NaK_max = 0.59933874;   // picoA_per_picoF (in sodium_potassium_pump)
	Ca_o = 1.8;   // millimolar (in standard_ionic_concentrations)
	K_o = 5.4;   // millimolar (in standard_ionic_concentrations)
	Na_o = 140.0;   // millimolar (in standard_ionic_concentrations)
	g_K1 = 0.09;   // nanoS_per_picoF (in time_independent_potassium_current)
	tau_tr = 180.0;   // millisecond (in transfer_current_from_NSR_to_JSR)
	K_Q10 = 3.0;   // dimensionless (in transient_outward_K_current)
	g_to = 0.1652;   // nanoS_per_picoF (in transient_outward_K_current)
	stim_amplitude = -1900.0;   // picoA (in membrane)

	//---------------------------------------------------------------------------
	// Computed variables
	//---------------------------------------------------------------------------

	V_rel = 0.0048*V_cell;
	tau_u = 8.0;
	tau_f_Ca = 2.0;
	V_i = V_cell*0.68;
	V_up = 0.0552*V_cell;
	sigma = 1.0/7.0*(exp(Na_o/67.3)-1.0);

	V_C_ = 1370; // micrometer^3
}

//------------------------------------------------------------------------------
// Computation
//------------------------------------------------------------------------------

/*
Fibroblast (MacCannell model)
Based on Ashihara et al. (2012) Supplemental Material
*/
void Fibroblast_compute(double dTime, int stim, int Nf)
{
    // I : pA/pF
	V_m_ = YY[0];
	I_f_to_m_ = YY[1];
	E_K_ = YY[2];
	K_i_ = YY[3];
	I_Kv_ = YY[4];
	r_ = YY[5];
	s_ = YY[6];
	r_inf_ = YY[7];
	tau_r_ = YY[8];
	s_inf_ = YY[9];
	tau_s_ = YY[10];
	alpha_K1_ = YY[11];
	beta_K1_ = YY[12];
	I_K1_ = YY[13];
	I_ion_ = YY[14];
	I_NaK_ = YY[15];
	Na_i_ = YY[16];
	E_Na_ = YY[17];
	I_Nab_ = YY[18];

	// Voltage-dependent K+ current (I_Kv_)
	tau_r_ = 20.3 + 138.0*exp(-pow((V_m_+20.0)/25.9, 2.0));
	tau_s_ = 1574.0 + 5268.0*exp(-pow((V_m_+23.0)/22.7, 2.0));
	r_inf_ = 1.0/(1.0+exp(-(V_m_+20.0)/11.0));
	s_inf_ = 1.0/(1.0+exp((V_m_+23.0)/7.0));
	r_ = r_inf_ + (r_-r_inf_)*exp(-dTime/tau_r_);
	s_ = s_inf_ + (s_-s_inf_)*exp(-dTime/tau_s_);
	E_K_ = log(K_o_/K_i_)*(R*T/F);
	I_Kv_ = g_Kv_*r_*s_*(V_m_-E_K_);

	// Inward-rectifying K+ current (I_K1_)
	alpha_K1_ = 0.1/(1.0+exp(0.06*(V_m_-E_K_-200.0)));
	beta_K1_ = (3.0*exp(2.0e-4*(V_m_-E_K_+100.0)) + exp(0.1*(V_m_-E_K_-10.0)))/(1.0 + exp(-0.5*(V_m_-E_K_)));
	I_K1_ = g_K1_*(alpha_K1_/(alpha_K1_+beta_K1_))*(V_m_-E_K_);

	// Na+ - K+ pump current (I_NaK_)
	I_NaK_ = P_NaK_*(K_o_/(K_o+K_mK_))*(pow(Na_i_, 1.5)/(pow(Na_i_, 1.5)+pow(K_mNa_, 1.5)))*(V_m_-V_rev_)/(V_m_+200.0);

	// Na+ backgound current (I_Nab_)
	E_Na_ = log(Na_o_/Na_i_)*(R*T/F);
	I_Nab_ = G_Nab_*(V_m_-E_Na_);

	// Ion concentration
	Na_i_ -= C_m_*(3.0*I_NaK_ + I_Nab_)/(V_C_*F)*dTime;
	K_i_ -= C_m_*(I_Kv_ + I_K1_ - 2.0*I_NaK_)/(V_C_*F)*dTime;

	// Coupling
	I_ion_ = I_Kv_ + I_K1_ + I_NaK_ + I_Nab_;

	V_m_ -= (I_ion_ + G_gap_*(V_m_-Y[14])/C_m_)/beta_sf_*dTime;

	I_f_to_m_ = G_gap_*(Y[14]-V_m_)*Nf;

	YY[0] = V_m_;
	YY[1] = I_f_to_m_;
	YY[2] = E_K_;
	YY[3] = K_i_;
	YY[4] = I_Kv_;
	YY[5] = r_;
	YY[6] = s_;
	YY[7] = r_inf_;
	YY[8] = tau_r_;
	YY[9] = s_inf_;
	YY[10] = tau_s_;
	YY[11] = alpha_K1_;
	YY[12] = beta_K1_;
	YY[13] = I_K1_;
	YY[14] = I_ion_;
	YY[15] = I_NaK_;
	YY[16] = Na_i_;
	YY[17] = E_Na_;
	YY[18] = I_Nab_;
}

// Atrial myocyte (Courtemanche model)
void CRN_compute(double dTime, int stim)
{
	Ca_CMDN = CMDN_max*Y[9]/(Y[9]+Km_CMDN);
	Ca_TRPN = TRPN_max*Y[9]/(Y[9]+Km_TRPN);
	Ca_CSQN = CSQN_max*Y[10]/(Y[10]+Km_CSQN);
	i_up_leak = I_up_max*Y[11]/Ca_up_max;
	i_rel = K_rel*pow(Y[0], 2.0)*Y[1]*Y[2]*(Y[10]-Y[9]);
	i_Ca_L = Cm*g_Ca_L*Y[3]*Y[5]*Y[4]*(Y[14]-65.0);
	i_NaCa = Cm*I_NaCa_max*(exp(_gamma*F*Y[14]/(R*T))*pow(Y[13], 3.0)*Ca_o-exp((_gamma-1.0)*F*Y[14]/(R*T))*pow(Na_o, 3.0)*Y[9])/((pow(K_mNa, 3.0)+pow(Na_o, 3.0))*(K_mCa+Ca_o)*(1.0+K_sat*exp((_gamma-1.0)*Y[14]*F/(R*T))));
	Fn = 1.0e3*(1.0e-15*V_rel*i_rel-1.0e-15/(2.0*F)*(0.5*i_Ca_L-0.2*i_NaCa));
	u_infinity = pow(1.0+exp(-(Fn-3.4175e-13)/13.67e-16), -1.0);
	Y[0] = u_infinity+(Y[0]-u_infinity)*exp(-dTime/tau_u); // dY[0] = (u_infinity-Y[0])/tau_u;
	tau_v = 1.91+2.09*pow(1.0+exp(-(Fn-3.4175e-13)/13.67e-16), -1.0);
	v_infinity = 1.0-pow(1.0+exp(-(Fn-6.835e-14)/13.67e-16), -1.0);
	Y[1] = v_infinity+(Y[1]-v_infinity)*exp(-dTime/tau_v); // dY[1] = (v_infinity-Y[1])/tau_v;

	if (fabs(Y[14]-7.9) < 1.0e-10)
		tau_w = 6.0*0.2/1.3;
	else
		tau_w = 6.0*(1.0-exp(-(Y[14]-7.9)/5.0))/((1.0+0.3*exp(-(Y[14]-7.9)/5.0))*1.0*(Y[14]-7.9));

	w_infinity = 1.0-pow(1.0+exp(-(Y[14]-40.0)/17.0), -1.0);
	Y[2] = w_infinity+(Y[2]-w_infinity)*exp(-dTime/tau_w); // dY[2] = (w_infinity-Y[2])/tau_w;
	i_up = I_up_max/(1.0+K_up/Y[9]);
	d_infinity = pow(1.0+exp((Y[14]+10.0)/-8.0), -1.0);

	if (fabs(Y[14]+10.0) < 1.0e-10)
		tau_d = 4.579/(1.0+exp((Y[14]+10.0)/-6.24));
	else
		tau_d = (1.0-exp((Y[14]+10.0)/-6.24))/(0.035*(Y[14]+10.0)*(1.0+exp((Y[14]+10.0)/-6.24)));

	Y[3] = d_infinity+(Y[3]-d_infinity)*exp(-dTime/tau_d); // dY[3] = (d_infinity-Y[3])/tau_d;
	f_Ca_infinity = pow(1.0+Y[9]/0.00035, -1.0);
	Y[4] = f_Ca_infinity+(Y[4]-f_Ca_infinity)*exp(-dTime/tau_f_Ca); // dY[4] = (f_Ca_infinity-Y[4])/tau_f_Ca;
	f_infinity = exp(-(Y[14]+28.0)/6.9)/(1.0+exp(-(Y[14]+28.0)/6.9));
	tau_f = 9.0*pow(0.0197*exp(-pow(0.0337, 2.0)*pow(Y[14]+10.0, 2.0))+0.02, -1.0);
	Y[5] = f_infinity+(Y[5]-f_infinity)*exp(-dTime/tau_f); // dY[5] = (f_infinity-Y[5])/tau_f;
	E_Ca = R*T/(2.0*F)*log(Ca_o/Y[9]);
	E_Na = R*T/F*log(Na_o/Y[13]);
	i_B_Na = Cm*g_B_Na*(Y[14]-E_Na);
	i_B_Ca = Cm*g_B_Ca*(Y[14]-E_Ca);
	E_K = R*T/F*log(K_o/Y[12]);
	i_B_K = Cm*g_B_K*(Y[14]-E_K);
	i_Na = Cm*g_Na*pow(Y[8], 3.0)*Y[6]*Y[7]*(Y[14]-E_Na);

	if (Y[14] < -40.0)
		alpha_h = 0.135*exp((Y[14]+80.0)/-6.8);
	else
		alpha_h = 0.0;

	if (Y[14] < -40.0)
		beta_h = 3.56*exp(0.079*Y[14])+3.1e5*exp(0.35*Y[14]);
	else
		beta_h = 1.0/(0.13*(1.0+exp((Y[14]+10.66)/-11.1)));

	h_inf = alpha_h/(alpha_h+beta_h);
	tau_h = 1.0/(alpha_h+beta_h);
	Y[6] = h_inf+(Y[6]-h_inf)*exp(-dTime/tau_h); // dY[6] = (h_inf-Y[6])/tau_h;

	if (Y[14] < -40.0)
		alpha_j = (-1.2714e5*exp(0.2444*Y[14])-3.474e-5*exp(-0.04391*Y[14]))*(Y[14]+37.78)/(1.0+exp(0.311*(Y[14]+79.23)));
	else
		alpha_j = 0.0;

	if (Y[14] < -40.0)
		beta_j = 0.1212*exp(-0.01052*Y[14])/(1.0+exp(-0.1378*(Y[14]+40.14)));
	else
		beta_j = 0.3*exp(-2.535e-7*Y[14])/(1.0+exp(-0.1*(Y[14]+32.0)));

	j_inf = alpha_j/(alpha_j+beta_j);
	tau_j = 1.0/(alpha_j+beta_j);
	Y[7] = j_inf+(Y[7]-j_inf)*exp(-dTime/tau_j); // dY[7] = (j_inf-Y[7])/tau_j;

	if (Y[14] == -47.13)
		alpha_m = 3.2;
	else
		alpha_m = 0.32*(Y[14]+47.13)/(1.0-exp(-0.1*(Y[14]+47.13)));

	beta_m = 0.08*exp(-Y[14]/11.0);
	m_inf = alpha_m/(alpha_m+beta_m);
	tau_m = 1.0/(alpha_m+beta_m);
	Y[8] = m_inf+(Y[8]-m_inf)*exp(-dTime/tau_m); // dY[8] = (m_inf-Y[8])/tau_m;
	f_NaK = pow(1.0+0.1245*exp(-0.1*F*Y[14]/(R*T))+0.0365*sigma*exp(-F*Y[14]/(R*T)), -1.0);
	i_NaK = Cm*i_NaK_max*f_NaK*1.0/(1.0+pow(Km_Na_i/Y[13], 1.5))*K_o/(K_o+Km_K_o);
	dY[13] = (-3.0*i_NaK-(3.0*i_NaCa+i_B_Na+i_Na))/(V_i*F);
	i_K1 = Cm*g_K1*(Y[14]-E_K)/(1.0+exp(0.07*(Y[14]+80.0)));
	i_to = Cm*g_to*pow(Y[17], 3.0)*Y[18]*(Y[14]-E_K);
	g_Kur = 0.005+0.05/(1.0+exp((Y[14]-15.0)/-13.0));
	i_Kur = Cm*g_Kur*pow(Y[19], 3.0)*Y[20]*(Y[14]-E_K);
	i_Kr = Cm*g_Kr*Y[15]*(Y[14]-E_K)/(1.0+exp((Y[14]+15.0)/22.4));
	i_Ks = Cm*g_Ks*pow(Y[16], 2.0)*(Y[14]-E_K);
	dY[12] = (2.0*i_NaK-(i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_K))/(V_i*F);
	i_CaP = Cm*i_CaP_max*Y[9]/(0.0005+Y[9]);
	B1 = (2.0*i_NaCa-(i_CaP+i_Ca_L+i_B_Ca))/(2.0*V_i*F)+(V_up*(i_up_leak-i_up)+i_rel*V_rel)/V_i;
	B2 = 1.0+TRPN_max*Km_TRPN/pow(Y[9]+Km_TRPN, 2.0)+CMDN_max*Km_CMDN/pow(Y[9]+Km_CMDN, 2.0);
	dY[9] = B1/B2;
	i_tr = (Y[11]-Y[10])/tau_tr;
	dY[11] = i_up-(i_up_leak+i_tr*V_rel/V_up);
	dY[10] = (i_tr-i_rel)*pow(1.0+CSQN_max*Km_CSQN/pow(Y[10]+Km_CSQN, 2.0), -1.0);

	// Stimulation
	if (stim == 1)
		i_st = stim_amplitude;
	else
		i_st = 0.0;

	dY[14] = -(i_Na+i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_Na+i_B_Ca+i_NaK+i_CaP+i_NaCa+i_Ca_L+i_st+YY[1])/Cm;

	if (fabs(Y[14]+14.1) < 1.0e-10)
		alpha_xr = 0.0015;
	else
		alpha_xr = 0.0003*(Y[14]+14.1)/(1.0-exp((Y[14]+14.1)/-5.0));

	if (fabs(Y[14]-3.3328) < 1.0e-10)
		beta_xr = 3.7836118e-4;
	else
		beta_xr = 0.000073898*(Y[14]-3.3328)/(exp((Y[14]-3.3328)/5.1237)-1.0);

	tau_xr = pow(alpha_xr+beta_xr, -1.0);
	xr_infinity = pow(1.0+exp((Y[14]+14.1)/-6.5), -1.0);
	Y[15] = xr_infinity+(Y[15]-xr_infinity)*exp(-dTime/tau_xr); // dY[15] = (xr_infinity-Y[15])/tau_xr;

	if (fabs(Y[14]-19.9) < 1.0e-10)
		alpha_xs = 0.00068;
	else
		alpha_xs = 0.00004*(Y[14]-19.9)/(1.0-exp((Y[14]-19.9)/-17.0));

	if (fabs(Y[14]-19.9) < 1.0e-10)
		beta_xs = 0.000315;
	else
		beta_xs = 0.000035*(Y[14]-19.9)/(exp((Y[14]-19.9)/9.0)-1.0);

	tau_xs = 0.5*pow(alpha_xs+beta_xs, -1.0);
	xs_infinity = pow(1.0+exp((Y[14]-19.9)/-12.7), -0.5);
	Y[16] = xs_infinity+(Y[16]-xs_infinity)*exp(-dTime/tau_xs); // dY[16] = (xs_infinity-Y[16])/tau_xs;
	alpha_oa = 0.65*pow(exp((Y[14]-(-10.0))/-8.5)+exp((Y[14]-(-10.0)-40.0)/-59.0), -1.0);
	beta_oa = 0.65*pow(2.5+exp((Y[14]-(-10.0)+72.0)/17.0), -1.0);
	tau_oa = pow(alpha_oa+beta_oa, -1.0)/K_Q10;
	oa_infinity = pow(1.0+exp((Y[14]-(-10.0)+10.47)/-17.54), -1.0);
	Y[17] = oa_infinity+(Y[17]-oa_infinity)*exp(-dTime/tau_oa); // dY[17] = (oa_infinity-Y[17])/tau_oa;
	alpha_oi = pow(18.53+1.0*exp((Y[14]-(-10.0)+103.7)/10.95), -1.0);
	beta_oi = pow(35.56+1.0*exp((Y[14]-(-10.0)-8.74)/-7.44), -1.0);
	tau_oi = pow(alpha_oi+beta_oi, -1.0)/K_Q10;
	oi_infinity = pow(1.0+exp((Y[14]-(-10.0)+33.1)/5.3), -1.0);
	Y[18] = oi_infinity+(Y[18]-oi_infinity)*exp(-dTime/tau_oi); // dY[18] = (oi_infinity-Y[18])/tau_oi;
	alpha_ua = 0.65*pow(exp((Y[14]-(-10.0))/-8.5)+exp((Y[14]-(-10.0)-40.0)/-59.0), -1.0);
	beta_ua = 0.65*pow(2.5+exp((Y[14]-(-10.0)+72.0)/17.0), -1.0);
	tau_ua = pow(alpha_ua+beta_ua, -1.0)/K_Q10;
	ua_infinity = pow(1.0+exp((Y[14]-(-10.0)+20.3)/-9.6), -1.0);
	Y[19] = ua_infinity+(Y[19]-ua_infinity)*exp(-dTime/tau_ua); // dY[19] = (ua_infinity-Y[19])/tau_ua;
	alpha_ui = pow(21.0+1.0*exp((Y[14]-(-10.0)-195.0)/-28.0), -1.0);
	beta_ui = 1.0/exp((Y[14]-(-10.0)-168.0)/-16.0);
	tau_ui = pow(alpha_ui+beta_ui, -1.0)/K_Q10;
	ui_infinity = pow(1.0+exp((Y[14]-(-10.0)-109.45)/27.48), -1.0);
	Y[20] = ui_infinity+(Y[20]-ui_infinity)*exp(-dTime/tau_ui); // dY[20] = (ui_infinity-Y[20])/tau_ui;

	Y[9] += dY[9]*dTime;
	Y[10] += dY[10]*dTime;
	Y[11] += dY[11]*dTime;
	Y[12] += dY[12]*dTime;
	Y[13] += dY[13]*dTime;
	Y[14] += dY[14]*dTime;
}

void CRN_Fib_compute(double dTime, int stim, int Nf)
{
	Fibroblast_compute(dTime, stim, Nf);
	CRN_compute(dTime, stim);
}

//==============================================================================
// End of file
//==============================================================================
