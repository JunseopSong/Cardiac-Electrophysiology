//------------------------------------------------------------------------------
// State variables
//------------------------------------------------------------------------------

#define _NB_OF_STATE_VARIABLES_ 21
#define _NB_OF_STATE_VARIABLES_Fibroblast_ 19

extern double Y[_NB_OF_STATE_VARIABLES_];
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

extern double YY[_NB_OF_STATE_VARIABLES_Fibroblast_]; // Fibroblast
extern double dYY[_NB_OF_STATE_VARIABLES_Fibroblast_]; // Fibroblast

//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------

// Fibroblast
extern double K_o_;
extern double g_Kv_;
extern double P_NaK_;
extern double K_mK_;
extern double K_mNa_;
extern double V_rev_;
extern double G_Nab_;
extern double Na_o_;
extern double C_m_;
extern double beta_sf_; // 1.0 for fibroblast, 8.0 for myofibroblast
extern double G_gap_;
extern double S_mf_;
extern double S_m_;
extern double V_C_; // Cytoplasmic volume (micrometer^3)
extern double g_K1_;

extern double V_m_;
extern double I_f_to_m_;
extern double E_K_;
extern double K_i_; //129.4349
extern double I_Kv_;
extern double r_; //0
extern double s_; //1
extern double r_inf_;
extern double tau_r_;
extern double s_inf_;
extern double tau_s_;
extern double alpha_K1_;
extern double beta_K1_;
extern double I_K1_;
extern double I_NaK_;
extern double Na_i_; //8.5547
extern double E_Na_;
extern double I_Nab_;
extern double I_ion_;


extern double CMDN_max;   // millimolar (in Ca_buffers)
extern double CSQN_max;   // millimolar (in Ca_buffers)
extern double Km_CMDN;   // millimolar (in Ca_buffers)
extern double Km_CSQN;   // millimolar (in Ca_buffers)
extern double Km_TRPN;   // millimolar (in Ca_buffers)
extern double TRPN_max;   // millimolar (in Ca_buffers)
extern double Ca_up_max;   // millimolar (in Ca_leak_current_by_the_NSR)
extern double K_rel;   // per_millisecond (in Ca_release_current_from_JSR)
extern double I_up_max;   // millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)
extern double K_up;   // millimolar (in Ca_uptake_current_by_the_NSR)
extern double g_Ca_L;   // nanoS_per_picoF (in L_type_Ca_channel)
extern double I_NaCa_max;   // picoA_per_picoF (in Na_Ca_exchanger_current)
extern double K_mCa;   // millimolar (in Na_Ca_exchanger_current)
extern double K_mNa;   // millimolar (in Na_Ca_exchanger_current)
extern double K_sat;   // dimensionless (in Na_Ca_exchanger_current)
extern double _gamma;   // dimensionless (in Na_Ca_exchanger_current)
extern double g_B_Ca;   // nanoS_per_picoF (in background_currents)
extern double g_B_K;   // nanoS_per_picoF (in background_currents)
extern double g_B_Na;   // nanoS_per_picoF (in background_currents)
extern double g_Na;   // nanoS_per_picoF (in fast_sodium_current)
extern double V_cell;   // micrometre_3 (in intracellular_ion_concentrations)
extern double Cm;   // picoF (in membrane)
extern double F;   // coulomb_per_millimole (in membrane)
extern double R;   // joule_per_mole_kelvin (in membrane)
extern double T;   // kelvin (in membrane)
extern double g_Kr;   // nanoS_per_picoF (in rapid_delayed_rectifier_K_current)
extern double i_CaP_max;   // picoA_per_picoF (in sarcolemmal_calcium_pump_current)
extern double g_Ks;   // nanoS_per_picoF (in slow_delayed_rectifier_K_current)
extern double Km_K_o;   // millimolar (in sodium_potassium_pump)
extern double Km_Na_i;   // millimolar (in sodium_potassium_pump)
extern double i_NaK_max;   // picoA_per_picoF (in sodium_potassium_pump)
extern double Ca_o;   // millimolar (in standard_ionic_concentrations)
extern double K_o;   // millimolar (in standard_ionic_concentrations)
extern double Na_o;   // millimolar (in standard_ionic_concentrations)
extern double g_K1;   // nanoS_per_picoF (in time_independent_potassium_current)
extern double tau_tr;   // millisecond (in transfer_current_from_NSR_to_JSR)
extern double K_Q10;   // dimensionless (in transient_outward_K_current)
extern double g_to;   // nanoS_per_picoF (in transient_outward_K_current)
extern double stim_amplitude;   // picoA (in membrane)

//------------------------------------------------------------------------------
// Computed variables
//------------------------------------------------------------------------------

extern double Ca_CMDN;   // millimolar (in Ca_buffers)
extern double Ca_CSQN;   // millimolar (in Ca_buffers)
extern double Ca_TRPN;   // millimolar (in Ca_buffers)
extern double i_up_leak;   // millimolar_per_millisecond (in Ca_leak_current_by_the_NSR)
extern double tau_u;   // millisecond (in Ca_release_current_from_JSR_u_gate)
extern double u_infinity;   // dimensionless (in Ca_release_current_from_JSR_u_gate)
extern double tau_v;   // millisecond (in Ca_release_current_from_JSR_v_gate)
extern double v_infinity;   // dimensionless (in Ca_release_current_from_JSR_v_gate)
extern double tau_w;   // millisecond (in Ca_release_current_from_JSR_w_gate)
extern double w_infinity;   // dimensionless (in Ca_release_current_from_JSR_w_gate)
extern double Fn;   // dimensionless (in Ca_release_current_from_JSR)
extern double i_rel;   // millimolar_per_millisecond (in Ca_release_current_from_JSR)
extern double i_up;   // millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)
extern double d_infinity;   // dimensionless (in L_type_Ca_channel_d_gate)
extern double tau_d;   // millisecond (in L_type_Ca_channel_d_gate)
extern double f_Ca_infinity;   // dimensionless (in L_type_Ca_channel_f_Ca_gate)
extern double tau_f_Ca;   // millisecond (in L_type_Ca_channel_f_Ca_gate)
extern double f_infinity;   // dimensionless (in L_type_Ca_channel_f_gate)
extern double tau_f;   // millisecond (in L_type_Ca_channel_f_gate)
extern double i_Ca_L;   // picoA (in L_type_Ca_channel)
extern double i_NaCa;   // picoA (in Na_Ca_exchanger_current)
extern double E_Ca;   // millivolt (in background_currents)
extern double i_B_Ca;   // picoA (in background_currents)
extern double i_B_K;   // picoA (in background_currents)
extern double i_B_Na;   // picoA (in background_currents)
extern double alpha_h;   // per_millisecond (in fast_sodium_current_h_gate)
extern double beta_h;   // per_millisecond (in fast_sodium_current_h_gate)
extern double h_inf;   // dimensionless (in fast_sodium_current_h_gate)
extern double tau_h;   // millisecond (in fast_sodium_current_h_gate)
extern double alpha_j;   // per_millisecond (in fast_sodium_current_j_gate)
extern double beta_j;   // per_millisecond (in fast_sodium_current_j_gate)
extern double j_inf;   // dimensionless (in fast_sodium_current_j_gate)
extern double tau_j;   // millisecond (in fast_sodium_current_j_gate)
extern double alpha_m;   // per_millisecond (in fast_sodium_current_m_gate)
extern double beta_m;   // per_millisecond (in fast_sodium_current_m_gate)
extern double m_inf;   // dimensionless (in fast_sodium_current_m_gate)
extern double tau_m;   // millisecond (in fast_sodium_current_m_gate)
extern double E_Na;   // millivolt (in fast_sodium_current)
extern double i_Na;   // picoA (in fast_sodium_current)
extern double B1;   // millimolar_per_millisecond (in intracellular_ion_concentrations)
extern double B2;   // dimensionless (in intracellular_ion_concentrations)
extern double V_i;   // micrometre_3 (in intracellular_ion_concentrations)
extern double V_rel;   // micrometre_3 (in intracellular_ion_concentrations)
extern double V_up;   // micrometre_3 (in intracellular_ion_concentrations)
extern double i_st;   // picoA (in membrane)
extern double alpha_xr;   // per_millisecond (in rapid_delayed_rectifier_K_current_xr_gate)
extern double beta_xr;   // per_millisecond (in rapid_delayed_rectifier_K_current_xr_gate)
extern double tau_xr;   // millisecond (in rapid_delayed_rectifier_K_current_xr_gate)
extern double xr_infinity;   // dimensionless (in rapid_delayed_rectifier_K_current_xr_gate)
extern double i_Kr;   // picoA (in rapid_delayed_rectifier_K_current)
extern double i_CaP;   // picoA (in sarcolemmal_calcium_pump_current)
extern double alpha_xs;   // per_millisecond (in slow_delayed_rectifier_K_current_xs_gate)
extern double beta_xs;   // per_millisecond (in slow_delayed_rectifier_K_current_xs_gate)
extern double tau_xs;   // millisecond (in slow_delayed_rectifier_K_current_xs_gate)
extern double xs_infinity;   // dimensionless (in slow_delayed_rectifier_K_current_xs_gate)
extern double i_Ks;   // picoA (in slow_delayed_rectifier_K_current)
extern double f_NaK;   // dimensionless (in sodium_potassium_pump)
extern double i_NaK;   // picoA (in sodium_potassium_pump)
extern double sigma;   // dimensionless (in sodium_potassium_pump)
extern double E_K;   // millivolt (in time_independent_potassium_current)
extern double i_K1;   // picoA (in time_independent_potassium_current)
extern double i_tr;   // millimolar_per_millisecond (in transfer_current_from_NSR_to_JSR)
extern double alpha_oa;   // per_millisecond (in transient_outward_K_current_oa_gate)
extern double beta_oa;   // per_millisecond (in transient_outward_K_current_oa_gate)
extern double oa_infinity;   // dimensionless (in transient_outward_K_current_oa_gate)
extern double tau_oa;   // millisecond (in transient_outward_K_current_oa_gate)
extern double alpha_oi;   // per_millisecond (in transient_outward_K_current_oi_gate)
extern double beta_oi;   // per_millisecond (in transient_outward_K_current_oi_gate)
extern double oi_infinity;   // dimensionless (in transient_outward_K_current_oi_gate)
extern double tau_oi;   // millisecond (in transient_outward_K_current_oi_gate)
extern double i_to;   // picoA (in transient_outward_K_current)
extern double alpha_ua;   // per_millisecond (in ultrarapid_delayed_rectifier_K_current_ua_gate)
extern double beta_ua;   // per_millisecond (in ultrarapid_delayed_rectifier_K_current_ua_gate)
extern double tau_ua;   // millisecond (in ultrarapid_delayed_rectifier_K_current_ua_gate)
extern double ua_infinity;   // dimensionless (in ultrarapid_delayed_rectifier_K_current_ua_gate)
extern double alpha_ui;   // per_millisecond (in ultrarapid_delayed_rectifier_K_current_ui_gate)
extern double beta_ui;   // per_millisecond (in ultrarapid_delayed_rectifier_K_current_ui_gate)
extern double tau_ui;   // millisecond (in ultrarapid_delayed_rectifier_K_current_ui_gate)
extern double ui_infinity;   // dimensionless (in ultrarapid_delayed_rectifier_K_current_ui_gate)
extern double g_Kur;   // nanoS_per_picoF (in ultrarapid_delayed_rectifier_K_current)
extern double i_Kur;   // picoA (in ultrarapid_delayed_rectifier_K_current)

//------------------------------------------------------------------------------
// Procedures
//------------------------------------------------------------------------------

extern void CRN_compute(double dTime, int stim);
extern void Fibroblast_compute(double dTime, int stim, int Nf);
extern void CRN_Fib_init();
extern void CRN_Fib_compute(double dTime, int stim, int Nf);

//------------------------------------------------------------------------------

//==============================================================================
// End of file
//==============================================================================
