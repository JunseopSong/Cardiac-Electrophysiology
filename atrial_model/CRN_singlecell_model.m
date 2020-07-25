% Courtemanche-Ramirez-Nattel atrial cell model

function dY = CRN_singlecell_model(t, Y, I_stim)

global g_Na g_B_Na g_Ca_L g_B_Ca g_to g_K1 g_Kr g_Ks
global i_NaK_max I_NaCa_max K_mCa K_mNa gamma K_sat Ca_up_max I_up_max K_rel
global K_up i_CaP_max GKur

dY = zeros(size(Y));  % dY/dt

CMDN_max = 0.05;   % millimolar (in Ca_buffers)
CSQN_max = 10.0;   % millimolar (in Ca_buffers)
Km_CMDN = 0.00238;   % millimolar (in Ca_buffers)
Km_CSQN = 0.8;   % millimolar (in Ca_buffers)
Km_TRPN = 0.0005;   % millimolar (in Ca_buffers)
TRPN_max = 0.07;   % millimolar (in Ca_buffers)
%Ca_up_max = 15.0;   % millimolar (in Ca_leak_current_by_the_NSR)
%K_rel = 30.0;   % per_millisecond (in Ca_release_current_from_JSR)
%I_up_max = 0.005;   % millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)
%K_up = 0.00092;   % millimolar (in Ca_uptake_current_by_the_NSR)
%g_Ca_L = 0.12375;   % nanoS_per_picoF (in L_type_Ca_channel)
%I_NaCa_max = 1600.0;   % picoA_per_picoF (in Na_Ca_exchanger_current)
%K_mCa = 1.38;   % millimolar (in Na_Ca_exchanger_current)
%K_mNa = 87.5;   % millimolar (in Na_Ca_exchanger_current)
%K_sat = 0.1;   % dimensionless (in Na_Ca_exchanger_current)
%gamma = 0.35;   % dimensionless (in Na_Ca_exchanger_current)
%g_B_Ca = 0.001131;   % nanoS_per_picoF (in background_currents)
g_B_K = 0.0;   % nanoS_per_picoF (in background_currents)
%g_B_Na = 0.0006744375;   % nanoS_per_picoF (in background_currents)
%g_Na = 7.8;   % nanoS_per_picoF (in fast_sodium_current)
V_cell = 20100.0;   % micrometre_3 (in intracellular_ion_concentrations)
Cm = 100.0;   % picoF (in membrane)

F = 96.4867;   % coulomb_per_millimole (in membrane)
R = 8.3143;   % joule_per_mole_kelvin (in membrane)
T = 310.0;   % kelvin (in membrane)

%g_Kr = 0.029411765;   % nanoS_per_picoF (in rapid_delayed_rectifier_K_current)
%i_CaP_max = 0.275;   % picoA_per_picoF (in sarcolemmal_calcium_pump_current)
%g_Ks = 0.12941176;   % nanoS_per_picoF (in slow_delayed_rectifier_K_current)
Km_K_o = 1.5;   % millimolar (in sodium_potassium_pump)
Km_Na_i = 10.0;   % millimolar (in sodium_potassium_pump)
%i_NaK_max = 0.59933874;   % picoA_per_picoF (in sodium_potassium_pump)
Ca_o = 1.8;   % millimolar (in standard_ionic_concentrations)
K_o = 5.4;   % millimolar (in standard_ionic_concentrations)
Na_o = 140.0;   % millimolar (in standard_ionic_concentrations)
%g_K1 = 0.09;   % nanoS_per_picoF (in time_independent_potassium_current)
tau_tr = 180.0;   % millisecond (in transfer_current_from_NSR_to_JSR)
K_Q10 = 3.0;   % dimensionless (in transient_outward_K_current)
%g_to = 0.1652;   % nanoS_per_picoF (in transient_outward_K_current)

% time (millisecond)
Ca_CMDN = CMDN_max*Y(10)/(Y(10)+Km_CMDN);
Ca_TRPN = TRPN_max*Y(10)/(Y(10)+Km_TRPN);
Ca_CSQN = CSQN_max*Y(11)/(Y(11)+Km_CSQN);
i_up_leak = I_up_max*Y(12)/Ca_up_max;
V_rel = 0.0048*V_cell;
i_rel = K_rel*Y(1)^2.0*Y(2)*Y(3)*(Y(11)-Y(10));

i_Ca_L = Cm*g_Ca_L*Y(4)*Y(6)*Y(5)*(Y(15)-65.0);
i_NaCa = Cm*I_NaCa_max*(exp(gamma*F*Y(15)/(R*T))*Y(14)^3.0*Ca_o-exp((gamma-1.0)*F*Y(15)/(R*T))*Na_o^3.0*Y(10))/((K_mNa^3.0+Na_o^3.0)*(K_mCa+Ca_o)*(1.0+K_sat*exp((gamma-1.0)*Y(15)*F/(R*T))));
Fn = 1.0e3*(1.0e-15*V_rel*i_rel-1.0e-15/(2.0*F)*(0.5*i_Ca_L-0.2*i_NaCa));
tau_u = 8.0;
u_infinity = (1.0+exp(-(Fn-3.4175e-13)/13.67e-16))^-1.0;
dY(1) = (u_infinity-Y(1))/tau_u;
tau_v = 1.91+2.09*(1.0+exp(-(Fn-3.4175e-13)/13.67e-16))^-1.0;
v_infinity = 1.0-(1.0+exp(-(Fn-6.835e-14)/13.67e-16))^-1.0;
dY(2) = (v_infinity-Y(2))/tau_v;

if (abs(Y(15)-7.9) < 1.0e-10)
    tau_w = 6.0*0.2/1.3;
else
    tau_w = 6.0*(1.0-exp(-(Y(15)-7.9)/5.0))/((1.0+0.3*exp(-(Y(15)-7.9)/5.0))*1.0*(Y(15)-7.9));
end;

w_infinity = 1.0-(1.0+exp(-(Y(15)-40.0)/17.0))^-1.0;
dY(3) = (w_infinity-Y(3))/tau_w;
i_up = I_up_max/(1.0+K_up/Y(10));
d_infinity = (1.0+exp((Y(15)+10.0)/-8.0))^-1.0;

if (abs(Y(15)+10.0) < 1.0e-10)
    tau_d = 4.579/(1.0+exp((Y(15)+10.0)/-6.24));
else
    tau_d = (1.0-exp((Y(15)+10.0)/-6.24))/(0.035*(Y(15)+10.0)*(1.0+exp((Y(15)+10.0)/-6.24)));
end;

dY(4) = (d_infinity-Y(4))/tau_d;
f_Ca_infinity = (1.0+Y(10)/0.00035)^-1.0;
tau_f_Ca = 2.0;
dY(5) = (f_Ca_infinity-Y(5))/tau_f_Ca;
f_infinity = exp(-(Y(15)+28.0)/6.9)/(1.0+exp(-(Y(15)+28.0)/6.9));
tau_f = 9.0*(0.0197*exp(-0.0337^2.0*(Y(15)+10.0)^2.0)+0.02)^-1.0;
dY(6) = (f_infinity-Y(6))/tau_f;
E_Ca = R*T/(2.0*F)*log(Ca_o/Y(10));
E_Na = R*T/F*log(Na_o/Y(14));
i_B_Na = Cm*g_B_Na*(Y(15)-E_Na);
i_B_Ca = Cm*g_B_Ca*(Y(15)-E_Ca);
E_K = R*T/F*log(K_o/Y(13));
i_B_K = Cm*g_B_K*(Y(15)-E_K);
i_Na = Cm*g_Na*Y(9)^3.0*Y(7)*Y(8)*(Y(15)-E_Na);

if (Y(15) < -40.0)
    alpha_h = 0.135*exp((Y(15)+80.0)/-6.8);
else
    alpha_h = 0.0;
end;

if (Y(15) < -40.0)
    beta_h = 3.56*exp(0.079*Y(15))+3.1e5*exp(0.35*Y(15));
else
    beta_h = 1.0/(0.13*(1.0+exp((Y(15)+10.66)/-11.1)));
end;

h_inf = alpha_h/(alpha_h+beta_h);
tau_h = 1.0/(alpha_h+beta_h);
dY(7) = (h_inf-Y(7))/tau_h;

if (Y(15) < -40.0)
    alpha_j = (-1.2714e5*exp(0.2444*Y(15))-3.474e-5*exp(-0.04391*Y(15)))*(Y(15)+37.78)/(1.0+exp(0.311*(Y(15)+79.23)));
else
    alpha_j = 0.0;
end;

if (Y(15) < -40.0)
    beta_j = 0.1212*exp(-0.01052*Y(15))/(1.0+exp(-0.1378*(Y(15)+40.14)));
else
    beta_j = 0.3*exp(-2.535e-7*Y(15))/(1.0+exp(-0.1*(Y(15)+32.0)));
end;

j_inf = alpha_j/(alpha_j+beta_j);
tau_j = 1.0/(alpha_j+beta_j);
dY(8) = (j_inf-Y(8))/tau_j;

if (Y(15) == -47.13)
    alpha_m = 3.2;
else
    alpha_m = 0.32*(Y(15)+47.13)/(1.0-exp(-0.1*(Y(15)+47.13)));
end;

beta_m = 0.08*exp(-Y(15)/11.0);
m_inf = alpha_m/(alpha_m+beta_m);
tau_m = 1.0/(alpha_m+beta_m);
dY(9) = (m_inf-Y(9))/tau_m;
V_i = V_cell*0.68;
V_up = 0.0552*V_cell;
sigma = 1.0/7.0*(exp(Na_o/67.3)-1.0);
f_NaK = (1.0+0.1245*exp(-0.1*F*Y(15)/(R*T))+0.0365*sigma*exp(-F*Y(15)/(R*T)))^-1.0;
i_NaK = Cm*i_NaK_max*f_NaK*1.0/(1.0+(Km_Na_i/Y(14))^1.5)*K_o/(K_o+Km_K_o);


dY(14) = (-3.0*i_NaK-(3.0*i_NaCa+i_B_Na+i_Na))/(V_i*F);
i_K1 = Cm*g_K1*(Y(15)-E_K)/(1.0+exp(0.07*(Y(15)+80.0)));
i_to = Cm*g_to*Y(18)^3.0*Y(19)*(Y(15)-E_K);

%% IKur
g_Kur = 0.005+0.05/(1.0+exp((Y(15)-15.0)/-13.0));
i_Kur = Cm*GKur*g_Kur*Y(20)^3.0*Y(21)*(Y(15)-E_K);

i_Kr = Cm*g_Kr*Y(16)*(Y(15)-E_K)/(1.0+exp((Y(15)+15.0)/22.4));
i_Ks = Cm*g_Ks*Y(17)^2.0*(Y(15)-E_K);
dY(13) = (2.0*i_NaK-(i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_K))/(V_i*F);
i_CaP = Cm*i_CaP_max*Y(10)/(0.0005+Y(10));
B1 = (2.0*i_NaCa-(i_CaP+i_Ca_L+i_B_Ca))/(2.0*V_i*F)+(V_up*(i_up_leak-i_up)+i_rel*V_rel)/V_i;
B2 = 1.0+TRPN_max*Km_TRPN/(Y(10)+Km_TRPN)^2.0+CMDN_max*Km_CMDN/(Y(10)+Km_CMDN)^2.0;
dY(10) = B1/B2;
i_tr = (Y(12)-Y(11))/tau_tr;
dY(12) = i_up-(i_up_leak+i_tr*V_rel/V_up);
dY(11) = (i_tr-i_rel)*(1.0+CSQN_max*Km_CSQN/(Y(11)+Km_CSQN)^2.0)^-1.0;

dY(15) = -(i_Na+i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_Na+i_B_Ca+i_NaK+i_CaP+i_NaCa+i_Ca_L+I_stim)/Cm;

if (abs(Y(15)+14.1) < 1.0e-10)
    alpha_xr = 0.0015;
else
    alpha_xr = 0.0003*(Y(15)+14.1)/(1.0-exp((Y(15)+14.1)/-5.0));
end;

if (abs(Y(15)-3.3328) < 1.0e-10)
    beta_xr = 3.7836118e-4;
else
    beta_xr = 0.000073898*(Y(15)-3.3328)/(exp((Y(15)-3.3328)/5.1237)-1.0);
end;

tau_xr = (alpha_xr+beta_xr)^-1.0;
xr_infinity = (1.0+exp((Y(15)+14.1)/-6.5))^-1.0;
dY(16) = (xr_infinity-Y(16))/tau_xr;

if (abs(Y(15)-19.9) < 1.0e-10)
    alpha_xs = 0.00068;
else
    alpha_xs = 0.00004*(Y(15)-19.9)/(1.0-exp((Y(15)-19.9)/-17.0));
end;

if (abs(Y(15)-19.9) < 1.0e-10)
    beta_xs = 0.000315;
else
    beta_xs = 0.000035*(Y(15)-19.9)/(exp((Y(15)-19.9)/9.0)-1.0);
end;

tau_xs = 0.5*(alpha_xs+beta_xs)^-1.0;
xs_infinity = (1.0+exp((Y(15)-19.9)/-12.7))^-0.5;
dY(17) = (xs_infinity-Y(17))/tau_xs;
alpha_oa = 0.65*(exp((Y(15)-(-10.0))/-8.5)+exp((Y(15)-(-10.0)-40.0)/-59.0))^-1.0;
beta_oa = 0.65*(2.5+exp((Y(15)-(-10.0)+72.0)/17.0))^-1.0;
tau_oa = (alpha_oa+beta_oa)^-1.0/K_Q10;
oa_infinity = (1.0+exp((Y(15)-(-10.0)+10.47)/-17.54))^-1.0;
dY(18) = (oa_infinity-Y(18))/tau_oa;
alpha_oi = (18.53+1.0*exp((Y(15)-(-10.0)+103.7)/10.95))^-1.0;
beta_oi = (35.56+1.0*exp((Y(15)-(-10.0)-8.74)/-7.44))^-1.0;
tau_oi = (alpha_oi+beta_oi)^-1.0/K_Q10;
oi_infinity = (1.0+exp((Y(15)-(-10.0)+33.1)/5.3))^-1.0;
dY(19) = (oi_infinity-Y(19))/tau_oi;
alpha_ua = 0.65*(exp((Y(15)-(-10.0))/-8.5)+exp((Y(15)-(-10.0)-40.0)/-59.0))^-1.0;
beta_ua = 0.65*(2.5+exp((Y(15)-(-10.0)+72.0)/17.0))^-1.0;
tau_ua = (alpha_ua+beta_ua)^-1.0/K_Q10;
ua_infinity = (1.0+exp((Y(15)-(-10.0)+20.3)/-9.6))^-1.0;
dY(20) = (ua_infinity-Y(20))/tau_ua;
alpha_ui = (21.0+1.0*exp((Y(15)-(-10.0)-195.0)/-28.0))^-1.0;
beta_ui = 1.0/exp((Y(15)-(-10.0)-168.0)/-16.0);
tau_ui = (alpha_ui+beta_ui)^-1.0/K_Q10;
ui_infinity = (1.0+exp((Y(15)-(-10.0)-109.45)/27.48))^-1.0;
dY(21) = (ui_infinity-Y(21))/tau_ui;
return