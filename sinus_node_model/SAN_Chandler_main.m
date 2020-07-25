clear all; clc;

global CMDN_max CSQN_max Km_CMDN Km_CSQN Km_TRPN TRPN_max Ca_up_max K_rel I_up_max K_up g_Ca_L;
global I_NaCa_max K_mCa K_mNa K_sat gamma g_B_Ca g_B_K g_B_Na g_Na V_cell Cm F R T g_Kr i_CaP_max;
global g_Ks Km_K_o Km_Na_i i_NaK_max Ca_o K_o Na_o g_K1 tau_tr K_Q10 g_to G_Kur G_f E_f G_CaT E_CaT;

CMDN_max = 0.05;   % millimolar (in Ca_buffers)
CSQN_max = 10.0;   % millimolar (in Ca_buffers)
Km_CMDN = 0.00238;   % millimolar (in Ca_buffers)
Km_CSQN = 0.8;   % millimolar (in Ca_buffers)
Km_TRPN = 0.0005;   % millimolar (in Ca_buffers)
TRPN_max = 0.07;   % millimolar (in Ca_buffers)
Ca_up_max = 15.0;   % millimolar (in Ca_leak_current_by_the_NSR)
K_rel = 30.0*0.31;   % per_millisecond (in Ca_release_current_from_JSR)
I_up_max = 0.005*0.30;   % millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)
K_up = 0.00092;   % millimolar (in Ca_uptake_current_by_the_NSR)
g_Ca_L = 0.12375*0.68;   % nanoS_per_picoF (in L_type_Ca_channel)
I_NaCa_max = 1600.0*0.5;   % picoA_per_picoF (in Na_Ca_exchanger_current)
K_mCa = 1.38;   % millimolar (in Na_Ca_exchanger_current)
K_mNa = 87.5;   % millimolar (in Na_Ca_exchanger_current)
K_sat = 0.1;   % dimensionless (in Na_Ca_exchanger_current)
gamma = 0.35;   % dimensionless (in Na_Ca_exchanger_current)
g_B_Ca = 0.001131;   % nanoS_per_picoF (in background_currents)
g_B_K = 0.0;   % nanoS_per_picoF (in background_currents)
g_B_Na = 0.0006744375;   % nanoS_per_picoF (in background_currents)
g_Na = 7.8*0.06;   % nanoS_per_picoF (in fast_sodium_current)
V_cell = 10050.0;   % micrometre_3 (in intracellular_ion_concentrations)
Cm = 50.0;   % picoF (in membrane)

F = 96.4867;   % coulomb_per_millimole (in membrane)
R = 8.3143;   % joule_per_mole_kelvin (in membrane)
T = 310.0;   % kelvin (in membrane)

g_Kr = 0.029411765*0.45;   % nanoS_per_picoF (in rapid_delayed_rectifier_K_current)
i_CaP_max = 0.275;   % picoA_per_picoF (in sarcolemmal_calcium_pump_current)
g_Ks = 0.12941176*0.69;   % nanoS_per_picoF (in slow_delayed_rectifier_K_current)
Km_K_o = 1.5;   % millimolar (in sodium_potassium_pump)
Km_Na_i = 10.0;   % millimolar (in sodium_potassium_pump)
i_NaK_max = 0.59933874;   % picoA_per_picoF (in sodium_potassium_pump)
Ca_o = 1.8;   % millimolar (in standard_ionic_concentrations)
K_o = 5.4;   % millimolar (in standard_ionic_concentrations)
Na_o = 140.0;   % millimolar (in standard_ionic_concentrations)
g_K1 = 0.09*0.37;   % nanoS_per_picoF (in time_independent_potassium_current)
tau_tr = 180.0;   % millisecond (in transfer_current_from_NSR_to_JSR)
K_Q10 = 3.0;   % dimensionless (in transient_outward_K_current)
g_to = 0.1652*0.40;   % nanoS_per_picoF (in transient_outward_K_current)
G_Kur = 0.26; % dimensionless

G_f = 3.76; % nS (If conductance)
E_f = -22.0; % mV (If reversal potential)
G_CaT = 11.0; % nS (I_CaT conductance)
E_CaT = 45.0; % mV (I_CaT reversal potential)


%%
Y0 = [ 0.0000    1.0000    0.9926    0.0156    0.3899    0.4100    0.0024    0.0017    0.4848    0.0005    4.9738    4.9837  139.7844    9.4330  -43.3145 0.4148    0.1769    0.2238    0.2092    0.2233    0.9814    0.0038    0.0561    0.0335    0.1619];

[time,variables] = ode15s(@SAN_Chandler_model,[0,4000],Y0);

variables(end,:);


%% If
I_f = G_f*variables(:,22).*(variables(:,15)-E_f);

%% I_CaL
I_CaL = Cm*g_Ca_L*variables(:,4).*variables(:,6).*variables(:,5).*(variables(:,15)-65.0);

%% I_CaT
I_CaT = G_CaT*variables(:,23).*variables(:,24).*(variables(:,15)-E_CaT);

%% I_NaCa
I_NaCa = Cm*I_NaCa_max*(exp(gamma*F*variables(:,15)/(R*T)).*variables(:,14).^3.0*Ca_o-exp((gamma-1.0)*F*variables(:,15)/(R*T)).*Na_o^3.0.*variables(:,10))/((K_mNa^3.0+Na_o^3.0)*(K_mCa+Ca_o)*(1.0+K_sat*exp((gamma-1.0)*variables(:,15)*F/(R*T))));

%% Print
figure(1);
hold on
subplot(2,1,1);
plot(time,variables(:,15));
xlabel('Time (ms)');
ylabel('V_{m} (mV)');

subplot(2,1,2);
plot(time, I_f, 'b');
xlabel('Time (ms)');
ylabel('I (pA)');

hold on
plot(time, I_CaL, 'r');
hold on
plot(time, I_CaT, 'g');
hold on
plot(time, I_NaCa, 'k');