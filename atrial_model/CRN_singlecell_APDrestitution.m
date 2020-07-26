% Courtemanche-Ramirez-Nattel atrial cell model
% Calculate APD restitution curve (S1S2 protocol)
% Requirement: CRN_singlecell_model.m


clear; clc;

%%
global g_Na g_B_Na g_Ca_L g_B_Ca g_to g_K1 g_Kr g_Ks
global i_NaK_max I_NaCa_max K_mCa K_mNa gamma K_sat Ca_up_max I_up_max K_rel
global K_up i_CaP_max GKur

g_Na = 7.8;   % nanoS_per_picoF (in fast_sodium_current)
g_K1 = 0.09;   % nanoS_per_picoF (in time_independent_potassium_current)
g_to = 0.1652;   % nanoS_per_picoF (in transient_outward_K_current)
g_Kr = 0.029411765;   % nanoS_per_picoF (in rapid_delayed_rectifier_K_current)
g_Ca_L = 0.123759;   % nanoS_per_picoF (in L_type_Ca_channel)
GKur = 1.0;

Ca_up_max = 15.0;   % millimolar (in Ca_leak_current_by_the_NSR)
K_rel = 30.0;   % per_millisecond (in Ca_release_current_from_JSR)
K_up = 0.00092;   % millimolar (in Ca_uptake_current_by_the_NSR)
K_mCa = 1.38;   % millimolar (in Na_Ca_exchanger_current)
K_mNa = 87.5;   % millimolar (in Na_Ca_exchanger_current)
K_sat = 0.1;   % dimensionless (in Na_Ca_exchanger_current)
gamma = 0.35;   % dimensionless (in Na_Ca_exchanger_current)
g_B_Ca = 0.001131;   % nanoS_per_picoF (in background_currents)
g_B_Na = 0.0006744375;   % nanoS_per_picoF (in background_currents)
I_NaCa_max = 1600.0;   % picoA_per_picoF (in Na_Ca_exchanger_current)
I_up_max = 0.005;   % millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)
i_CaP_max = 0.275;   % picoA_per_picoF (in sarcolemmal_calcium_pump_current)
g_Ks = 0.12941176;   % nanoS_per_picoF (in slow_delayed_rectifier_K_current)
i_NaK_max = 0.59933874;   % picoA_per_picoF (in sodium_potassium_pump)

%%
I_stim = -2200;
stimDur = 2;
BCL = 300
stimN = 30;

%% S1
y2 = [2.35e-112, 1.0, 0.9992, 1.367e-4, 7.755e-1, 9.996e-1, 9.649e-1, 9.775e-1, 2.908e-3, 1.013e-4, 1.488, 1.488, 1.39e2, 1.117e1, -81.18, 3.296e-5, 1.869e-2, 3.043e-2, 9.992e-1, 4.966e-3, 9.986e-1];
for i = 1:stimN
    [t1, y1] = ode15s(@(t,y) CRN_singlecell_model(t,y,I_stim), [0 stimDur], y2(end,:));
    [t2, y2] = ode15s(@(t,y) CRN_singlecell_model(t,y,0), [stimDur BCL], y1(end,:));
    t = [t1; t2(2:end,:)]; y = [y1; y2(2:end,:)];
end

T = min(t):0.01:max(t); V = interp1(t, y(:,15), T);
dV = [diff(V) 0];
depolT = T(find(dV == max(dV), 1));
repolV = V(1) + (max(V) - V(1))*0.1;
repolT = T(find(dV<0 & V<=repolV, 1));
APD_BCL = repolT - depolT

[~, y1] = ode15s(@(t,y) CRN_singlecell_model(t,y,I_stim), [0 stimDur], y2(end,:));
[~, y2] = ode15s(@(t,y) CRN_singlecell_model(t,y,0), [stimDur APD_BCL], y1(end,:));
y0 = y2(end,:);

%% S2
DI = []; APD = [];
for DI_ = [100:-10:51 50:-0.5:1]
    [~, y2] = ode15s(@(t,y) CRN_singlecell_model(t,y,0), [0 DI_], y0);
    [t1, y1] = ode15s(@(t,y) CRN_singlecell_model(t,y,I_stim), [0 stimDur], y2(end,:));
    [t2, y2] = ode15s(@(t,y) CRN_singlecell_model(t,y,0), [stimDur 300], y1(end,:));
    t = [t1; t2(2:end,:)]; y = [y1; y2(2:end,:)];
    
    T = min(t):0.01:max(t); V = interp1(t, y(:,15), T);
    dV = [diff(V) 0];
    depolT = T(find(dV == max(dV), 1));
    repolV = V(1) + (max(V) - V(1))*0.1;
    repolT = T(find(dV<0 & V<=repolV, 1));
    APD_ = repolT - depolT;
    
    if max(V) < 0
        break;
    end
    
    DI = [DI; DI_]; APD = [APD; APD_];
end

figure(1); plot(DI, APD, 'k'); xlabel('DI (ms)'); ylabel('APD_{90} (ms)');

