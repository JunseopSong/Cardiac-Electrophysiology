% Courtemanche-Ramirez-Nattel atrial cell model

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
stimCL = [500 400]; stimNum = [10 10];
for i = 1:length(stimCL)
    if i ~= 1
        stimTime = [stimTime stimTime(end)+(1:stimNum(i))*stimCL(i)];
    else
        stimTime = (1:stimNum(i))*stimCL(i);
    end
end

stimN = length(stimTime);
I_stim = -2200;
stimDur = 2;

%%
y = [2.35e-112, 1.0, 0.9992, 1.367e-4, 7.755e-1, 9.996e-1, 9.649e-1, 9.775e-1, 2.908e-3, 1.013e-4, 1.488, 1.488, 1.39e2, 1.117e1, -81.18, 3.296e-5, 1.869e-2, 3.043e-2, 9.992e-1, 4.966e-3, 9.986e-1];
[~, y] = ode15s(@(t,y) CRN_singlecell_model(t,y,0), [0 10000], y);
y = y(end,:);

if stimTime(1) ~= 0
    [t, y] = ode15s(@(t,y) CRN_singlecell_model(t,y,0), [0 stimTime(1)], y);
else
    t = [0]; y = [y];
end

for i = 1:stimN
    if i == stimN
        nextTime = stimTime(end) + 2000;
    else
        nextTime = stimTime(i+1);
    end
    
    [tt yy] = ode15s(@(t,y) CRN_singlecell_model(t,y,I_stim), [stimTime(i) stimTime(i)+stimDur], y(end,:));
    t = [t; tt(2:end,:)]; y = [y; yy(2:end,:)];
    
    [tt yy] = ode15s(@(t,y) CRN_singlecell_model(t,y,0), [stimTime(i)+stimDur nextTime], y(end,:));
    t = [t; tt(2:end,:)]; y = [y; yy(2:end,:)];
end

T = min(t):0.1:max(t);
for i = 1:size(y,2)
    Y(:,i) = interp1(t, y(:,i), T);
end
clear t y tt yy;

%%
figure(1); plot(T, Y(:,15), 'k');
