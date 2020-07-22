clear
clc
close all

state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
state(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
component = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
component(2:17,1) = {'Ram Recovery';'Inlet Actual';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
design = {'Parameter','Value'};
design(2:8,1) = {'alpha';'beta';'epsilon1';'epsilon2';'PtoL';'PtoH';'h_PR'};
inputs = {'Parameter','Value'};
inputs(2:8,1) = {'Altitude (m)','Mach Number','F/mdot','Mass Flow Rate (kg/s)','SFC (kg/s/N)','Max Burner Temp (K)','Po9/P9'};

M0 = 1.451;
alt = 36000/3.281; %altitude [m from feet]
To4 = 3200*.555556;
To7 = 3600*.555556;
pi_n = .97;

component(2:16,2) = {3.4211; .9354; 3.9000; 3.9000; 5.1282; []; [];    .4231; []; [];    .4831; []; .9779; []; pi_n};
component(2:16,4) = {1.4211; [];    1.5479; 1.5479; 1.6803; []; .9673; .8381; []; .9731; .8598; []; .8404; []; []};

mdot = 187.45*0.45359237;
alpha = .449;
beta = .01;
ep1 = .05;
ep2 = .05;
f = .03070;

state{2,5} = mdot;
state{9,4} = f;
state{9,3} = To4;
state{16,3} = To7;

inputs{2,2} = alt;
inputs{3,2} = M0;

design{2,2} = alpha;
design{3,2} = beta;
design{4,2} = ep1;
design{5,2} = ep2;


[state,design] = derived_parameters_performance(state,inputs,design,component);

%% State 0
[T0, ~, ~, ~] = atmosisa(alt);
state(2,3) = {T0};
[state] = unFAIR3(state,2);

%% State 0o
tau_r = component{2,4};
T_o0 = T0*tau_r;
state(3,3) = {T_o0};
[state] = unFAIR3(state,3);
P_o0 = state{3,2};

%% State 2o
pi_d = component{3,2};
P_o2 = P_o0*pi_d;
state(4,2) = {P_o2};
[state] = unFAIR3(state,4);
T_o2 = state{4,3};

%% State 1.3o
tau_f = component{4,4};
T_o13 = T_o2*tau_f;
state(5,3) = {T_o13};
[state] = unFAIR3(state,5);

%% State 2.5o
tau_cL = component{5,4};
T_o25 = T_o2*tau_cL;
state{6,3} = T_o25;
[state] = unFAIR3(state,6);

%% State 3o
tau_cH = component{6,4};
T_o3 = T_o25*tau_cH;
state{7,3} = T_o3;
[state] = unFAIR3(state,7);

%% State 3.1o
T_o31 = T_o3;
state{8,3} = T_o31;
[state] = unFAIR3(state,8);

%% State 4o
state{9,3} = To4;
[state] = unFAIR3(state,9);

%% State 4.1o
tau_m1 = component{8,4};
T_o41 = To4*tau_m1;
state{10,3} = T_o41;
[state] = unFAIR3(state,10);

%% State 4.4o
tau_tH = component{9,4};
T_o44 = T_o41*tau_tH;
state{11,3} = T_o44;
[state] = unFAIR3(state,11);

%% State 4.5o
tau_m2 = component{11,4};
T_o45 = T_o44*tau_m2;
state{12,3} = T_o45;
[state] = unFAIR3(state,12);

%% State 5o
tau_tL = component{12,4};
T_o5 = T_o45*tau_tL;
state{13,3} = T_o5;
[state] = unFAIR3(state,13);

%% State 6o (stop here)




function [state,design] = derived_parameters_performance(state,inputs,design,component)
%Derived parameters for performance model, w/ changing values

mdot = state{2,5};
alpha = design{2,2};
beta = design{3,2};
ep1 = design{4,2};
ep2 = design{5,2};
f = state{9,4};
f0 = 0;

% MASS FLOW AND FUEL TO AIR CALCULATIONS

mdot25 = mdot/(1+alpha); %after bypass leaves
mdot13 = mdot - mdot25; % bypass mass flow

mdot31 = mdot25*(1-beta - ep1 -ep2); %after bleed and coolant leaves
mdotbeta = mdot25*beta; %bleed air
mdotep1 = mdot25*ep1; %coolant air 1
mdotep2 = mdot25*ep2; %coolant air 2
mdotep = mdotep1+mdotep2; %combined coolant air

mdot_f = mdot31 * f; %mass flow per fuel/air ratio
mdot4 = mdot31 + mdot_f; %mass flow rate post-burner
f4 = mdot_f / mdot31; %fuel/air ratio post-burner

mdot41 = mdot4 + mdotep1; %mass flow rate after addtion of cooling air 1
f41 = f4*mdot4 / mdot41; %fuel/air ratio after addtion of cooling air 1

mdot45 = mdot41 + mdotep2; %mass flow rate after addtion of cooling air 2
f45 = f41*mdot41 / mdot45; %fuel/air ratio after addtion of cooling air 2

mdot6A = mdot45 + mdot13; %mass flow rate after addtion of bypass air
f6A = f45*mdot45 / mdot6A; %fuel/air ratio after addtion of bypass air

% Store all values
state(2:8,4) = {f0};
state(2:4,5) = {mdot};
state(5,5) = {mdot13};
state(6:7,5) = {mdot25};
state(8,5) = {mdot31};
state(9,4) = {f4};
state(9,5) = {mdot4};
state(10:11,4) = {f41};
state(10:11,5) = {mdot41};
state(12:13,4) = {f45};
state(12:13,5) = {mdot45};
state(14:18,4) = {f6A};
state(14:18,5) = {mdot6A};
state(19,5) = {mdotbeta};
state(20,5) = {mdotep};
state(21,5) = {mdotep1};
state(22,5) = {mdotep2};

end

