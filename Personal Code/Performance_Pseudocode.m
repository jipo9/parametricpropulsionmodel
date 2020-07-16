clear
clc
close all



% Independent variables: M0, P0, T0, Tt4, Tt7, P0/P9
% Constants: beta, pi_d = f(M0), eta or pi of all components, epsilons, 
% ?? :M of 4 & 4.5 = 1??, A of 6, 6A, 16, 8dry
% Dependent variables: pi and tau of all components, f, mdot, alpha, M 6,16,6A,8,9
% 
% Assumption changes:
% - A lot of misc stuff
% - They assume calorically perfect???
%% Appendix I

clc
clear
close all
state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
state(2:18,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9'};
component = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
component(2:17,1) = {'Ram Recovery';'Inlet(Ideal,Actual)';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
design = {'Parameter','Value'};
design(2:8,1) = {'alpha';'beta';'epsilon1';'epsilon2';'PtoL';'PtoH';'h_PR'};

%% Preliminary computations

% run on design analysis!!!

% just use ??? [state, component,v0] = ambient(state,component,alt,M0)
% Fair(f=0,T0)
v0 = M0*a0;
ho0 = V0^2 / 2;
% Fair(f=0,ho0)
tau_r = h0/ho0;
pi_r = Pro0/Pr0;

% use function [state,component] = inlet(state,component,M0) to find pi_d
pi_ABdry = 1 - (1-pi_ABR)/2; %for us just assume that this is 1???

%% Set initial values
% pi_f = component{4,2};
% pi_cL = component{5,2};
% pi_cH = component{6,2};
% pi_tH = component{9,2};
% pi_tL = component{12,2};


% tau_f = component{4,4};
% tau_cL = component{5,4};
% tau_cH = component{6,4};
tau_m1 = component{8,4};
tau_tH = component{9,4};
tau_m2 = component{11,4};
tau_tL = component{12,4};

mdot4 = state{9,5};
mdot45 = state{12,5};

f = state{9,4};

M4 = 1;
M45 = 1;
%M6A
%M8


%Fair(f,Tt4)
ho45 = ho4*tau_m1*tau_tH*tau_m2;
f45 = f*mdot4/mdot45;
%Fair(f45,ho45)
ho5 = ho45*tau_tL;
%Fair(f45,ho5)
ho6A = ho5*tau_M;
f6A = state{15,4};
%Fair (f6A,ho6A)

%% "" Loop "" 1 (will be in its own function)

tau_cL = component{5,4};
tau_cH = component{6,4};
ho0 = state{3,8};
To4 = state{8,3};
alpha = design{2,2};

mdot0 = state{2,5};
mdot4 = state{9,5}; %make sure to change these w/ new f
mdot45 = state{12,5}; %make sure to change these w/ new f

f = state{9,4};



ho3 = ho0*tau_cL*tau_cH;
%Fair (0,ho3)
alpha_prime = alpha* (mdot0/mdot45);
%Fair (0,To4)
f45 = f*mdot4/mdot45;

%TurbC and turb