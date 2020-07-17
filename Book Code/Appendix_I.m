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
state = {'Station','Relative Pessure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
state(2:18,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9'};
component = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
component(2:17,1) = {'Ram Recovery';'Inlet(Ideal,Actual)';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
component_ref = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
component_ref(2:17,1) = {'Ram Recovery';'Inlet(Ideal,Actual)';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
design = {'Parameter','Value'};
design(2:8,1) = {'alpha';'beta';'epsilon1';'epsilon2';'PtoL';'PtoH';'h_PR'};
inputs = {'Parameter','Value'};
inputs(2:8,1) = {'Altitude (m)','Mach Number','F/mdot','Mass Flow Rate (kg/s)','SFC (kg/s/N)','Max Burner Temp (K)','Po9/P9'};


%% INPUTS

pi_dmax = .96;
pi_M_max = .97;
pi_AB = .95;
pi_b = .95;
pi_n = .97;

eta_cL = .8674;
eta_mL = .995;
eta_f = .8674;
eta_AB = .99;
eta_cH = .8751;
eta_mH = .995;
eta_b = .999;
eta_tH = .8995;
eta_PL = 1;
eta_PH = 1;
eta_tL = .9074;



PtoL = 0;
PtoH = 281900; % W
h_PR = 18400*2326; %J/kg

beta = .01;
ep1 = .05;
ep2 = .05;


% gamma_c = 1.4;
% gamma_t = 1.3;
% gamma_AB = 1.3;
% cp_c = .24* 4186.8; %J/kg K
% cp_t = .295* 4186.8; %J/kg K
% cp_AB = .295* 4186.8; %J/kg K

M0 = 1.8;
alt = 40000/3.281; %altitude [m from feet]

%Control limits
T_t4 = 3200*.5556; %max burner temperature [R to K]
pi_c = 20;



design(3,2) = {beta};
design(7,2) = {PtoH};
design(6,2) = {PtoL};
design(8,2) = {h_PR};

state(9,3) = {T_t4};

component(3,2) = {pi_dmax}; %store values in component
component(14,2) = {pi_M_max};
component(7,2) = {pi_b};
component(15,2) = {pi_AB};
component(16,2) = {pi_n};

component(4,5) = {eta_f};
component(5,5) = {eta_cL};
component(6,5) = {eta_cH};
component(8,5) = {eta_mH};
component(9,5) = {eta_tH};
component(10,5) = {eta_PH}; 
component(11,5) = {eta_mL};
component(12,5) = {eta_tL};
component(13,5) = {eta_PL}; 

inputs(2,2) = {alt};
inputs(3,2) = {M0};


%% Peliminary computations

% run on design analysis!!!
state(2:18,4) = {0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0};
[state, component,v0] = ambient(state,component,inputs);
[state,component] = inlet(state,component,inputs);%Pressures are calculated a bit weird here...)
pi_ABR = component{15,2}; %??? is this right?/?
pi_ABdry = 1 - (1-pi_ABR)/2;

%% Set initial values
pi_fR = 3.9;
pi_cLR = 3.9;
pi_cHR = 5.1282;
pi_tHR = .4231;
pi_tLR = .4831;

tau_fR = 1.5479;
tau_cLR = 1.5479;
tau_cHR = 1.6803;
tau_m1R = .9673;
tau_tHR = .8381;
tau_m2R = .9731;
tau_tLR = .8598;

component_ref{4,2} = pi_fR;
component_ref{5,2} = pi_cLR;
component_ref{6,2} = pi_cHR;
component_ref{9,2} = pi_tHR;
component_ref{12,2} = pi_tLR;

component_ref{4,4} = tau_fR;
component_ref{5,4} = tau_cLR;
component_ref{6,4} = tau_cHR;
component_ref{8,4} = tau_m1R;
component_ref{9,4} = tau_tHR;
component_ref{11,4} = tau_m2R;
component_ref{12,4} = tau_tLR;

fR = .03070;
tau_MR = .8404;
alphaR = .449;
M6AR = .4188;
M8R = 1; %???



pi_f = pi_fR;
pi_cL = pi_cLR;
pi_cH = pi_cHR;
pi_tH = pi_tHR;
pi_tL = pi_tLR;

tau_f = tau_fR;
tau_cL = tau_cLR;
tau_cH = tau_cHR;
tau_m1 = tau_m1R;
tau_tH = tau_tHR;
tau_m2 = tau_m2R;
tau_tL = tau_tLR;

component{4,2} = pi_f;
component{5,2} = pi_cL;
component{6,2} = pi_cH;
component{9,2} = pi_tH;
component{12,2} = pi_tL;

component{4,4} = tau_f;
component{5,4} = tau_cL;
component{6,4} = tau_cH;
component{8,4} = tau_m1;
component{9,4} = tau_tH;
component{11,4} = tau_m2;
component{12,4} = tau_tL;

% mdot4 = state{9,5};
% mdot45 = state{12,5};

f = fR;
state(9,4) = {f};

M4 = 1;
M45 = 1;
M6A = M6AR;
M8 = M8R;



%% Initial calcs

[state] = unFAIR3(state,9);
ho4 = state{9,8};

ho45 = ho4*tau_m1*tau_tH*tau_m2;
state(12,8) = {ho45};
f45 = f*(1-beta-ep1-ep2)/(1-beta);
state(12,4) = {f45};
[state] = unFAIR3(state,12);

ho5 = ho45*tau_tL;
state(13,8) = {ho5};
state(13,4) = {f45};
[state] = unFAIR3(state,13);

ho6A = ho5*tau_MR; %Is this the right eqn?
f6A = f45*(1-beta)/(1+alphaR-beta); %Is this the right eqn?
state(15,8) = {ho6A};
state(15,4) = {f6A};
[state] = unFAIR3(state,15);


%% Current checkpoint
    % Rewrite RGCompr, Massfp, turb, turbc
    % then proceed
%% "" Loop "" 1 (will be in its own function) - Calculates states of comPessors

tau_cL = component{5,4};
tau_cH = component{6,4};
ho0 = state{3,8};
To4 = state{8,3};
alpha = design{2,2};

mdot0 = state{2,5};
mdot4 = state{9,5}; %make sure to change these w/ new f
mdot45 = state{12,5}; %make sure to change these w/ new f
mdotep1 = state{21,5};
mdotep2 = state{22,5};

f = state{9,4};



ho3 = ho0*tau_cL*tau_cH;
%Fair (0,ho3)
alpha_Pime = alpha* (mdot0/mdot45);
%Fair (0,To4)
f45 = f*mdot4/mdot45;

%TurbC and turb???
tau_lambda = ho4/h0;
tau_f = 1 ...
    + ( (1-tau_tL)*eta_mL* ((mdot4*tau_lambda*tau_tH/tau_r  +  (mdotep1*tau_tH + mdotep2)*tau_cL*tau_cH)/mdot0)   -   (1+alpha)*PtoL/(tau_r*eta_mPL*mdot0*h0)) ...
    / ((tau_cL - 1)/(tau_f - 1) + alpha); %where did this come from???
tau_cL = 1 + (tau_f -1)*(tau_cL - 1)/(tau_f - 1);
tau_cH = 1 ...
    + (1-tau_tH)*eta_mH* ((mdot4*tau_lambda/(tau_r*tau_cL) + mdotep1*tau_cH)/mdot0)...
    - (1+alpha)*PtoH/(tau_r*tau_cL*eta_mPH*mdot0*h0);
ho2 = ho0;
Po2 = Po0;
ho13 = ho2*tau_f;
ho25 = ho2*tau_cL;
ho3 = ho25*tau_cH;
ho13i = ho2*(1+eta_f*(tau_f-1));
ho25i = ho2*(1+eta_cL*(tau_cL-1));
ho3i = ho25*(1+eta_cH*(tau_cH-1));
%Fair(f=0,ho25)
%Fair(f=0,ho13)
%Fair(f=0,ho3)
%Fair(f=0,ho13i)
%Fair(f=0,ho25i)
%Fair(f=0,ho3i)

pi_f = Po13i / Po2;
pi_cL = Po25i / Po2;
pi_cH = Po3i / Po25;
pi_c = pi_cH*pi_cL;
tau_c = tau_cL*tau_cH;

%% %% "" Loop "" 2 (will be in its own function) - Calculates fuel to air ratio

f_temp = f;
%Fair(f,To4)
f = (ho4 - ho3) / (h_PR*eta_b - ho4);
% if norm(f-f_temp) > .0001 goto 2

%% Section 2.5 Mixer pt 1

mdot31 = state{8,5};
mdotep1 = state{21,5};
mdotep2 = state{22,5};

ho5 = ho5;
To6 = To5;
ho16 = ho13;
To16 = To13;
Po6 = P0*pi_r*pi_d*pi_cL*pi_cH*pi_b*pi_tH*pi_tL;
f45 = f*mdot31 / (mdot31+mdotep1+mdotep2);
%RGCOMP(1,To6,M6)
P6 = Po6/Po6_P6;
T6 = To6/To6_T6;
Po16 = P0*pi_r*pi_d*pi_f;
Po16_P16 = Po16/P6;
%RGCOMP(1,To16,M=1)???

% if Po16_P16 > PtP then M6 =- .01 and go to 1
% if Po16_P16 < 1 then M6 =- .01 and go to 1

%% Section 2.75 Mixer pt 2
%RGCOMP(1,To16,M16)
T16 = To16/To16_T16;
alpha_prime_new = Po16*A16*MFP16/sqrt(To16) ...
    / (Po6*A6*MFP6/sqrt(To6));
alpha_prime_error = norm((alpha_prime_new-alpha_prime)/alpha_prime);
alpha = alpha_prime_new *(mdot5/mdot0);

% if alpha_prime_error > .001
    alpha_prime = alpha_prime + alpha_prime / (alpha_prime_new -alpha_prime); %newtons iteration method??
% Goto 1

%% Section 2.875 More mixer
%Fair(f45,T6)
%Fair(f = 0,T16)
ho6A = (ho6 +alpha_prime*ho16) / (1+alpha_prime);
tau_M = ho6A/ho6;
f6A = f45*(1- beta)/(1 + alpha - betta);
% fair(f6a,ho6A)
constant = 1/(1+alpha_prime)...
    *(sqrt(R6*T6/gamma6)*(1+gamma6*M6^2)/M6...
    + alpha_prime * sqrt(R16*T16/gamma16)*(1+gamma16*M16^2)/M16);
%% Section 3 Even more mixer
% RGcompr (To6a,M6a,f6A)
T6A = To6A / To6A_T6A;
% Fair(f6A,T6A)
M6A_new = sqrt(R6A*T6A/gamma6A)*(1+gamma6A*M6A^2)/constant;
M6A_error = norm(M6A_new-M6A);

%if M6_error > .001 then M6A = M6A_new and go to 3
%% Section 3.5 Even more more mixer

pi_M_ideal = sqrt(To6A / To6) * MFP6/MFP6A * (1+alpha_prime)/(1+A16/A6);
pi_M = pi_M_max*pi_M_ideal;
Po9_P0 = pi_r*pi_d*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_M*pi_ABdry*pi_n;
% RGcompr (To6a,M9,f6A)
if M9 > 1
    M8 = 1;
else
    M8 = M9;
end
% RGcompr (To6a,M8,f6A)
MFP6 = MFP8 * (pi_M*pi_ABdry)/(1+alpha_prime) * (A8/A6) * sqrt(To6/To6A);
% RGcompr (To6,M6new,f45);
M6_error = norm(M6 - M6_new);
if M6_error > .0005
    if M6>M6_new
        M6 = M6 - .0001;
    else
        M6 = M6 + .0002;
    end
    % go to 1
end

%% Section 3.75 Overall mass flow
% RGcompr (To4,M4,f);
mdot0_new = mdot0 ...
    * (1+fR)*sqrt(To4R)/ (P0R*(1+alphaR)*pi_rR*pi_dR*pi_cLR*pi_cHR*MFP4R) ...
    / (1+f)*sqrt(To4)/ (P0*(1+alpha)*pi_R*pi_d*pi_cL*pi_cH*MFP4);
mdot0_error = norm((mdot0_new - mdot0 )/ mdot0R);
if mdot0_error > .001
    mdot0 = mdot0_new;
    %goto 1
end
f7 = f6A;
%% Section 4 Afterburner fuel
%Fair(f7,Tt7)
fAB =  (ho7 - ho6A) / (eta_AB*h_PR - ho6A);
f7_new = f6A + fAB;
f7_error = norm(f7_new-f7);
if f7_error > .00001
    f7 =f7_new;
    %goto 4
end
%% Section 4.5 Afterburner performance
percent_AB = 100*(To7 - To6A)/(To7R - To6AR);
pi_AB = pi_ABdry + .01*percent_AB*(pi_ABR - pi_ABdry)
%% Section 4.75 Nozzle conditions
Po9_P0 = pi_r*pi_d*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_M*pi_AB*pi_n;
Po9_P9 = Po9_P0/P0_P9;
To9 = To7;
%RGCOMPR(3,To9,M9,f7)

mdot9 = mdot0*(1+f7)*(1-betta/(1+alpha));
Po9 = P0*Po9_P0;
A9 = mdot9*sqrt(To9)/(Po9*MFP9);
T9 = To9/To9_T9;
%Fair(f7,T9)
v9 = M9*a9;
f0 = (f*(1-betta-ep1-ep2) + fAB*(1+alpha-beta))  /  (1+alpha);
%% Section 5 overall performance





F_mdot = (1+f_0-(beta/(1+alpha)))*v9     -   v0  +   (1+f_0-(beta/(1+alpha)))*R9*T9*(1-Pr0/Pr9)/(R0*v9*gamma0);
S = f_0 / F_mdot;



% Remainder of stuff surpressed for now
% 
% 
% 
% 
% 
% eta_TH = ((v0^2/2*((1+f_0-(beta/(1+alpha)))*(v9/v0)^2 - 1) + (PtoL + PtoH)/mdot0))/...
%     (f_0*h_PR);
% %eta_P = 2/(1+v9/v0); Simplified case, neglected for now
% eta_P = (2*F_mdot/v0)/...
%     ((1+f_0-beta/(1+alpha))*(v9/v0)^2 - 1);
% eta_o = eta_TH*eta_P;
% % eta_o = v0/(h_PR*S); Alternate eqn























%go through and change everything to either normal values (f) or iteration? (f_R)















function [state, component,v0] = ambient(state,component,inputs)
alt = inputs{2,2};
M0 = inputs{3,2};
[T0, ~, ~, ~] = atmosisa(alt); %obtain standard atmospheric conditions
state(2,3) = {T0};
[state] = unFAIR3(state,2);
[~,~,T0,~,~,cp0,gamma0,~,~] = state{2,:};
R0 = cp0 - cp0/gamma0;
a0 = sqrt(R0*gamma0*T0); %[m/s]
v0 = M0*a0; %[m/s]

T_o0 = T0*(1+((M0^2)*((gamma0-1)/2))); %find total temperature using isentropic
state(3,3) = {T_o0};
[state] = unFAIR3(state,3);


P0 = state{2,2};
Po0 = state{3,2};
pi_r = Po0/P0;
component{2,2} = pi_r;

h0 = state{2,8};
ho0 = state{3,8};
tau_r = ho0/h0;
component{2,4} = tau_r;
end
function [state,component] = inlet(state,component,inputs)
M0 = inputs{3,2};
pi_dmax = component{3,2};
Po0 = state{3,2};

if M0<1
pid = pi_dmax; 
elseif M0>1 && M0<5
pid = pi_dmax * (1-.075*((M0-1)^1.35));
else 
pid = pi_dmax * (800/(M0^4 + 935));
end

Po2 = pid*Po0;

component{3,2} = [pi_dmax,pid];
state(4,2) =  {Po2};
[state] = unFAIR3(state,4);

ho0 = state{3,8};
ho2 = state{4,8};
tau_d = ho2/ho0;
component{3,4} = tau_d;
end