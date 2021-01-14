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
design = {'Parameter','Value'};
design(2:8,1) = {'alpha';'beta';'epsilon1';'epsilon2';'PtoL';'PtoH';'h_PR'};

%% Peliminary computations

% run on design analysis!!!

% just use ??? [state, component,v0] = ambient(state,component,alt,M0)
% Fair(f=0,T0)
v0 = M0*a0;
ho0 = V0^2 / 2;
% Fair(f=0,ho0)
tau_r = h0/ho0;
pi_r = Po0/P0;

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
MFP6 = MFP8 * (pi_M?*pi_ABdry)/(1+alpha_prime) * (A8/A6) * sqrt(To6/To6A);
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
    mdot0 = mdot0_new
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


%% Section 5 overall performance




























%go through and change everything to either normal values (f) or iteration? (f_R)