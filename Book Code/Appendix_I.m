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







% If everything else is correct and using traditional method for
% calculating turbines, we get:
    % .1% enthalpy error in m1
    % 2.1% enthalpy error in tH
    % .3% enthalpy error in m2
    % 1.3% enthalpy error in tC
%% On Design Analysis


clc
clear
close all
stateR = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
stateR(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
componentR = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
componentR(2:17,1) = {'Ram Recovery';'Inlet(Ideal,Actual)';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
designR = {'Parameter','Value'};
designR(2:8,1) = {'alpha';'beta';'epsilon1';'epsilon2';'PtoL';'PtoH';'h_PR'};
inputsR = {'Parameter','Value'};
inputsR(2:8,1) = {'Altitude (m)','Mach Number','F/mdot','Mass Flow Rate (kg/s)','SFC (kg/s/N)','Max Burner Temp (K)','Po9/P9'};

% gamma_c = 1.4;
% gamma_t = 1.3;
% gamma_AB = 1.3;
% cp_c = .24* 4186.8; %J/kg K
% cp_t = .295* 4186.8; %J/kg K
% cp_AB = .295* 4186.8; %J/kg K

alt = 36000/3.281; %altitude [m from feet]
M0 = 1.451; %freestream mach number
F_mdot = 62.493*9.806655; %thrust/mdot [N/kg/s from lbf/(lbm/s)]
mdot = 200*0.45359237; %freestream mass flow rate [kg/s from lbm/s]
S = 1.0862*((.453592/3600)/4.44822); %specific fuel consuption[kg/s/N from lbm/(lbf/s)]
T_t4 = 3200*.555556; %max burner temperature [R to K]
Po9_P9 = 11.327; %stagnation to static pressure ratio at the nozzle (note: in future analysis calculated by total pressure ratio)

alpha = .4; %bypass ratio (-1?)
beta = .01; %bleed ratio
PtoH = 300.61*10^3; %power takeoff high spool [W]
PtoL = 0.000000001; %power takeoff low spool [W]
h_PR = 18400*2326; %fuel heating value for a (CH2)n propellant [J/kg]

pi_dmax = .96; %diffuser pressure ratio
pif = 3.9; %fan pressure ratio
ef = .89; %fan polytropic efficiency
picL = 3.9; %low pressure compressor pressure ratio
ecL = .89; %low pressure polytropic efficiency
    picH = 5.1282; %high pressure compressor pressure ratio
ecH = .9; %high pressure polytropic efficiency
eta_b = .999; %burner efficiency
pi_b = .95; % burner pressure ratio
etH = .89; %high pressure turbine polytropic efficiency
etL = .9; %low pressure turbine polytropic efficiency
etamH = .995; %high pressure shaft mechanical efficiency
etamPH = 1; %high pressure shaft power takeoff mechancal efficiency
etamL = .995; %low pressure shaft mechanical efficiency
etamPL = 1; %low pressure shaft power takeoff mechancal efficiency
pi_M_max = .97; %mixer pressure ratio
pin = .97; %nozzle pressure ratio


% %Control limits
% pi_c = 20;


% Store all values
designR(2,2) = {alpha}; %store values in design
designR(3,2) = {beta};
designR(7,2) = {PtoH};
designR(6,2) = {PtoL};
designR(8,2) = {h_PR};
stateR(9,3) = {T_t4};
componentR(3,2) = {pi_dmax}; %store values in component
componentR(4,2:3) = {pif,ef};
componentR(5,2:3) = {picL,ecL};
componentR(6,2:3) = {picH,ecH};
componentR(7,2) = {pi_b};
componentR(7,5) = {eta_b};
componentR(9,5) = {etamH};
componentR(10,5) = {etamPH};
componentR(9,3) = {etH};
componentR(12,5) = {etamL};
componentR(13,5) = {etamPL};
componentR(12,3) = {etL};
componentR(14,2) = {pi_M_max};
componentR(16,2) = {pin};
componentR(2:3,3) = {1};
componentR(7:8,3) = {1};
componentR(11,3) = {1};
componentR(14:16,3) = {1};
inputsR(2,2) = {alt};
inputsR(3,2) = {M0};
inputsR(4,2) = {F_mdot};
inputsR(5,2) = {mdot};
inputsR(6,2) = {S};
inputsR(7,2) = {T_t4};
inputsR(8,2) = {Po9_P9};

[stateR,designR] = derived_parameters_parametric(stateR,inputsR,designR,componentR);

prompt = 'Run analysis with combined compressors/turbines? Y/N : ';
str = input(prompt,'s');
if str == ('Y') || str == ('y')
    clc
    [stateR,componentR,performanceR] = component_combined(stateR,componentR,designR,inputsR);
elseif str == ('N') || str == ('n')
    clc
    [stateR,componentR,performanceR] = component_seperate(stateR,componentR,designR,inputsR);
else
    error('Invalid Input. Try again and only type Y or N you dummy')
end

err_T_mdot = performanceR{2,1} /F_mdot; %T/mdot error compared to book
err_s = performanceR{2,2} / S; %S error compared to book
err_efftherm = performanceR{2,3} / .5589; %thermal efficiency error compared to book
err_effprop =performanceR{2,4} / .6162; %propulsive efficiency compared to book

fprintf('%s%.3f%s\n','Thrust                    of this analysis is ',abs(100*(1-err_T_mdot)),'% off book solution.')
fprintf('%s%.3f%s\n','Specific Fuel Consumption of this analysis is ',abs(100*(1-err_s)),'% off book solution.')
fprintf('%s%.3f%s\n','Thermal Efficiency        of this analysis is ',abs(100*(1-err_efftherm)),'% off book solution.')
fprintf('%s%.3f%s\n','Propulsive Efficiency     of this analysis is ',abs(100*(1-err_effprop)),'% off book solution.')




%% Off Design Analysis

state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
state(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
component = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
component(2:17,1) = {'Ram Recovery';'Inlet(Ideal,Actual)';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
design = {'Parameter','Value'};
design(2:8,1) = {'alpha';'beta';'epsilon1';'epsilon2';'PtoL';'PtoH';'h_PR'};
inputs = {'Parameter','Value'};
inputs(2:8,1) = {'Altitude (m)','Mach Number','F/mdot','Mass Flow Rate (kg/s)','SFC (kg/s/N)','Max Burner Temp (K)','Po9/P9'};


%should be automatically input, not hardcoded in
eta_f = .8674;
eta_cL = .8674;
eta_cH = .8751;
eta_tH = .8995;
eta_tL = .9074;
component(:,5) = componentR(:,5);
    % eta_b = .999;
    % eta_PL = 1;
    % eta_PH = 1;
    % eta_AB = .99;
component(4,5) = {eta_f};
component(5,5) = {eta_cL};
component(6,5) = {eta_cH};
component(9,5) = {eta_tH};
component(12,5) = {eta_tL};

pi_b = .95;
pi_AB = .95;
component(3,2) = {componentR{3,2}(1)};
component(7,2) = {pi_b};
component(15,2) = {pi_AB};
component(16,2) = componentR(16,2);

alt = 4200/3.281; %altitude [m from feet]
M0 = 1.9; %freestream mach number
inputs(2,2) = {alt};
inputs(3,2) = {M0};

design = designR;
PtoH = 281.9*10^3; %power takeoff high spool [W]
PtoL = 0; %power takeoff high spool [W]
design(2,2) = {[]};
design(6,2) = {PtoL};
design(7,2) = {PtoH};














% mdot = state{2,5};
% alpha = design{2,2};
% beta = design{3,2};
% ep1 = design{4,2};
% ep2 = design{5,2};
% f = state{9,4};
% f0 = 0;




%% Peliminary computations


[state,design] = derived_parameters_performance(state,inputs,design,component);
[state, component,v0] = ambient(state,component,inputs);
[state,component] = inlet(state,component,inputs);%Pressures are calculated a bit weird here...)



pi_ABR = component{15,2}; %??? is this right?/? Its temporary anyways
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
tau_MR = .8404;


% mdot0R = 200*0.45359237; %kg/s
% alphaR = .4;

componentR{4,2} = pi_fR;
componentR{5,2} = pi_cLR;
componentR{6,2} = pi_cHR;
componentR{9,2} = pi_tHR;
componentR{12,2} = pi_tLR;

componentR{4,4} = tau_fR;
componentR{5,4} = tau_cLR;
componentR{6,4} = tau_cHR;
componentR{8,4} = tau_m1R;
componentR{9,4} = tau_tHR;
componentR{11,4} = tau_m2R;
componentR{12,4} = tau_tLR;
componentR{14,4} = tau_MR;


fR = .03069;
stateR(9,4) = {fR};

% alphaR = .449;
M6AR = .4188;
M8R = 1; %???

component(4,2) = componentR(4,2);
component(5,2) = componentR(5,2);
component(6,2) = componentR(6,2);
component(9,2) = componentR(9,2);
component(12,2) = componentR(12,2);

component(4,4) = componentR(4,4);
component(5,4) = componentR(5,4);
component(6,4) = componentR(6,4);
component(8,4) = componentR(8,4);
component(9,4) = componentR(9,4);
component(11,4) = componentR(11,4);
component(12,4) = componentR(12,4);
component(14,4) = componentR(14,4);


state(2,5) = stateR(2,5); %mdot
state(9,4) = stateR(9,4); %f
state(9,3) = stateR(9,3); %To4
design(2,2) = designR(2,2); %alpha

% mdot4 = state{9,5};
% mdot45 = state{12,5};

M4 = 1;
M45 = 1;
M6A = M6AR;
M8 = M8R;

[state,design] = derived_parameters_performance(state,inputs,design,component);

%% Initial calcs
tau_m1 = component{8,4};
tau_m2 = component{11,4}; 
tau_tH = component{9,4};
tau_tL = component{12,4}; 
tau_MR = component{14,4};
f = state{9,4};
alpha = design{2,2};
beta = design{3,2};
ep1 = design{4,2};
ep2 = design{5,2};


[state] = unFAIR3(state,9);
ho4 = state{9,8};

ho45 = ho4*tau_m1*tau_tH*tau_m2;
state(12,8) = {ho45};
[state] = unFAIR3(state,12);

ho5 = ho45*tau_tL;
state(13,8) = {ho5};
f45 = state{12,4};
[state] = unFAIR3(state,13);

ho6A = ho5*tau_MR; %Is this the right eqn?
state(15,8) = {ho6A};
[state] = unFAIR3(state,15);


%% "" Loop "" 1 (will be in its own function) - Calculates states of comPessors

%adjust values at every time loop happens
[state,design] = derived_parameters_performance(state,inputs,design,component); %update w/ every new f or mdot or alpha

tau_cL = component{5,4};
tau_cH = component{6,4};
ho0 = state{3,8};
To4 = state{8,3};
alpha = design{2,2};
f = state{9,4};

ho3 = ho0*tau_cL*tau_cH;
state(7,8) = {ho3};
state(7,2) = {[]};
state(7,3) = {[]};
[state] = unFAIR3(state,7);
state(8,8) = state(7,8);
state(8,2) = {[]};
state(8,3) = {[]};
[state] = unFAIR3(state,8);

alpha_Pime = alpha / ((1+f)*(1-beta-ep1-ep2)+ep1+ep2);

state(9,2) = {[]};
state(9,8) = {[]};
[state] = unFAIR3(state,9); %reset burner w/ new f
f45 = state{12,4};




% M45 = .3;
% M6 = .4;
% [a] = turb_R(stateR,M45,M6)


% M4 = 1;
% 
% 
% A4 = 1;
% A45 = 1;
To45R = state{12,3};
eta_t = component{9,5};



%% Current checkpoint
    % Fix turb and turb c w/ precalcs
    % Eventually rewrite RGCompr
    % Proceed
    
% First find A4, A45, and A6 from reference code (M6 must be known)! Need to try to implement!!!
% [state] = turbc(state,M4,M45,A4,A45,To45R,eta_t); % HP turbine, driven by mach and area constant
% turb needs to be added



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














% 
% function [state, component,v0] = ambient(state,component,inputs)
% alt = inputs{2,2};
% M0 = inputs{3,2};
% [T0, ~, ~, ~] = atmosisa(alt); %obtain standard atmospheric conditions
% state(2,3) = {T0};
% [state] = unFAIR3(state,2);
% [~,~,T0,~,~,cp0,gamma0,~,~] = state{2,:};
% R0 = cp0 - cp0/gamma0;
% a0 = sqrt(R0*gamma0*T0); %[m/s]
% v0 = M0*a0; %[m/s]
% 
% T_o0 = T0*(1+((M0^2)*((gamma0-1)/2))); %find total temperature using isentropic
% state(3,3) = {T_o0};
% [state] = unFAIR3(state,3);
% 
% 
% P0 = state{2,2};
% Po0 = state{3,2};
% pi_r = Po0/P0;
% component{2,2} = pi_r;
% 
% h0 = state{2,8};
% ho0 = state{3,8};
% tau_r = ho0/h0;
% component{2,4} = tau_r;
% end
% function [state,component] = inlet(state,component,inputs)
% M0 = inputs{3,2};
% pi_dmax = component{3,2};
% Po0 = state{3,2};
% 
% if M0<1
% pid = pi_dmax; 
% elseif M0>1 && M0<5
% pid = pi_dmax * (1-.075*((M0-1)^1.35));
% else 
% pid = pi_dmax * (800/(M0^4 + 935));
% end
% 
% Po2 = pid*Po0;
% 
% component{3,2} = [pi_dmax,pid];
% state(4,2) =  {Po2};
% [state] = unFAIR3(state,4);
% 
% ho0 = state{3,8};
% ho2 = state{4,8};
% tau_d = ho2/ho0;
% component{3,4} = tau_d;
% end
% 









%% Functions
%Parametric Functions
function [state,component,performance] = component_seperate(state,component,design,inputs)
% Runs an engine analysis w/ a seperated LP and HP spools
[state, component,v0] = ambient(state,component,inputs);
[state,component] = inlet(state,component,inputs);
[state,component] = fan(state,component);
[state,component] = LPcomp(state,component);
[state,component] = HPcomp(state,component);
[state,component] = burner(state,component);
[state,component] = HPturb(state,component,design);
[state,component] = LPturb(state,component,design);
[state,component] = mixer(state,component);
[state,component,performance] = nozzle(state,component,inputs,v0,design);
fprintf('%s\n\n','This analysis was completed using SEPERATE high and low spools.')
end
function [state,component,performance] = component_combined(state,component,design,inputs)
% Runs an engine analysis w/ a combined LP and HP spools
[state, component,v0] = ambient(state,component,inputs);
[state,component] = inlet(state,component,inputs);
[state,component] = fan(state,component);
[state,component] = combinedcomp(state,component);
[state,component] = burner(state,component);
[state,component] = combinedturb(state,component,design);
[state,component] = mixer(state,component);
[state,component,performance] = nozzle(state,component,inputs,v0,design);
fprintf('%s\n\n','This analysis was completed using COMBINED high and low spools.')
end
function [state,design] = derived_parameters_parametric(state,inputs,design,component)
%% Derived Parameters
% COOLING AIR CALCULATIONS
% Found in Aircraft Engine Design - Mattingly but
% unable to find rationale...

alt = inputs{2,2};
M0 = inputs{3,2};
F_mdot = inputs{4,2};
mdot = inputs{5,2};
S = inputs{6,2};
T_t4 = inputs{7,2};
alpha = design{2,2};
beta = design{3,2};
eta_b = component{7,5};


if T_t4 > 2400*.55556
    ep1 = (T_t4/.5556-2400)/(16000);
    ep2 = ep1;
else
    ep1 = 0;
    ep2 = 0;
end

design(4,2) = {ep1};
design(5,2) = {ep2};


% MASS FLOW AND FUEL TO AIR CALCULATIONS
mdot_f = S*F_mdot*mdot /eta_b; %mass flow per fuel/air ratio

f0 = 0; %freestream fuel/air ratio

mdot25 = mdot/(1+alpha); %after bypass leaves
mdot13 = mdot - mdot25; % bypass mass flow

mdot31 = mdot25*(1-beta - ep1 -ep2); %after bleed and coolant leaves
mdotbeta = mdot25*beta; %bleed air
mdotep1 = mdot25*ep1; %coolant air 1
mdotep2 = mdot25*ep2; %coolant air 2
mdotep = mdotep1+mdotep2; %combined coolant air

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
function [state,component] = fan(state,component)
pif = component{4,2};
ef = component{4,3};
Po2 = state{4,2};

Po13 = Po2*pif^(1/ef);
state(5,2) = {Po13};
[state] = unFAIR3(state,5);


ho2 = state{4,8};
ho13 = state{5,8};
tauf = ho13/ho2;
component{4,4} = tauf;
end
function [state,component] = LPcomp(state,component)

picl = component{5,2};
ecl = component{5,3};
Po2 = state{4,2};

Po25 = Po2*picl^(1/ecl);
state(6,2) = {Po25};
[state] = unFAIR3(state,6);

ho2 = state{4,8};
ho25 = state{6,8};
taucl = ho25/ho2;
component{5,4} = taucl;
end
function [state,component] = HPcomp(state,component)
pich = component{6,2};
ech = component{6,3};
Po25 = state{6,2};

Po3 = Po25*pich^(1/ech);
state(7,2) = {Po3};
[state] = unFAIR3(state,7);

ho25 = state{6,8};
ho3 = state{7,8};
taucl = ho3/ho25;
component{6,4} = taucl;
end
function [state,component] = combinedcomp(state,component)
picl = component{4,2};
pich = component{5,2};
pic = picl*pich;

ecl = component{4,3};
ech = component{5,3};
ec = (ecl+ech)/2;

Po2 = state{4,2};

Po3 = Po2*pic^(1/ec);
state(7,2) = {Po3};
[state] = unFAIR3(state,7);
end
function [state,component] = burner(state,component)
state(8,2:3) = state(7,2:3);
state(8,6:12) = state(7,6:12);
[state] = unFAIR3(state,9);

Pro31 = state{8,2};
Pro4  = state{9,2};
pi_b_total = Pro4 / Pro31;
component(7,2) = { pi_b_total};

ho31 = state{7,8};
ho4 = state{9,8};
taub = ho4/ho31;
component{7,4} = taub;
end
function [state,component] = HPturb(state,component,design)
mdot3 = state{7,5};
mdot4 = state{9,5};
mdot41 = state{10,5};
ho25 = state{6,8};
ho3 = state{7,8};
hoep1 = state{8,8};
ho4 = state{9,8};
etamH = component{9,5};
etamPH = component{10,5};
PtoH = design{7,2};
mdotep1 = state{21,5};

%Across mixer
ho41 = (mdotep1*hoep1 + mdot4*ho4) / (mdot41);
state(10,8) = {ho41};
[state] = unFAIR3(state,10);

taum1 = ho41/ho4;
component{8,4} = taum1;
pim1 = state{10,2} / state{9,2};
component{8,2} = pim1;

%Across turbine
% etH = .89;
fun = @(ho44) mdot41*(ho41-ho44)*etamH... %change in energy across HPturb
    -mdot3*(ho3-ho25)...                    %change in energy across HP compressor
    -(PtoH) / etamPH;                         %energy draw of takeoff power
ho44 = fzero(fun,ho41);

state(11,8) = {ho44};
[state] = unFAIR3(state,11);

tauth = ho44/ho41;
component{9,4} = tauth;
pitH = state{11,2} / state{10,2};
component{9,2} = pitH;
end
function [state,component] = LPturb(state,component,design)
mdot13 = state{5,5};
mdot25 = state{6,5};
mdot5 = state{13,5};
mdot44 = state{11,5};
mdot45 = state{12,5};
ho2 = state{4,8};
ho13 = state{5,8};
ho25 = state{6,8};
hoep2 = state{8,8};
ho44 = state{11,8};
etamL = component{12,5};
etamPL = component{13,5};
PtoL = design{6,2};
mdotep2 = state{22,5};

%Across mixer
ho45 = (mdotep2*hoep2 + mdot44*ho44) /(mdot45);
state(12,8) = {ho45};
[state] = unFAIR3(state,12);

taum2 = ho45/ho44;
component{11,4} = taum2;
pim2 = state{12,2} / state{11,2};
component{11,2} = pim2;

%Across Turbine
fun = @(ho5) mdot5*(ho45 - ho5)*etamL...    %change in energy across LP turb
    -mdot25*(ho25-ho2)...                   %change in energy across LP compressor
    -mdot13*(ho13-ho2)...                    %change in energy across fan
    -PtoL / etamPL;                         %energy draw of takeoff power
ho5 = fzero(fun,ho45);

state(13,8) = {ho5};
[state] = unFAIR3(state,13);

tautl = ho5/ho45;
component{12,4} = tautl;
pitL = state{13,2} / state{12,2};
component{12,2} = pitL;
end
function [state,component] = combinedturb(state,component,design)
hoep = state{8,8};
ho4 = state{9,8};
ho25 = state{6,8};
ho13 = state{5,8};
ho2 = state{4,8};
mdot13 = state{5,5};
mdot25 = state{6,5};
mdot4 = state{9,5};
mdot45 = state{12,5};
mdotep = state{20,5};

ho45 = (mdotep*hoep + mdot4*ho4) /(mdot45);
state(12,8) = {ho45};
[state] = unFAIR3(state,12);


etamL = component{12,5};
etamPL = component{13,5};
PtoL = design{6,2};
etamH = component{9,5};
etamPH = component{10,5};
PtoH = design{7,2};


etam = (etamH+etamL)/2;
etamP = (etamPH+etamPL)/2;

mdot3 = state{7,5};
ho3 = state{7,8};
ho2 = state{4,8};
Pto = PtoH + PtoL;

fun = @(ho5) mdot45*(ho4-ho5)*etam... %change in energy across turbine
    -mdot3*(ho3-ho13)...              %change in energy across comp
    -mdot13*(ho13-ho2)...             %change in energy across fan
    -(PtoL) / etamPL...               %energy draw off LP spool
    -(PtoH) / etamPH;                 %energy draw off HP spool

ho5 = fzero(fun,ho45);
state(13,8) = {ho5};
[state] = unFAIR3(state,13);

taut = ho5/ho45;
component{12,4} = taut;
component{9,4} = taut;
pit = state{13,2} / state{12,2};
component{12,2} = pit;
component{9,2} = pit;
end
function [state,component] = mixer(state,component)
%Assume ideal pressure ratio (mach independent)
%Assume perfect polytropic efficiency

piM = component{14,2};
hoalpha = state{5,8};
ho5 = state{13,8};
mdotalpha = state{5,5};
mdot5 = state{13,5};
mdot6 = state{14,5};

ho6 = (mdotalpha*hoalpha + mdot5*ho5) / (mdot6);
state(14,8) = {ho6};
[state] = unFAIR3(state,14);

Po6 = state{14,2};
Po6A = Po6*piM;
state(15,2) = {Po6A};
[state] = unFAIR3(state,15);

tauM = ho6/ho5;
component{14,4} = tauM;
end
function [state,component] = afterburner(state,component)
[state] = unFAIR3(state,16);

Pro31 = state{8,2};
Pro4  = state{9,2};
pi_ab_total = Pro4 / Pro31;
component(7,2) = { pi_ab_total};

ho31 = state{7,8};
ho4 = state{9,8};
tauab = ho4/ho31;
component{7,4} = tauab;
end
function [state,component,performance] = nozzle(state,component,inputs,v0,design)
alpha = design{2,2};
beta = design{3,2};
h_PR = design{8,2};
PtoL = design{6,2};
PtoH = design{7,2};
Po9_P9 = inputs{8,2};

state(16,2:end) = state(15,2:end); %assume no afterburner

%Calculate pressure drop across nozzle
Pro7 = state{16,2};
pin = component{16,2};
Pro9 = Pro7*pin;
state(17,2) = {Pro9};
[state] = unFAIR3(state,17);

%Calculate static conditions of exhaust
Pr9 = Pro9 / Po9_P9;
state(18,2) = {Pr9};
component{16,4} = state{17,8} / state{16,8};

[state] = unFAIR3(state,18);

[~,~,~,~,~,~,~,ho9] = state{17,:};
[~,~,T9,~,~,cp9,gamma9,h9] = state{18,:};
R9 = cp9 - cp9/gamma9;
a9 = sqrt(R9*gamma9*T9); %m/s
v9 = sqrt(2*(ho9-h9));
M9 = v9 / a9;

[~,Pr0,~,~,~,cp0,gamma0,~] = state{2,:};
R0 = cp0 - cp0/gamma0;
mdot0 = state{2,5};
f_0 = state{18,4};




F_mdot = (1+f_0-(beta/(1+alpha)))*v9     -   v0  +   (1+f_0-(beta/(1+alpha)))*R9*T9*(1-Pr0/Pr9)/(R0*v9*gamma0);
S = f_0 / F_mdot;
eta_TH = ((v0^2/2*((1+f_0-(beta/(1+alpha)))*(v9/v0)^2 - 1) + (PtoL + PtoH)/mdot0))/...
    (f_0*h_PR);
%eta_P = 2/(1+v9/v0); Simplified case, neglected for now
eta_P = (2*F_mdot/v0)/...
    ((1+f_0-beta/(1+alpha))*(v9/v0)^2 - 1);
eta_o = eta_TH*eta_P;
% eta_o = v0/(h_PR*S); Alternate eqn
performance(1,:) = {'Thrust','Specific Fuel Consumption','Propulsive Efficiency','Thermal Efficiency','Overall Efficiency','Exhaust Mach'};
performance(2,:) = {F_mdot,S,eta_TH,eta_P,eta_o,M9};
end

%Performance Functions
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
