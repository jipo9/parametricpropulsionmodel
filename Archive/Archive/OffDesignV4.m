%% Off Design V4
clc
clear all

%% On Design Analysis
[stateR,componentR,designR,performanceR,inputsR] = ondesign();


%% Initial Guesses
f = stateR{9,4};
alpha = designR{2,2};
alpha = .53;
beta = designR{3,2};
ep1 = designR{4,2};
ep2 = designR{5,2};
alphap = alpha/((1-beta-ep1-ep2)*(1+f) + ep1 + ep2);
M6 = .3835;
M8 = .5;
mdot0 = stateR{2,5};
mdot0 = 85.602;

%% Independent Variables

invar = {'Parameter','Value'};
invar(2:7,1) = {'M0','T0','P0','To4','To7','P0/P9'};
invar{2,2} = inputsR{3,2};
invar{3,2} = stateR{2,3};
invar{4,2} = stateR{2,2};
invar{5,2} = stateR{9,3};
invar{6,2} = stateR{16,3};
invar{7,2} = 1;

%% Constant Declarations
beta = designR{3,2};
pi_dmax = .96;
eta_f = componentR{4,5};
eta_cL = componentR{5,5};
eta_cH = componentR{6,5};
eta_b = componentR{7,5};
ep1 = designR{4,2};
ep2 = designR{5,2};
eta_tH = componentR{9,5};
eta_tL = componentR{12,5};
pi_Mmax = componentR{14,2};
A6 = 1; %placeholder
A16 = 1; %placeholder
A6A = 1; %placeholder
pi_AB = .95;
pi_n = componentR{16,2};
A8dry = 1; %placeholder
eta_mPL = 1; %not sure where to get these
eta_mPH = 1; %""

%constants from on-design
eta_cL = .8674;
eta_mL = .995;
eta_f = .8674;
pi_b = .95;
Cp_c = 1004.83; 
Cp_t = 1235.11;
eta_AB = .99;
eta_cH = .8751;
eta_mH = .995;
eta_b = .999;
gamma_c = 1.4;
Cp_AB = 1235.11;
eta_tH = .8995;
eta_PL = 1;
pi_c = 20;
pi_n = .97;
gamma_t = 1.3;
gamma_AB = 1.3;
eta_tL = .9074;
eta_PH = 1;
PtoL = designR{6,2};
PtoH = designR{7,2};




%% Off Design Start

[state,component,design,inputs] = offdesign();

alt = inputs{2,2};
M0 = inputs{3,2};
[T0, ~, P0, ~] = atmosisa(alt); %obtain standard atmospheric conditions
state(2,3) = {T0};
state{2,4} = 0;
[state] = unFAIR3(state,2);
V0 = M0*sqrt(state{2,7}*state{2,10}*T0);

ho0 = state{2,8}+(V0^2/2);
state{3,8} = ho0;
state{3,4} = 0;
[state] = unFAIR3(state,3);

h0 = state{2,8};
tau_r = ho0/h0;
Pro0 = state{3,2};
Pr0 = state{2,2};
pi_r = Pro0/Pr0;

if M0<1
    pi_d = pi_dmax;
elseif M0>1 && M0<5
    pi_d = pi_dmax * (1-.075*((M0-1)^1.35));
else
    pi_d = pi_dmax * (800/(M0^4 + 935));
end

state{4,8} = state{3,8};
state{4,2} = state{3,2};
pi_ABdry = 1-(1-pi_AB)/2;



%% Set Initial Values

pi_tH = componentR{9,2};
pi_tL = componentR{12,2};
tau_tH = componentR{9,4};
tau_tL = componentR{12,4};
pi_f = componentR{4,2};
pi_cL = componentR{5,2};
pi_cH = componentR{6,2};
tau_f = componentR{4,4};
tau_cL = componentR{5,4};
tau_cH = componentR{6,4};
tau_m1 = componentR{8,4};
tau_m2 = componentR{11,4};
f = stateR{9,4};
M4 = 1;
M45 = 1;
M6A = .42; %M6A = M6AR, need to get reference machs
M8 = .9;  %M8 = M8R

state{9,4} = f;
[state] = unFAIR3(state,9);
ho4 = state{9,8};
ho45 = ho4*tau_m1*tau_tH*tau_m2;
state{12,8} = ho45;
f45 = f*(1-beta-ep1-ep2)/(1-beta);
state{12,4} = f45;
[state] = unFAIR3(state,12);
ho5 = ho45*tau_tL;
state{13,8} = ho5;
state{13,4} = f45;
[state] = unFAIR3(state,13);

tau_mR = componentR{14,4};
ho6A = ho5*tau_mR;
f6A = f45*((1-beta)/(1+alpha-beta));
state{15,4} = f6A;
state{15,8} = ho6A;
[state] = unFAIR3(state,15);

%1
ho3 = ho0*tau_cL*tau_cH;
state{7,8} = ho3;
state{7,4} = 0;
[state] = unFAIR3(state,7);

alphap = alpha/((1+f)*(1-beta-ep1-ep2)+ep1+ep2);

state{9,2} = [];
state{9,4} = f;
state{9,8} = [];
[state] = unFAIR3(state,9);
f45 = f*(1-beta-ep1-ep2)/(1-beta);
state{12,4} = f45;

% Insert TURB AND TURBC FXs

tau_cLR = componentR{5,2};
tau_fR = componentR{4,4};


tau_lambda = ho4/h0;
tau_f = 1 ...
    + ( (1-tau_tL)*eta_mL* (((1-beta-ep1-ep2)*(1+f)*tau_lambda*tau_tH/tau_r  +  (ep1*tau_tH + ep2)*tau_cL*tau_cH))   -   (1+alpha)*PtoL/(tau_r*eta_mPL*mdot0*h0)) ...
    / ((tau_cLR - 1)/(tau_fR - 1) + alpha) 
tau_cL = 1 + (tau_f -1)*(tau_cLR - 1)/(tau_fR - 1)
tau_cH = 1 ...
    + (1-tau_tH)*eta_mH* (((1-beta-ep1-ep2)*(1+f)*tau_lambda/(tau_r*tau_cL) + ep1*tau_cH))...
    - (1+alpha)*PtoH/(tau_r*tau_cL*eta_mPH*mdot0*h0)

ho2 = ho0;
state{4,8} = state{3,8};
Pro2 = Pro0;
state{4,2} = state{3,2};
ho13 = ho2*tau_f;
state{5,8} = ho13;
ho25 = ho2*tau_cL;
state{6,8} = ho25;
ho3 = ho25*tau_cH;
state{7,8} = ho3;

state{5,4} = 0;
state{6,4} = 0;
state{7,4} = 0;

[state] = unFAIR3(state,5);
[state] = unFAIR3(state,6);
state{7,2} = [];
state{7,3} = [];
[state] = unFAIR3(state,7);


statei = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
statei(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
ho13i = ho2*(1+eta_f*(tau_f-1));
ho25i = ho2*(1+eta_cL*(tau_cL-1));
ho3i = ho25*(1+eta_cH*(tau_cH-1));

statei{5,8} = ho13i;
statei{6,8} = ho25i;
statei{7,8} = ho3i;
statei{5,4} = 0;
statei{6,4} = 0;
statei{7,4} = 0;

[statei] = unFAIR3(statei,5);
[statei] = unFAIR3(statei,6);
[statei] = unFAIR3(statei,7);

Pro13i = statei{5,2};
Pro25i = statei{6,2};
Pro3i = statei{7,2};
Pro2 = state{4,2};
Pro25 = state{6,2};

pi_f = Pro13i/Pro2;
pi_cL = Pro25i/Pro2;
pi_cH = Pro3i/Pro25;
pi_c = pi_cL*pi_cH;
tau_c = tau_cL*tau_cH;

%2
ftemp = f;





%% Functions
function [state,component,design,inputs] = offdesign()
state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
state(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
component = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
component(2:17,1) = {'Ram Recovery';'Inlet(Ideal,Actual)';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
design = {'Parameter','Value'};
design(2:8,1) = {'alpha';'beta';'epsilon1';'epsilon2';'PtoL';'PtoH';'h_PR'};
inputs = {'Parameter','Value'};
inputs(2:8,1) = {'Altitude (m)','Mach Number','F/mdot','Mass Flow Rate (kg/s)','SFC (kg/s/N)','Max Burner Temp (K)','Po9/P9'};

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
To4 = 3200*.555556; %max burner temperature [R to K]
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
design(2,2) = {alpha}; %store values in design
design(3,2) = {beta};
design(7,2) = {PtoH};
design(6,2) = {PtoL};
design(8,2) = {h_PR};
state(9,3) = {To4};
component(3,2) = {pi_dmax}; %store values in component
component(4,2:3) = {pif,ef};
component(5,2:3) = {picL,ecL};
component(6,2:3) = {picH,ecH};
component(7,2) = {pi_b};
component(7,5) = {eta_b};
component(9,5) = {etamH};
component(10,5) = {etamPH};
component(9,3) = {etH};
component(12,5) = {etamL};
component(13,5) = {etamPL};
component(12,3) = {etL};
component(14,2) = {pi_M_max};
component(16,2) = {pin};
component(2:3,3) = {1};
component(7:8,3) = {1};
component(11,3) = {1};
component(14:16,3) = {1};
inputs(2,2) = {alt};
inputs(3,2) = {M0};
inputs(4,2) = {F_mdot};
inputs(5,2) = {mdot};
inputs(6,2) = {S};
inputs(7,2) = {To4};
inputs(8,2) = {Po9_P9};

%[state,design] = derived_parameters_parametric(state,inputs,design,component);


end

function [stateR,componentR,designR,performanceR,inputsR] = ondesign()
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
To4 = 3200*.555556; %max burner temperature [R to K]
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
stateR(9,3) = {To4};
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
inputsR(7,2) = {To4};
inputsR(8,2) = {Po9_P9};

[stateR,designR] = derived_parameters_parametric(stateR,inputsR,designR,componentR);

%prompt = 'Run analysis with combined compressors/turbines? Y/N : ';
%str = input(prompt,'s');
%if str == ('Y') || str == ('y')
    %clc
    %[stateR,componentR,performanceR] = component_combined(stateR,componentR,designR,inputsR);
%elseif str == ('N') || str == ('n')
    %clc
    [stateR,componentR,performanceR] = component_seperate(stateR,componentR,designR,inputsR);
%else
   % error('Invalid Input. Try again and only type Y or N you dummy')
%end

err_T_mdot = performanceR{2,1} /F_mdot; %T/mdot error compared to book
err_s = performanceR{2,2} / S; %S error compared to book
err_efftherm = performanceR{2,3} / .5589; %thermal efficiency error compared to book
err_effprop =performanceR{2,4} / .6162; %propulsive efficiency compared to book

%fprintf('%s%.3f%s\n','Thrust                    of this analysis is ',abs(100*(1-err_T_mdot)),'% off book solution.')
%fprintf('%s%.3f%s\n','Specific Fuel Consumption of this analysis is ',abs(100*(1-err_s)),'% off book solution.')
%fprintf('%s%.3f%s\n','Thermal Efficiency        of this analysis is ',abs(100*(1-err_efftherm)),'% off book solution.')
%fprintf('%s%.3f%s\n','Propulsive Efficiency     of this analysis is ',abs(100*(1-err_effprop)),'% off book solution.')
end
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
[T0, ~, P00, ~] = atmosisa(alt); %obtain standard atmospheric conditions
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
state{23,1} = P00;
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

statet = state;
Po13i = Po2*pif; %mechanical efficiency
state(5,2) = {Po13i};
state(5,3) = {[]};
state(5,8) = {[]};
[state] = unFAIR3(state,5);
ho13i = state{5,8};
state = statet;
etaf = (ho13i-ho2)/(ho13-ho2);
component{4,5} = etaf;
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

statet = state;
Po25i = Po2*picl; %mechanical efficiency
state(6,2) = {Po25i};
state(6,3) = {[]};
state(6,8) = {[]};
[state] = unFAIR3(state,6);
ho25i = state{6,8};
state = statet;
etacL = (ho25i-ho2)/(ho25-ho2);
component{5,5} = etacL;
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

statet = state;
Po3i = Po25*pich; %mechanical efficiency
state(7,2) = {Po3i};
state(7,3) = {[]};
state(7,8) = {[]};
[state] = unFAIR3(state,7);
ho3i = state{7,8};
state = statet;
etach = (ho3i-ho25)/(ho3-ho25);
component{6,5} = etach;
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
Po41 = state{10,2};

%Across turbine
% etH = .89;
fun = @(ho44) mdot41*(ho41-ho44)*etamH... %change in energy across HPturb
    -mdot3*(ho3-ho25)...                    %change in energy across HP compressor
    -(PtoH) / etamPH;                         %energy draw of takeoff power
ho44 = fzero(fun,ho41);

state(11,8) = {ho44};
[state] = unFAIR3(state,11);
Po44 = state{11,2};

tauth = ho44/ho41;
component{9,4} = tauth;
eth = component{9,3};
pitH = (state{11,2} / state{10,2})^(1/eth);
component{9,2} = pitH;

statet = state;
Po44i = Po41*pitH; %mechanical efficiency
state(11,2) = {Po44i};
state(11,3) = {[]};
state(11,8) = {[]};
[state] = unFAIR3(state,11);
ho44i = state{11,8};
state = statet;
etatH = (ho41-ho44)/(ho41-ho44i);
component{9,5} = etatH;

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
Po45 = state{12,2};

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
Po5 = state{13,2};

tautl = ho5/ho45;
component{12,4} = tautl;
etl = component{12,3};
pitL = (state{13,2} / state{12,2})^(1/etl);
component{12,2} = pitL;

Po5i = Po45*pitL; %mechanical efficiency
state(13,2) = {Po5i};
state(13,3) = {[]};
state(13,8) = {[]};
[state] = unFAIR3(state,13);
ho5i = state{13,8};
state(13,2) = {Po5};
state(13,3) = {[]};
state(13,8) = {[]};
[state] = unFAIR3(state,13);
etatH = (ho45-ho5)/(ho45-ho5i);
component{12,5} = etatH;
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



