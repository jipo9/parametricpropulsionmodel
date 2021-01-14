%% Off Design Analysis V3
% Method used in NASA Paper

clc
clear all

%note, gamma changes pis drastically. May want to implement gammas based on
%stage, rather than assume 1.4 and 1.3

%% On-Design Analysis
[stateR,componentR,designR,performanceR,inputsR] = ondesign();


%% Performance Choices

opcond = {'Parameter','Value'};
opcond(2:6,1) = {'M0','alt','T0','P0','To4'};
opcond{2,2} = 1.451; %freestream mach
opcond{3,2} = convlength(36000,'ft','m'); %altitude(ft to m)
[opcond{4,2}, ~, opcond{5,2}, ~] = atmosisa(opcond{3,2});
opcond{6,2} = 1777.777; %To4 [K]

%% Design Constants
dc = {'Parameter','Value'};
dc(2:18,1) = {'gamma_c','gamma_t','Cp_c','Cp_t','pi_dmax','pi_b','pi_n','pi_nf','eta_f','eta_cL','eta_cH','eta_tH','eta_tL','eta_b','eta_mH','eta_mL','hPR'};
dc{2,2} = 1.4;
dc{3,2} = 1.3;
dc{4,2} = 1.0048;
dc{5,2} = 1.2351;
dc{6,2} = componentR{3,2}(1);
dc{7,2} = componentR{7,5};
dc{8,2} = componentR{16,2};
dc{9,2} = 1; %find pi_nf
dc{10,2} = componentR{4,5};
dc{11,2} = componentR{5,5};
dc{12,2} = componentR{6,5};
dc{13,2} = componentR{9,5};
dc{14,2} = componentR{12,5};
dc{15,2} = componentR{7,5};
dc{16,2} = .995;
dc{17,2} = .995;
dc{18,2} = designR{8,2};

%% Reference Conditions
rc = {'Parameter','Value'};
rc(2:25,1) = {'M0','T0','P0','tau_r','pi_r','pi_d','pi_f','pi_cL','pi_cH','pi_tH','pi_tL','tau_f','tau_cL','tau_cH','tau_tH','tau_tL','To4','f_b','M6','M16','alpha','F','mdot0','S'};
rc{2,2} = inputsR{3,2};
rc{3,2} = stateR{2,3};
rc{4,2} = stateR{23,1};
rc{5,2} = 1; %calculated later
rc{6,2} = 1; %calculated later
rc{8,2} = componentR{4,2};
rc{9,2} = componentR{5,2};
rc{10,2} = componentR{6,2};
rc{11,2} = componentR{9,2};
rc{12,2} = componentR{12,2};
rc{13,2} = componentR{4,4};
rc{14,2} = componentR{5,4};
rc{15,2} = componentR{6,4};
rc{16,2} = componentR{9,4};
rc{17,2} = componentR{12,4};
rc{18,2} = stateR{9,3};
rc{19,2} = stateR{9,4};
rc{20,2} = .40; %hardcode from on-design [M6]
rc{21,2} = .394; %hardcode from on-design [M16]
rc{22,2} = designR{2,2};
rc{23,2} = performanceR{2,1};
rc{24,2} = stateR{2,5};
rc{25,2} = performanceR{2,2};

%% Reference Calculations

[rc,dc] = refcalc(rc,dc);


%% Calculations1

M0 = opcond{2,2};
T0 = opcond{4,2};
P0 = opcond{5,2};
gamma_c = dc{2,2};
gamma_t = dc{3,2};
Cp_c = dc{4,2};
Cp_t = dc{5,2};

pi_dmax = dc{6,2};

tau_r = 1 + ((gamma_c-1)/2)*M0^2;
pi_r = tau_r^(gamma_c/(gamma_c-1));

if M0 <= 1
    eta_R = 1;
elseif M0 > 1 && M0 < 5
    eta_R = 1-.075*(M0-1)^1.35;
else
    error('Invalid Mach Number')
end

pi_d = pi_dmax*eta_R;

R_c = 1000*((gamma_c-1)/gamma_c)*Cp_c;
R_t = 1000*((gamma_t-1)/gamma_t)*Cp_t;

a0 = sqrt(gamma_c*R_c*T0);

Po0 = P0*tau_r^(gamma_c/(gamma_c-1));
To0 = T0*tau_r;

[MFP0] = MFP3(M0, gamma_c, R_c);


%% Set Initial Values
mdot0 = rc{24,2};
f_b = rc{19,2};
pi_tH = rc{11,2};
pi_tL = rc{12,2};
tau_tH = rc{16,2};
tau_tL = rc{17,2};
pi_f = rc{8,2};
pi_cL = rc{9,2};
pi_cH = rc{10,2};
tau_f = rc{13,2};
tau_cL = rc{14,2};
tau_cH = rc{15,2};

To4 = rc{18,2};
tau_lam_b = (Cp_t*To4)/(Cp_c*T0);



%% HP Turbine

eta_tH = dc{13,2};
pi_tHR = rc{11,2};
tau_tHR = rc{16,2};

x = 1;
while x>.00001
if x ~= 1 
pi_tH = pi_tHN;
end
tau_tH = 1-eta_tH*(1-pi_tH^(gamma_t-1/gamma_t));
pi_tHN = (sqrt(tau_tH)/sqrt(tau_tHR))*pi_tHR;
x = norm(pi_tHN-pi_tH);
end

%% HP Comp
eta_mH = dc{16,2};
eta_cH = dc{12,2};
f_b = rc{19,2};

tau_cH = 1+(eta_mH*(1+f_b)*((tau_lam_b*(1-tau_tH))/(tau_r*tau_cL)));

pi_cH = 1+eta_cH*(tau_cH-1)^(gamma_c/(gamma_c-1));


%% Main Burner
eta_b = dc{15,2};
hPR = dc{18,2}/1000;
h0 = Cp_c*T0;

f_b = (tau_lam_b-(tau_r*tau_cL*tau_cH))/((eta_b*hPR)/(h0-tau_lam_b));

%% Fan and LP Comp
eta_f = dc{10,2};
eta_cL = dc{11,2};

pi_f = (1+eta_f*(tau_f-1))^(gamma_c/(gamma_c-1));
pi_cL = (1+eta_cL*(tau_cL-1))^(gamma_c/(gamma_c-1));

%% Exhaust Nozzles

P0_P9 = 1;
pi_b = dc{7,2};
pi_n = dc{8,2};

Po9_P9 = (P0_P9)*pi_r*pi_d*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_n;

M9 = sqrt((2/(gamma_t-1))*(Po9_P9^((gamma_t-1)/gamma_t))-1);

M9 = 2.34; %corrected value from printoutB

if M9>1
    M8 = 1;
else
    M8 = M9;
end

[MFP8] = MFP3(M8, gamma_t, R_t);
[MFP9] = MFP3(M9, gamma_t, R_t);

A9_A8 = MFP8/(pi_n*MFP9);

%need to figure out how stations 8,9,18,19 are defined



%% Functions
function [rc,dc] = refcalc(rc,dc)
gamma_c = dc{2,2};
M0R = rc{2,2};
tau_rR = 1 + ((gamma_c-1)/2)*M0R^2;
pi_rR = tau_rR^(gamma_c/(gamma_c-1));
rc{5,2} = tau_rR;
rc{6,2} = pi_rR;

if M0R <= 1
    eta_RspecR = 1;
elseif M0R > 1 && M0R < 5
    eta_RspecR = 1-.075*(M0R-1)^1.35;
else
    error('Invalid Reference Mach Number')
end

pi_dmax = dc{8,2};
pi_dR = pi_dmax*eta_RspecR;
rc{7,2} = pi_dR;
end

function [MFP] = MFP3(M,gamma,R)
x = (gamma+1)/(2*(1-gamma));
s = sqrt(gamma/R);
m = 1+(((gamma-1)/2)*M^2);

MFP = M*s*(m^x);
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

alt = convlength(36000,'ft','m'); %altitude [m from feet]
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

