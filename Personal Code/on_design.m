function [state,component,design, inputs, performance] = on_design(alt,M0,pi_f,pi_cL,pi_cH,alpha,beta,PtoH,PtoL,A0,year)
%% Read Me
% This function takes in 11 inputs for turbofan efficiency parameters and
% calculates component efficiencies and on-design performance

%% Setup
% Establish cells
state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
state(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o19';'19';'~';'o9';'9';'beta';'eptot';'ep1';'ep2'};
component = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
component(2:17,1) = {'Ram Recovery';'Inlet(Ideal,Actual)';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
design = {'Parameter','Value'};
design(2:8,1) = {'alpha';'beta';'epsilon1';'epsilon2';'PtoL';'PtoH';'h_PR'};
inputs = {'Parameter','Value'};
inputs(2:8,1) = {'Altitude (m)','Mach Number','F/mdot','Mass Flow Rate (kg/s)','SFC (kg/s/N)','Max Burner Temp (K)','Inlet Area (m^2)'};

% Assume values come from level of tech
[pi_dmax, e_c, e_f, pi_b, eta_b, e_t, pi_n, T_t4] = LOT(year);

% Calculate coolant air
ep1 = (T_t4 - 1333.333333333)/8888.8888888888889;
if ep1 < 0
    ep1 = 0;
end
ep2 = ep1;


% Assume 100% mechanical efficiency on shaft and shaft takeoff power
etamH = 1;
etamPH = 1;
etamL = 1;
etamPL = 1;

% Assume a CH2n propellant
h_PR = 18400*2326; %J/kg, for a (CH2)n propellant

%Store values
state(9,3) = {T_t4};
state(9,3) = {1294};
component(3,2) = {pi_dmax}; %store values in component
component(4,2:3) = {pi_f,e_f};
component(5,2:3) = {pi_cL,e_c};
component(6,2:3) = {pi_cH,e_c};
component(7,2) = {pi_b};
component(7,5) = {eta_b};
component(9,5) = {etamH};
component(10,5) = {etamPH};
component(9,3) = {e_t};
component(12,5) = {etamL};
component(13,5) = {etamPL};
component(12,3) = {e_t};
component(16,2) = {pi_n};
component(2:3,3) = {1};
component(7:8,3) = {1};
component(11,3) = {1};
component(14:16,3) = {1};
design(2,2) = {alpha}; %store values in design
design(3,2) = {beta};
design(4,2) = {ep1};
design(5,2) = {ep2};
design(6,2) = {PtoL};
design(7,2) = {PtoH};
design(8,2) = {h_PR};
inputs(2,2) = {alt};
inputs(3,2) = {M0};
inputs(7,2) = {T_t4};
inputs(8,2) = {A0};

%% Calculations
% Step through the engine station by station modeling the change of state
[state, component,v0] = ambient(state,component,design,inputs);
[state,component] = inlet(state,component,inputs);
[state,component] = fan(state,component);
[state,component] = LPcomp(state,component);
[state,component] = HPcomp(state,component);
[state,component] = burner(state,component,design);
[state,component] = HPturb(state,component,design);
[state,component] = LPturb(state,component,design);
[state,component,performance] = nozzle(state,component,v0,design);

end

%% Functions
function [state,design] = mdot(state,design)
% Calculates mass flow throughout engine sections

% Input Parameters
mdot0 = state{2,5};
f0 = 0; %freestream fuel/air ratio
alpha = design{2,2};
beta = design{3,2};
ep1 = design{4,2};
ep2 = design{5,2};

% Calculate mass flow and fuel to air ratio at every station
f = state{9,4};
mdot4 = state{9,5};
mdot31 = state{8,5};
if size(f,1) == 1
    mdot_f = mdot4 - mdot31;
else
    mdot_f = 0;
end

mdot25 = mdot0/(1+alpha); %after bypass leaves
mdot13 = mdot0 - mdot25; % bypass mass flow

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

% mdot6A = mdot45 + mdot13; %mass flow rate after addtion of bypass air
% f6A = f45*mdot45 / mdot6A; %fuel/air ratio after addtion of bypass air

% Store all values
state(2:8,4) = {f0};
state(2:4,5) = {mdot0};
state(5,5) = {mdot13};
state(6:7,5) = {mdot25};
state(8,5) = {mdot31};
state(9,4) = {f4};
state(9,5) = {mdot4};
state(10:11,4) = {f41};
state(10:11,5) = {mdot41};
state(12:13,4) = {f45};
state(12:13,5) = {mdot45};
state(14:15,4) = {f0};
state(14:15,5) = {mdot13};
state(17:18,4) = {f45};
state(17:18,5) = {mdot45};

state(19,5) = {mdotbeta};
state(20,5) = {mdotep};
state(21,5) = {mdotep1};
state(22,5) = {mdotep2};
end
function [state, component,v0] = ambient(state,component,design,inputs)
% Calculates the thermodynamic state at freestream ambient and stagnation
% ambient

% Input Parameters
alt = inputs{2,2};
M0 = inputs{3,2};
A0 = inputs{8,2};

% Model freestream ambient thermo state based on altitude at standard atmospheric model
[T0, ~, ~, rho0] = atmosisa(alt);
state(2,3) = {T0};
state(2,4) = {0};
[state] = unFAIR3(state,2);

% Calculate inlet velocity
[~,~,T0,~,~,cp0,gamma0,~,~] = state{2,:};
R0 = cp0 - cp0/gamma0;
a0 = sqrt(R0*gamma0*T0); %[m/s]
v0 = M0*a0; %[m/s]

% Calculate mass flow rate through the engine and calculate the subsequent
% mass flow and fuel to air ratio at every station that follows. Note: This
% assumes a wrong fuel to air ratio at the burner and is overriden whenever
% the burner fuel to air ratio is found.
mdot0 = A0 * v0 * rho0;
state{2,5} = mdot0;
[state,design] = mdot(state,design);

% Calculate stagnation ambient thermo conditions
T_o0 = T0*(1+((M0^2)*((gamma0-1)/2))); %find total temperature using isentropic
state(3,3) = {T_o0};
[state] = unFAIR3(state,3);

% Find and store ram pressure ratio
P0 = state{2,2};
Po0 = state{3,2};
pi_r = Po0/P0;
component{2,2} = pi_r;

% Find and store ram enthalpy ratio
h0 = state{2,8};
ho0 = state{3,8};
tau_r = ho0/h0;
component{2,4} = tau_r;
end
function [state,component] = inlet(state,component,inputs)
% Calculates the thermodynamic state at stagnation post-inlet

% Input Parameters
pi_d = component{3,2};
Po0 = state{3,2};

% Calculate stagnation post-inlet thermo conditions
Po2 = pi_d*Po0;
state(4,2) =  {Po2};
[state] = unFAIR3(state,4);

% Find and store ram enthalpy ratio
ho0 = state{3,8};
ho2 = state{4,8};
tau_d = ho2/ho0;
component{3,4} = tau_d;
end
function [state,component] = fan(state,component)
% Calculates the thermodynamic state at stagnation post-fan

% Input Parameters
pi_f = component{4,2};
e_f = component{4,3};
Po2 = state{4,2};

% Calculate stagnation post-fan thermo conditions
Po13 = Po2*pi_f^(1/e_f);
state(5,2) = {Po13};
[state] = unFAIR3(state,5);

% Find and store fan enthalpy ratio
ho2 = state{4,8};
ho13 = state{5,8};
tauf = ho13/ho2;
component{4,4} = tauf;

% Find and store fan efficiency
state_ideal = state;
Po13i = Po2*pi_f; %mechanical efficiency
state_ideal(5,2) = {Po13i};
state_ideal(5,3) = {[]};
state_ideal(5,8) = {[]};
[state_ideal] = unFAIR3(state_ideal,5);
ho13i = state_ideal{5,8};
etaf = (ho13i-ho2)/(ho13-ho2);
component{4,5} = etaf;
end
function [state,component] = LPcomp(state,component)
% Calculates the thermodynamic state at stagnation post-LowPressureCompressor

% Input Parameters
picl = component{5,2};
ecl = component{5,3};
Po2 = state{4,2};

% Calculate stagnation post-LPComp thermo conditions
Po25 = Po2*picl^(1/ecl);
state(6,2) = {Po25};
[state] = unFAIR3(state,6);

% Find and store LPComp enthalpy ratio
ho2 = state{4,8};
ho25 = state{6,8};
taucl = ho25/ho2;
component{5,4} = taucl;

% Find and store LPComp efficiency
state_ideal = state;
Po25i = Po2*picl;
state_ideal(6,2) = {Po25i};
state_ideal(6,3) = {[]};
state_ideal(6,8) = {[]};
[state_ideal] = unFAIR3(state_ideal,6);
ho25i = state_ideal{6,8};
etacL = (ho25i-ho2)/(ho25-ho2);
component{5,5} = etacL;
end
function [state,component] = HPcomp(state,component)
% Calculates the thermodynamic state at stagnation post-HighPressureCompressor

% Input Parameters
pich = component{6,2};
ech = component{6,3};
Po25 = state{6,2};

% Calculate stagnation post-HPComp thermo conditions
Po3 = Po25*pich^(1/ech);
state(7,2) = {Po3};
[state] = unFAIR3(state,7);

% Find and store HPComp enthalpy ratio
ho25 = state{6,8};
ho3 = state{7,8};
tauch = ho3/ho25;
component{6,4} = tauch;

% Find and store HPComp efficiency
state_ideal = state;
Po3i = Po25*pich;
state_ideal(7,2) = {Po3i};
state_ideal(7,3) = {[]};
state_ideal(7,8) = {[]};
[state_ideal] = unFAIR3(state_ideal,7);
ho3i = state_ideal{7,8};
etach = (ho3i-ho25)/(ho3-ho25);
component{6,5} = etach;
end
function [state,component] = burner(state,component,design)
% Calculates the thermodynamic state at stagnation pre burner and
% stagnation post burner

% Model state after bypass and bleed air is removed
state(8,2:3) = state(7,2:3);
state(8,6:12) = state(7,6:12);
[state] = unFAIR3(state,9);

% Input Parameters
eta_b = component{7,5};
mdot31 = state{8,5};
ho31 = state{7,8};
ho4 = state{9,8};
h_PR = design{8,2};
f_i = .03;
error = 1;

% Iterate until fuel to air ratio is found. Since fuel to air ratio affects
% the enthalpy ratio for equivalent temperature ratios, this process takes
% several iterations
while error > .001
    state{9,2} = [];
    state{9,4} = f_i;
    state{9,8} = [];
    [state] = unFAIR3(state,9);
    ho4 = state{9,8};

    mdotf = mdot31*(ho4 - ho31) / (h_PR*eta_b - ho4);
    mdot4 = mdot31 + mdotf;
    state{9,5} = mdot4;
    f = mdotf/mdot31;
    error = (f - f_i)/f_i;
    f_i = f;
end

% Find and store burner enthalpy ratio
taub = ho4/ho31;
component{7,4} = taub;
[state,design] = mdot(state,design);
end
function [state,component] = HPturb(state,component,design)
% Calculates the thermodynamic state at stagnation pre HighPressureTurbine and
% stagnation post HighPressureTurbine

% Input Parameters
mdot3 = state{7,5};
mdot4 = state{9,5};
mdot41 = state{10,5};
ho25 = state{6,8};
ho3 = state{7,8};
hoep1 = state{8,8};
ho4 = state{9,8};
eth = component{9,3};
etamH = component{9,5};
etamPH = component{10,5};
PtoH = design{7,2};
mdotep1 = state{21,5};

% Calculate stagnation pre-HPTurb/post-mixer1 thermo conditions
ho41 = (mdotep1*hoep1 + mdot4*ho4) / (mdot41);
state(10,8) = {ho41};
[state] = unFAIR3(state,10);

% Find and store mixer1 enthalpy ratio
taum1 = ho41/ho4;
component{8,4} = taum1;
Po41 = state{10,2};

% Calculate stagnation post-HPTurb thermo conditions. This involves a power
% balance on the HP spool to find enthalpy change
fun = @(ho44) mdot41*(ho41-ho44)*etamH...       %change in energy across HPturb
    -mdot3*(ho3-ho25)...                        %change in energy across HP compressor
    -(PtoH) / etamPH;                           %energy draw of takeoff power
ho44 = fzero(fun,ho41);
state(11,8) = {ho44};
[state] = unFAIR3(state,11);

% Find and store HPTurb enthalpy ratio
tauth = ho44/ho41;
component{9,4} = tauth;

% Find and store HPTurb pressure ratio
Po44 = state{11,2};
pitH = (Po44 / Po41)^eth;
component{9,2} = pitH;

% Find and store HPTurb efficiency
state_ideal = state;
Po44i = Po41*pitH; %mechanical efficiency
state_ideal(11,2) = {Po44i};
state_ideal(11,3) = {[]};
state_ideal(11,8) = {[]};
[state_ideal] = unFAIR3(state_ideal,11);
ho44i = state_ideal{11,8};
etatH = (ho41-ho44i)/(ho41-ho44);
component{9,6} = etatH;
end
function [state,component] = LPturb(state,component,design)
% Calculates the thermodynamic state at stagnation pre LowPressureTurbine and
% stagnation post LowPressureTurbine

% Input Parameters
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

% Calculate stagnation pre-LPTurb/post-mixer2 thermo conditions
ho45 = (mdotep2*hoep2 + mdot44*ho44) /(mdot45);
state(12,8) = {ho45};
[state] = unFAIR3(state,12);
Po45 = state{12,2};

% Find and store mixer2 enthalpy ratio
taum2 = ho45/ho44;
component{11,4} = taum2;

% Calculate stagnation post-HPTurb thermo conditions. This involves a power
% balance on the HP spool to find enthalpy change
fun = @(ho5) mdot5*(ho45 - ho5)*etamL...    %change in energy across LP turb
    -mdot25*(ho25-ho2)...                   %change in energy across LP compressor
    -mdot13*(ho13-ho2)...                   %change in energy across fan
    -PtoL / etamPL;                         %energy draw of takeoff power
ho5 = fzero(fun,ho45);
state(13,8) = {ho5};
[state] = unFAIR3(state,13);
Po5 = state{13,2};

% Find and store HPTurb enthalpy ratio
tautl = ho5/ho45;
component{12,4} = tautl;

% Find and store HPTurb pressure ratio
etl = component{12,3};
pitL = (Po5 / Po45)^etl;
component{12,2} = pitL;

% Find and store HPTurb efficiency
state_ideal = state;
Po5i = Po45*pitL;
state_ideal(13,2) = {Po5i};
state_ideal(13,3) = {[]};
state_ideal(13,8) = {[]};
[state_ideal] = unFAIR3(state_ideal,13);
ho5i = state_ideal{13,8};
etatH = (ho45-ho5i)/(ho45-ho5);
component{12,6} = etatH;
end
function [state,component,performance] = nozzle(state,component,v0,design)
% Calculates the thermodynamic state at stagnation bypass-exhaust,
% stagnation core-exhaust, static bypass-exhaust, and static core-exhaust 

% Input Parameters
h_PR = design{8,2};
PtoL = design{6,2};
PtoH = design{7,2};
mdot0 = state{2,5};
[~,pi_r,pi_d,pi_f,pi_cL,pi_cH,pi_b,~,pi_tH,~,~,pi_tL,~,~,~,pi_n,~] = component{:,2};

% Calculate stagnation bypass-exhaust thermo conditions
Pro13 = state{5,2};
Pro19 = Pro13*pi_n;
state(14,2) = {Pro19};
[state] = unFAIR3(state,14);

% Calculate static bypass-exhaust thermo conditions
            Po19_P19 = pi_r*pi_d*pi_f*pi_n;
Pr19 = Pro19 / Po19_P19;
state(15,2) = {Pr19};
[state] = unFAIR3(state,15);

% Calculate bypass-exhaust velocity
[~,~,T19,f19,mdot19,~,gamma19,~,~,R19,~,~] = state{15,:};
M19 = sqrt((( Po19_P19^((gamma19-1)/gamma19)  -  1 ) / ( (gamma19 - 1)/2)));
a19 = sqrt(R19*gamma19*T19); %m/s
v19 = M19*a19;


% Calculate stagnation core-exhaust thermo conditions
Pro5 = state{13,2};
Pro9 = Pro5*pi_n;
state(17,2) = {Pro9};
[state] = unFAIR3(state,17);

% Calculate static core-exhaust thermo conditions
            Po9_P9 = pi_r*pi_d*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_n;
Pr9 = Pro9 / Po9_P9;
state(18,2) = {Pr9};
[state] = unFAIR3(state,18);

% Calculate core-exhaust velocity
[~,~,T9,f9,mdot9,~,gamma9,~,~,R9,~,~] = state{18,:};
M9 = sqrt((( Po9_P9^((gamma9-1)/gamma9)  -  1 ) / ((gamma9 - 1)/2)));
a9 = sqrt(R9*gamma9*T9); %m/s
v9 = M9*a9;

%Calculate fuel flow
mdot4 =  state{9,5};
f = state{9,4};
mdotf = mdot4*f/(1+f);

%Calculate performance parameters
F = mdot19*v19 +mdot9*v9 - mdot0*v0; %Thrust
S = mdotf / F; %Specific fuel consumption
mech_power = -.5*mdot0*v0^2 + .5*mdot19*v19^2 + .5*mdot9*v9^2 + PtoH + PtoL; %Mechanical power consumed
thrust_power = F*v0; %Thrust power generated
chem_power = mdotf*h_PR; %Chemical power released
eta_o = thrust_power/chem_power; %Overall engine efficiency
eta_th = mech_power/chem_power; %Thermal engine efficiency
eta_p = thrust_power/mech_power; %Propulsive engine efficiency

% Alternative equations
% eta_o = v0/(h_PR*S);
% % eta_o = F*v0 / ((state{9,5} - state{8,5})*h_PR);
% eta_p = (F*v0) / (.5*mdot19*v19^2 + .5*mdot9*v9^2 - .5*mdot0*v0^2 + PtoH + PtoL);
% eta_th = eta_o/eta_p;

% Store parameters
performance(1,:) = {'Thrust (N)','Specific Fuel Consumption (kg/N-s)','Propulsive Efficiency','Thermal Efficiency','Overall Efficiency','Bypass Exhaust Mach','Core Flow Mach'};
performance(2,:) = {F,S,eta_th,eta_p,eta_o,M19,M9};
end
