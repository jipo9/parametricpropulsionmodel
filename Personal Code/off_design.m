function [state,component,design,inputs,performance] = off_design(state,component,design,inputs,M0,alt,A0,A45_9)
%% Read Me
% This function takes in the results from the on-design analysis and
% returns the engine performance for any given flight condition

%% Calculations
% Step through the engine station by station modeling the change of state
% Iteration is used because the off design enthalpy ratios are related to
% each other and take several iterations to converge
inputs{2,2} = alt;
inputs{3,2} = M0;
inputs{8,2} = A0;

tauf_i = component{4,4};
error = 1;
while error > .0001
[state, component,v0] = off_ambient(state,component,design,inputs);
[state,component] = off_inlet(state,component,inputs);
[state,component] = off_fan(state,component,componentR,design);
[state,component] = off_comp(state,component,design, componentR);
[state,component] = off_burner(state,component,design);
[state,component] = off_turb(state,component,componentR);
[state,component,performance] = off_nozzle(state,component,v0,design);
tauf = component{4,4};
error = norm(tauf-tauf_i)/tauf;
tauf_i = tauf;
end
    
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
function [state, component,v0] = off_ambient(state,component,design,inputs)

% Input Parameters
alt = inputs{2,2};
M0 = inputs{3,2};
A0 = inputs{8,2};

% Model freestream ambient thermo state based on altitude at standard atmospheric model
[T0, ~, ~, rho0] = atmosisa(alt); %obtain standard atmospheric conditions
state(2,3) = {T0};
state(2,4) = {0};
state{2,2} = [];
state{2,8} = [];
[state] = unFAIR3(state,2);

% Calculate inlet velocity
[~,~,T0,~,~,cp0,gamma0,~,~] = state{2,:};
R0 = cp0 - cp0/gamma0;
a0 = sqrt(R0*gamma0*T0); %[m/s]
v0 = M0*a0; %[m/s]

% Calculate mass flow rate through the engine and calculate the subsequent
% mass flow and fuel to air ratio at every station that follows, this
% process is iterated and changes w/ f
mdot0 = A0 * v0 * rho0;
state{2,5} = mdot0;
[state,design] = mdot(state,design);

% Calculate stagnation ambient thermo conditions
To0 = T0*(1+((M0^2)*((gamma0-1)/2))); %find total temperature using isentropic
state(3,3) = {To0};
state{3,2} = [];
state{3,8} = [];
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
function [state,component] = off_inlet(state,component,inputs)
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
function [state,component] = off_fan(state,component,componentR,design)
% Calculates the thermodynamic state at stagnation post-fan

% Input Parameters
ho4 = state{9,8};
h0 = state{2,8};
tau_tL = component{12,4};
tau_tH = component{9,4};
tau_cL = component{5,4};
tau_cLR = componentR{5,4};
tau_cH = component{6,4};
eta_mL = .995;
alpha = design{2,2};
beta = design{3,2};
ep1 = design{4,2};
ep2 = design{5,2};
f = state{9,4};
ho0 = state{3,8};
tau_r = ho0/h0;
PtoL = design{6,2};
eta_mPL = 1;
mdot0 = state{2,5};
tau_fR = componentR{4,4};
eta_f = componentR{4,5};

% Calculate the enthalpy ratio from ambient static to burner stagnation
tau_lambda = ho4/h0;

% Find and store fan enthalpy ratio
tau_f = 1 ...
    + ( (1-tau_tL)*eta_mL* (((1-beta-ep1-ep2)*(1+f)*tau_lambda*tau_tH/tau_r  +  (ep1*tau_tH + ep2)*tau_cL*tau_cH))   -   (1+alpha)*PtoL/(tau_r*eta_mPL*mdot0*h0)) ...
    / ((tau_cLR - 1)/(tau_fR - 1) + alpha);
component{4,4} = tau_f;

% Find and store fan pressure ratio
ho2 = state{4,8};
ho13i = ho2*(1+eta_f*(tau_f-1));
statei{5,8} = ho13i;
statei{5,4} = 0;
[statei] = unFAIR3(statei,5);
Pro13i = statei{5,2};
Pro2 = state{4,2};
pi_f = Pro13i/Pro2;
component{4,2} = pi_f;

% Calculate stagnation post-fan thermo conditions
state{5,2} = [];
state{5,3} = [];
ho13 = ho2*tau_f;
state{5,8} = ho13;
[state] = unFAIR3(state,5);
end
function [state,component] = off_comp(state,component,design, componentR)
% Calculates the thermodynamic state at stagnation post-LowPressureCompressor
% and stagnation post-LowPressureCompressor

%% LP Compressor
% Input Parameters
tau_f = component{4,4};
tau_fR = componentR{4,4};
tau_cLR = componentR{5,4};
eta_cL = component{5,5};

% Find and store LPComp enthalpy ratio
tau_cL = 1 + (tau_f -1)*(tau_cLR - 1)/(tau_fR - 1);
component{5,4} = tau_cL;

% Find and store LPComp pressure ratio
ho2 = state{4,8};
ho25i = ho2*(1+eta_cL*(tau_cL-1));
statei{6,8} = ho25i;
statei{6,4} = 0;
[statei] = unFAIR3(statei,6);
Pro25i = statei{6,2};
Pro2 = state{4,2};
pi_cL = Pro25i/Pro2;
component{5,2} = pi_cL;

% Calculate stagnation post-LPComp thermo conditions
state{6,2} = [];
state{6,3} = [];
ho25 = ho2*tau_cL;
state{6,8} = ho25;
[state] = unFAIR3(state,6);

%% HP Compressor
% Input Parameters
tau_r = component{2,4};
tau_tH = component{9,4};
tau_cL = component{5,4};
tau_cH = component{6,4};
eta_cH = component{6,5};
etamH = component{9,5};
etamPH = component{10,5};
alpha = design{2,2};
beta = design{3,2};
ep1 = design{4,2};
ep2 = design{5,2};
PtoH = design{7,2};
f = state{9,4};
h0 = state{2,8};
ho4 = state{9,8};
mdot0 = state{2,5};

% Find and store HPComp enthalpy ratio
tau_lambda = ho4/h0;
tau_cH = 1 ...
    + (1-tau_tH)*etamH* (((1-beta-ep1-ep2)*(1+f)*tau_lambda/(tau_r*tau_cL) + ep1*tau_cH))...
    - (1+alpha)*PtoH/(tau_r*tau_cL*etamPH*mdot0*h0);
component{6,4} = tau_cH;

% Find and store HPComp pressure ratio
ho25 = state{6,8};
ho3i = ho25*(1+eta_cH*(tau_cH-1));
statei{7,8} = ho3i;
statei{7,4} = 0;
[statei] = unFAIR3(statei,7);
Pro3i = statei{7,2};
Pro25 = state{6,2};
pi_cH = Pro3i/Pro25;
component{6,2} = pi_cH;

% Calculate stagnation post-HPComp thermo conditions
state{7,2} = [];
state{7,3} = [];
ho3 = ho25*tau_cL;
state{7,8} = ho3;
[state] = unFAIR3(state,7);
end
function [state,component] = off_burner(state,component,design)
% Calculates the thermodynamic state at stagnation pre burner and
% stagnation post burner

% Model state after bypass and bleed air is removed
state(8,2:3) = state(7,2:3);
state(8,6:12) = state(7,6:12);

% Input Parameters
mdot31 = state{8,3};
h_PR = design{8,2};
f_i = state{9,4};
error = 1;

% Iterate until fuel to air ratio is found. Since fuel to air ratio affects
% the enthalpy ratio for equivalent temperature ratios, this process takes
% several iterations
while error > .001
    state{9,2} = [];
    state{9,8} = [];
    [state] = unFAIR3(state,9);

    ho31 = state{7,8};
    ho4 = state{9,8};

    mdotf = mdot31*(ho4 - ho31) / h_PR; %plus eta_burner
    f = mdotf/mdot31;
    state{9,4} = f;
    [state,design] = mdot(state,design);
    error = (f - f_i)/f_i;
    f_i = f;
end

% Find and store burner enthalpy ratio
taub = ho4/ho31;
component{7,4} = taub;
[state,design] = mdot(state,design);
end
function [state,component] = off_turb(state,component,componentR)
% Calculates the thermodynamic state at stagnation pre and post both
% turbines
% CURRENT ASSUMPTION: CONSTANT TURBINE PERFORMANCE

% Input Parameters
tau_m1 = component{8,4};
tau_tH = componentR{9,4};
tau_tL = componentR{12,4};
tau_m2 = component{11,4};
ho4 = state{9,8};

% Calculate stagnation pre-HPTurb/post-mixer1 thermo conditions
state{10,2} = [];
state{10,3} = [];
ho41 = ho4*tau_m1;
state{10,8} = ho41;
[state] = unFAIR3(state,10);

% Calculate stagnation post-HPTurb thermo conditions
state{11,2} = [];
state{11,3} = [];
ho44 = ho41*tau_tH;
state{11,8} = ho44;
[state] = unFAIR3(state,11);

% Calculate stagnation pre-LPTurb/post-mixer2 thermo conditions
state{12,2} = [];
state{12,3} = [];
ho45 = ho44*tau_m2;
state{12,8} = ho45;
[state] = unFAIR3(state,12);

% Calculate stagnation post-HPTurb thermo conditions
state{13,2} = [];
state{13,3} = [];
ho5 = ho45*tau_tL;
state{13,8} = ho5;
[state] = unFAIR3(state,13);


% Alternative method, use constant area of chocked flow at station 4.5 to find changine HP turbine performance 
% pi_tH = componentR{9,2};
% pi_tL = componentR{12,2};
% gamma = 1.3;
% eta_tL = component{12,6};
% X = A45_9 * 1/M9 * (2/(gamma+1)    *  (1+(gamma-1)/2*M9^2)  )^  ((gamma+1)/(2*(gamma-1)));
% 
% error = 1;
% up = 10;
% low = 0;
% while error > .0000001
%     pi_tL = (up+low)/2;
%     tau_tL = 1 - eta_tL*(1-pi_tL^((gamma-1)/gamma));
%     X_i = pi_tL/sqrt(tau_tL);
%     error = norm(X_i - X)/X;
%     if X_i > X
%        up =  pi_tL;
%     else
%        low =  pi_tL;
%     end
% end     
%         component{8,4} = tau_m1;
%         component{9,4} = tau_tH;
%         component{11,4} = tau_m2;

% component{12,2} = pi_tL;
% component{12,4} = tau_tL;
end
function [state,component,performance] = off_nozzle(state,component,v0,design)
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
state{14,3} = [];
state{14,8} = [];
[state] = unFAIR3(state,14);

% Calculate static bypass-exhaust thermo conditions
Po19_P19 = pi_r*pi_d*pi_f*pi_n;
Pr19 = Pro19 / Po19_P19;
state(15,2) = {Pr19};
state(15,3) = {[]};
state(15,8) = {[]};
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
state{17,3} = [];
state{17,8} = [];
[state] = unFAIR3(state,17);

% Calculate static core-exhaust thermo conditions
Po9_P9 = pi_r*pi_d*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_n;
Pr9 = Pro9 / Po9_P9;
state(18,2) = {Pr9};
state{18,3} = [];
state{18,8} = [];
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
function [To4] = thetabreak(state,inputs)
To4max = 2000*.555555556; %input max To4
gamma = state{2,7};
T0 = state{2,3};
M0 = inputs{3,2};

To4 = To4max;

[Tstd,~,~,~] = atmosisa(0);
To0 = T0*(1+((M0^2)*((gamma-1)/2)));

theta0 = To0/Tstd;

To4sls = To4/theta0;
TR = To4max/To4sls;

M0break = sqrt((TR-1)*(2/(gamma-1)));
while isreal(M0break) == 0
    To4 = To4*.999;
    To4sls = To4/theta0;
    TR = To4max/To4sls;
    M0break = sqrt((TR-1)*(2/(gamma-1)));
end
    
end
