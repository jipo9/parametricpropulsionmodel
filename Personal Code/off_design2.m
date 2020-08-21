function [state,component,design,inputs,performance] = off_design2(state,component,design,inputs,M0,alt,A0)
%% Read Me
% This function takes in the results from the on-design analysis and
% returns the engine performance for any given flight condition

%% Calculations
alphaR = design{2,2};
To4R = state{9,3};
pi_cHR = component{6,2};
pi_fR = component{4,2};
[C3] = constant(0,0,alphaR,To4R,pi_cHR,pi_fR)
componentR = component;
error = 1;
alpha_i = alphaR;
[state, component,v0] = off_ambient(state,component,design,alt,M0,A0);
[state,component] = off_inlet(state,component);
while error > .0001
[state,component] = off_fan(state,component,componentR,design);
[state,component] = off_comp(state,component,design, componentR);
% [To4] = thetabreak(state,inputs);
% state(9,3) = {To4};
[state,component] = off_burner(state,component,design);
[state,component] = off_turb(state,component);

[pi_cH] = component{6,2};
[tau_r,~,tau_f] = component{2:4,4};
tau_lambda = state{9,8}/state{2,8};
alpha = C3 / (pi_cH*sqrt(tau_r*tau_f/tau_lambda));
error = norm(alpha-alpha_i)/alpha;
alpha_i = alpha
design{2,2} = alpha_i;
end
[state,component,performance] = off_nozzle(state,component,M0,v0,design);
end   

%% Functions
function [C3] = constant(M0,alt,alpha,To4,pi_cH,p_if)
gamma = 1.4;
ef = .9;

[T0, ~, ~, ~] = atmosisa(alt);
tau_lambda = To4/T0 * 1.25; %constant to reflect changing Cp without complicated eqns
tau_r = 1 + .2*(M0)^2;
tau_f = p_if ^ (1/(gamma*ef/(gamma-1)));
C3 = alpha*pi_cH*sqrt(tau_r*tau_f/tau_lambda); %why const?
end

function [A] = corrmass(A)
[T0, ~, ~, rho0] = atmosisa(alt); %obtain standard atmospheric conditions
state(2,3) = {T0};
state(2,4) = {0};
state{2,2} = [];
state{2,8} = [];
[state] = unFAIR3(state,2);
[~,~,T0,~,~,cp0,gamma0,~,~] = state{2,:};
R0 = cp0 - cp0/gamma0;
a0 = sqrt(R0*gamma0*T0); %[m/s]
v0 = M0*a0; %[m/s]

To0 = T0*(1+((M0^2)*((gamma0-1)/2))); %find total temperature using isentropic
state(3,3) = {To0};
state{3,2} = [];
state{3,8} = [];
[state] = unFAIR3(state,3);

P0 = state{2,2};
Po0 = state{3,2};
pi_r = Po0/P0;
component{2,2} = pi_r;

h0 = state{2,8};
ho0 = state{3,8};
tau_r = ho0/h0;
component{2,4} = tau_r;


% mdot_on = ;
% To0_on = ;
% Po0_on = ;

mdot0 = A0 * v0 * rho0;
mdot0 = mdot_on * sqrt(To0_on/To0) * Po0/Po0_on
state{2,5} = mdot0;
[state,design] = mdot(state,design);
end
function [state,design] = mdot(state,design)
% Calculates mass flow throughout engine sections

mdot0 = state{2,5};
alpha = design{2,2};
beta = design{3,2};
ep1 = design{4,2};
ep2 = design{5,2};

% MASS FLOW AND FUEL TO AIR CALCULATIONS
f = state{9,4};
mdot4 = state{9,5};
mdot31 = state{8,5};
if size(f,1) == 1
    mdot_f = mdot4 - mdot31;
else
    mdot_f = 0;
end

f0 = 0; %freestream fuel/air ratio

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
function [state, component,v0] = off_ambient(state,component,design,alt,M0,A0)
[T0, ~, ~, rho0] = atmosisa(alt); %obtain standard atmospheric conditions
state(2,3) = {T0};
state(2,4) = {0};
state{2,2} = [];
state{2,8} = [];
[state] = unFAIR3(state,2);
[~,~,T0,~,~,cp0,gamma0,~,~] = state{2,:};
R0 = cp0 - cp0/gamma0;
a0 = sqrt(R0*gamma0*T0); %[m/s]
v0 = M0*a0; %[m/s]

To0 = T0*(1+((M0^2)*((gamma0-1)/2))); %find total temperature using isentropic
state(3,3) = {To0};
state{3,2} = [];
state{3,8} = [];
[state] = unFAIR3(state,3);

P0 = state{2,2};
Po0 = state{3,2};
pi_r = Po0/P0;
component{2,2} = pi_r;

h0 = state{2,8};
ho0 = state{3,8};
tau_r = ho0/h0;
component{2,4} = tau_r;

mdot0 = A0 * v0 * rho0;
state{2,5} = mdot0;
[state,design] = mdot(state,design);
end
function [state,component] = off_inlet(state,component)
pid = component{3,2};
Po0 = state{3,2};
Po2 = pid*Po0;

state(4,2) =  {Po2};
state{4,3} = [];
state{4,8} = [];
[state] = unFAIR3(state,4);

ho0 = state{3,8};
ho2 = state{4,8};
tau_d = ho2/ho0;
component{3,4} = tau_d;
end
function [state,component] = off_fan(state,component,componentR,design)
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

tau_lambda = ho4/h0;
tau_f = 1 ...
    + ( (1-tau_tL)*eta_mL* (((1-beta-ep1-ep2)*(1+f)*tau_lambda*tau_tH/tau_r  +  (ep1*tau_tH + ep2)*tau_cL*tau_cH))   -   (1+alpha)*PtoL/(tau_r*eta_mPL*mdot0*h0)) ...
    / ((tau_cLR - 1)/(tau_fR - 1) + alpha);


statei = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
statei(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
ho2 = state{4,8};
ho13i = ho2*(1+eta_f*(tau_f-1));
statei{5,8} = ho13i;
statei{5,4} = 0;
[statei] = unFAIR3(statei,5);
Pro13i = statei{5,2};
Pro2 = state{4,2};

pi_f = Pro13i/Pro2;

component{4,2} = pi_f;
component{4,4} = tau_f;

state{5,2} = [];
state{5,3} = [];
ho13 = ho2*tau_f;
state{5,8} = ho13;
[state] = unFAIR3(state,5);
end
function [state,component] = off_comp(state,component,design, componentR)
statei = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
statei(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};

% LP Compressor
tau_f = component{4,4};
tau_fR = componentR{4,4};
tau_cLR = componentR{5,4};
eta_cL = component{5,5};

tau_cL = 1 + (tau_f -1)*(tau_cLR - 1)/(tau_fR - 1);

ho2 = state{4,8};
ho25i = ho2*(1+eta_cL*(tau_cL-1));
statei{6,8} = ho25i;
statei{6,4} = 0;
[statei] = unFAIR3(statei,6);
Pro25i = statei{6,2};
Pro2 = state{4,2};
pi_cL = Pro25i/Pro2;

component{5,2} = pi_cL;
component{5,4} = tau_cL;

state{6,2} = [];
state{6,3} = [];
ho25 = ho2*tau_cL;
state{6,8} = ho25;
[state] = unFAIR3(state,6);



% HP Compressor
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

tau_lambda = ho4/h0;

tau_cH = 1 ...
    + (1-tau_tH)*etamH* (((1-beta-ep1-ep2)*(1+f)*tau_lambda/(tau_r*tau_cL) + ep1*tau_cH))...
    - (1+alpha)*PtoH/(tau_r*tau_cL*etamPH*mdot0*h0);

ho25 = state{6,8};
ho3i = ho25*(1+eta_cH*(tau_cH-1));
statei{7,8} = ho3i;
statei{7,4} = 0;
[statei] = unFAIR3(statei,7);
Pro3i = statei{7,2};
Pro25 = state{6,2};

pi_cH = Pro3i/Pro25;

component{6,2} = pi_cH;
component{6,4} = tau_cH;

state{7,2} = [];
state{7,3} = [];
ho3 = ho25*tau_cL;
state{7,8} = ho3;
[state] = unFAIR3(state,7);
end
function [state,component] = off_burner(state,component,design)
state(8,2:3) = state(7,2:3);
state(8,6:12) = state(7,6:12);
mdot31 = state{8,3};
h_PR = design{8,2};

f_i = state{9,4};
error = 1;
while error > .01
    state{9,2} = [];
    state{9,8} = [];
    [state] = unFAIR3(state,9);

    ho31 = state{7,8};
    ho4 = state{9,8};
    taub = ho4/ho31;
    component{7,4} = taub;

    mdotf = mdot31*(ho4 - ho31) / h_PR; %plus eta_burner
    f = mdotf/mdot31;
    state{9,4} = f;
    [state,design] = mdot(state,design);
    error = (f - f_i)/f_i;
    f_i = f;
end

end
function [state,component] = off_turb(state,component)
% Mixer 1
mdotep1 = state{21,5};
mdot4 = state{9,5};
mdot41 = state{10,5};
hoep1 = state{8,8};
ho4 = state{9,8};
ho41 = (mdotep1*hoep1 + mdot4*ho4) / (mdot41);
state{10,2} = [];
state{10,3} = [];
state(10,8) = {ho41};
[state] = unFAIR3(state,10);
taum1 = ho41/ho4;
component{8,4} = taum1;

% HP Turbine
tau_tH = component{9,4};
state{11,2} = [];
state{11,3} = [];
ho44 = ho41*tau_tH;
state{11,8} = ho44;
[state] = unFAIR3(state,11);

% Mixer 2
mdotep2 = state{22,5};
mdot44 = state{11,5};
mdot45 = state{12,5};
hoep2 = state{8,8};
ho44 = state{11,8};
ho45 = (mdotep2*hoep2 + mdot44*ho44) /(mdot45);
state{12,2} = [];
state{12,3} = [];
state(12,8) = {ho45};
[state] = unFAIR3(state,12);
taum2 = ho45/ho44;
component{11,4} = taum2;

% LP Turbine
tau_tL = component{12,4};
state{13,2} = [];
state{13,3} = [];
ho5 = ho45*tau_tL;
state{13,8} = ho5;
[state] = unFAIR3(state,13);
end
function [state,component,performance] = off_nozzle(state,component,M0,v0,design)
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



%Calculate Momentum
pinlet = mdot0*v0;
pcore = mdot9*(v9);
pbypass = mdot19*(v19);


%Calculate performance parameters
F = pcore + pbypass - pinlet; %Thrust
S = mdotf / F; %Specific fuel consumption
mech_power = -.5*mdot0*v0^2 + .5*mdot19*v19^2 + .5*mdot9*v9^2 + PtoH + PtoL; %Mechanical power consumed
thrust_power = F*v0; %Thrust power generated
chem_power = mdotf*h_PR; %Chemical power released
eta_o = thrust_power/chem_power; %Overall engine efficiency
eta_th = mech_power/chem_power; %Thermal engine efficiency
eta_p = thrust_power/mech_power; %Propulsive engine efficiency

% Alternative equations34
% eta_o = v0/(h_PR*S);
% % eta_o = F*v0 / ((state{9,5} - state{8,5})*h_PR);
% eta_p = (F*v0) / (.5*mdot19*v19^2 + .5*mdot9*v9^2 - .5*mdot0*v0^2 + PtoH + PtoL);
% eta_th = eta_o/eta_p;

% Installation bad stuff BIG APPROX
F = F* (1-M0);

% Store parameters
performance(1,:) = {'Thrust (N)','Specific Fuel Consumption (kg/N-s)','Propulsive Efficiency','Thermal Efficiency','Overall Efficiency','Bypass Exhaust Mach','Core Flow Mach','Inlet Momentum','Core Momentum','Bypass Momentum','Inlet Mass Flow','Core Mass Flow','Bypass Mass Flow','Inlet Velocity','Core Velocity','Bypass Velocity'};
performance(2,:) = {F,S,eta_th,eta_p,eta_o,M19,M9,pinlet,pcore,pbypass,mdot0,mdot9,mdot19,v0,v9,v19};
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