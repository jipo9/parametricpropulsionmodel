function [state,component,design,inputs,performance] = off_design(state,component,design,inputs,componentR,M0,alt,A0)
%% Read Me
% This function takes in the results from the on-design analysis and
% returns the engine performance for any given flight condition

%% Calculations
tauf_i = component{4,4};
error = 1;
while error > .0001
[state, component,v0] = off_ambient(state,component,design,alt,M0,A0);
[state,component] = off_inlet(state,component,inputs);
[state,component] = off_fan(state,component,componentR,design);
[state,component] = off_comp(state,component,design, componentR);
[state,component] = off_burner(state,component,design);
[state,component] = off_turb(state,componentR,component);
[state,component,performance] = off_nozzle(state,component,v0,design);

tauf = component{4,4};
error = norm(tauf-tauf_i)/tauf;
tauf_i = tauf;
end
    
end   

%% Functions
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

mdot6A = mdot45 + mdot13; %mass flow rate after addtion of bypass air
f6A = f45*mdot45 / mdot6A; %fuel/air ratio after addtion of bypass air

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
[state] = unFAIR3(state,2);
[~,~,T0,~,~,cp0,gamma0,~,~] = state{2,:};
R0 = cp0 - cp0/gamma0;
a0 = sqrt(R0*gamma0*T0); %[m/s]
v0 = M0*a0; %[m/s]

mdot0 = A0 * v0 * rho0;
state{2,5} = mdot0;
[state,design] = mdot(state,design);

To0 = T0*(1+((M0^2)*((gamma0-1)/2))); %find total temperature using isentropic
state(3,3) = {To0};
[state] = unFAIR3(state,3);


P0 = state{2,2};
Po0 = state{3,2};
pi_r = Po0/P0;
component{2,2} = pi_r;

h0 = state{2,8};
ho0 = state{3,8};
tau_r = ho0/h0;
component{2,4} = tau_r;



% mdot0 = (A0*Po)*sqrt(gamma0/(To0*R0))*M0*(1+((gamma0-1)/2)*M0^2)^(-1*(gamma0+1)/(2*(gamma0-1))); %compressible eqn
% state{2,5} = mdot0;
end
function [state,component] = off_inlet(state,component,inputs)
M0 = inputs{3,2};
pi_dmax = component{3,2}(1);
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
function [state,component] = off_turb(state,componentR,component)
 tau_m1 = componentR{8,4};
% tau_tH = component{9,4};
 tau_m2 = componentR{11,4};
 ho4 = state{9,8};
% 
% A45_6 = 1;
% M6 = .5;
% 
% gamma = 1.3;
% eta_tL = component{12,6};
% X = A45_6 * 1/M6 * (2/(gamma+1)    *  (1+(gamma-1)/2*M6^2)  )^  ((gamma+1)/(2*(gamma-1)));
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
% 
% 
% 
% %         tau_tL = 0.8598;
% %         pi_tL = 0.4831;
% %         
% %         tau_m1 = 0.9673;
% %         tau_tH = 0.8381;
% %         tau_m2 = 0.9731;
% %         
% %         component{8,4} = tau_m1;
% %         component{9,4} = tau_tH;
% %         component{11,4} = tau_m2;
%         

pi_tH = componentR{9,2};
tau_tH = componentR{9,4};
pi_tL = componentR{12,2};
tau_tL = componentR{12,4};

state{10,2} = [];
state{10,3} = [];
ho41 = ho4*tau_m1;
state{10,8} = ho41;
[state] = unFAIR3(state,10);

state{11,2} = [];
state{11,3} = [];
ho44 = ho41*tau_tH;
state{11,8} = ho44;
[state] = unFAIR3(state,11);

state{12,2} = [];
state{12,3} = [];
ho45 = ho44*tau_m2;
state{12,8} = ho45;
[state] = unFAIR3(state,12);

state{13,2} = [];
state{13,3} = [];
ho5 = ho45*tau_tL;
state{13,8} = ho5;
[state] = unFAIR3(state,13);

%         component{9,2} = (state{11,2} / state{10,2})^(1/.89);
%         component{9,2} = 0.4231;


component{12,2} = pi_tL;
component{12,4} = tau_tL;
component{9,2} = pi_tH;
component{9,4} = tau_tH;
end
function [state,component,performance] = off_nozzle(state,component,v0,design)
%Take in inputs    
alpha = design{2,2};
beta = design{3,2};
P0_P9 = 1;
h_PR = design{8,2};
PtoL = design{6,2};
PtoH = design{7,2};
[~,pi_r,pi_d,~,pi_cL,pi_cH,pi_b,~,pi_tH,~,~,pi_tL,~,~,~,pi_n,~] = component{:,2};

%Calculate pressure drop across nozzle
Pro5 = state{13,2};
pi_n = component{16,2};
Pro9 = Pro5*pi_n;
state(17,2) = {Pro9};
state{17,3} = [];
state{17,8} = [];
[state] = unFAIR3(state,17);

%Calculate static conditions of exhaust
Po9_P0 = pi_r*pi_d(2)*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_n;
Po9_P9 = Po9_P0/P0_P9;
Pr9 = Pro9 / Po9_P9;
state(18,2) = {Pr9};
[state] = unFAIR3(state,18);
component{16,4} = state{17,8} / state{13,8};

%Find exhause mach number
[~,~,~,~,~,~,~,ho9] = state{17,:};
[~,~,T9,~,~,cp9,gamma9,h9] = state{18,:};
R9 = cp9 - cp9/gamma9;
a9 = sqrt(R9*gamma9*T9); %m/s
v9 = sqrt(2*(ho9-h9));
M9 = v9 / a9;

%Find misc values
[~,Pr0,~,~,~,cp0,gamma0,~] = state{2,:};
R0 = cp0 - cp0/gamma0;
mdot0 = state{2,5};
f_0 = state{18,4};

%Calculate Performance
F_mdot = (1+f_0-(beta/(1+alpha)))*v9     -   v0  +   (1+f_0-(beta/(1+alpha)))*R9*T9*(1-Pr0/Pr9)/(R0*v9*gamma0);
S = f_0 / F_mdot;

eta_TH = ((v0^2/2*((1+f_0-(beta/(1+alpha)))*(v9/v0)^2 - 1) + (PtoL + PtoH)/mdot0))/...
    (f_0*h_PR);
eta_P = (2*F_mdot/v0)/...
    ((1+f_0-beta/(1+alpha))*(v9/v0)^2 - 1);
eta_o = eta_TH*eta_P;
performance(1,:) = {'Thrust','Specific Fuel Consumption','Propulsive Efficiency','Thermal Efficiency','Overall Efficiency','Exhaust Mach'};
performance(2,:) = {F_mdot,S,eta_TH,eta_P,eta_o,M9};
end
