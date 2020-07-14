function[state, component, performance, design,error] = Turbofan_Barebones_iteration(state, design, component, M0, alt, mdotep,F_mdot,S)
%% Ambient condition
[state, component,v0] = ambient(state,component,alt,M0);
%% Inlet 
[state,component] = inlet(state,component,M0);
%% Fan
[state,component] = fan(state,component);
%% Low Pressure Compressor
% [state,component] = LPcomp(state,component);
%% High Pressure Compressor
% [state,component] = HPcomp(state,component);
%% Combined Compressor
[state,component] = combinedcomp(state,component);
%% Burner
[state,component] = burner(state,component);
%% High Pressure Turbine
% [state,component] = HPturb(state,component,design,mdotep1);
%% Low Pressure Turbine
% [state,component] = LPturb(state,component,design,mdotep2);
%% Combined Turbine
[state,component] = combinedturb(state,component,design,mdotep);
%% Mixer
[state,component] = mixer(state,component);
% close enough approx, maybe make mixer inneficiencies some middle mach number?
%% Nozzle
[state,component,performance] = nozzle(state,component,v0,design);
%% Results

F_mdot= 637.533110510502;
S = 3.17869785560938e-05;
efftherm = 0.577305593707712;
effprop = 0.608344703158938;

err_T_mdot = performance{2,1} /F_mdot; %T/mdot error compared to book
err_s = performance{2,2} / S; %S error compared to book
err_efftherm = performance{2,3} / efftherm; %thermal efficiency error compared to book
err_effprop =performance{2,4} / effprop; %propulsive efficiency compared to book
error = [err_T_mdot,err_s,err_efftherm,err_effprop];

%% Engine Cycle

% [~,To2,To3,To4,To5,To6,To7,To8,To9,To10,To11,To12,To13,To14,To15,To16,~] = state{2:18,3};
% To = [To2,To3,To4,To5,To6,To7,To8,To9,To10,To11,To12,To13,To14,To15,To16];
% 
% [~,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,~] = state{2:18,9};
% s = [s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16];
% 
% figure
% plot(s,To)
% title('T-s Diagram for Turbofan Engine')
% xlabel('Entropy (s)')
% ylabel('Total Temperature (T.o)')
% grid('on')

%% Functions
function [state, component,v0] = ambient(state,component,alt,M0)
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

function [state,component] = inlet(state,component,M0)

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
state(8,6:9) = state(7,6:9);
[state] = unFAIR3(state,9);
%Add in pressure losses
% pi_b = component{7,2};
% Pro4_ideal  = state{9,2};
% Pro4_actual = 1*pi_b*Pro4_ideal;
% state(9,2:3) = {Pro4_actual,[]};
% state(9,8) = {[]};
% [state] = unFAIR3(state,9);

Pro31 = state{8,2};
Pro4  = state{9,2};
pi_b_total = Pro4 / Pro31;
component(7,2) = { pi_b_total};

ho31 = state{7,8};
ho4 = state{9,8};
taub = ho4/ho31;
component{7,4} = taub;
end

function [state,component] = HPturb(state,component,design,mdotep1)
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

%Across mixer
ho41 = (mdotep1*hoep1 + mdot4*ho4) / (mdot41);
state(10,8) = {ho41};
[state] = unFAIR3(state,10);

taum1 = ho41/ho4;
component{8,4} = taum1;

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
end

function [state,component] = LPturb(state,component,design,mdotep2)
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

%Across mixer
ho45 = (mdotep2*hoep2 + mdot44*ho44) /(mdot45);
state(12,8) = {ho45};
[state] = unFAIR3(state,12);

taum2 = ho45/ho44;
component{11,4} = taum2;

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
end

function [state,component] = combinedturb(state,component,design,mdotep)
hoep = state{8,8};
ho4 = state{9,8};
mdot4 = state{9,5};
mdot45 = state{12,5};

ho45 = (mdotep*hoep + mdot4*ho4) /(mdot45);
state(12,8) = {ho45};
[state] = unFAIR3(state,12);


etamL = component{12,5};
etamPL = component{13,5};
Pto = design{6,2};
etamH = component{9,5};
etamPH = component{10,5};


etam = (etamH+etamL)/2;
etamP = (etamPH+etamPL)/2;

mdot3 = state{7,5};
ho3 = state{7,8};
ho2 = state{4,8};

fun = @(h5) mdot45*(ho45 - h5)*etam...  %change in energy across turb
    -mdot3*(ho3-ho2)...                 %change in energy across HP compressor
    -Pto / etamP;                       %energy draw of takeoff power
ho5 = fzero(fun,ho45);
state(13,8) = {ho5};
[state] = unFAIR3(state,13);
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

function [state,component,performance] = nozzle(state,component,v0,design)
alpha = design{2,2};
beta = design{3,2};
Pto = design{6,2};
pitotal = design{7,2};
h_PR = design{8,2};
state(16,2:end) = state(15,2:end);%assume no afterburner

%Calculate pressure drop across nozzle
Pro7 = state{16,2};
pin = component{16,2};
Pro9 = Pro7*pin;
state(17,2) = {Pro9};
[state] = unFAIR3(state,17);

%Calculate static conditions of exhaust
[~,Pr0,~,~,~,cp0,gamma0,~] = state{2,:};
R0 = cp0 - cp0/gamma0;
mdot0 = state{2,5};
f_0 = state{18,4};

Po0_P0 = state{3,2} / state{2,2};
P0_P9 = Pr0/ state{2,2};
Po9_P9 = pitotal * P0_P9 / Po0_P0; 
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



F_mdot_ = (1+f_0-(beta/(1+alpha)))*v9     -   v0  +   (1+f_0-(beta/(1+alpha)))*R9*T9*(1-Pr0/Pr9)/(R0*v9*gamma0);
S_ = f_0 / F_mdot_;
eta_TH = ((v0^2/2*((1+f_0-(beta/(1+alpha)))*(v9/v0)^2 - 1) + (Pto)/mdot0))/...
    (f_0*h_PR);
%eta_P = 2/(1+v9/v0); Simplified case, neglected for now
eta_P = (2*F_mdot_/v0)/...
    ((1+f_0-beta/(1+alpha))*(v9/v0)^2 - 1);
eta_o = eta_TH*eta_P;
% eta_o = v0/(h_PR*S); Alternate eqn
performance(1,:) = {'Thrust','Specific Fuel Consumption','Propulsive Efficiency','Thermal Efficiency','Overall Efficiency','Exhaust Mach'};
performance(2,:) = {F_mdot_,S_,eta_TH,eta_P,eta_o,M9};
end
end

