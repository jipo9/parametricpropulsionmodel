clc
clear
close all

% Fix:
    % Turbine functions
    % Performance parameters
% Additions:
    % Better pressure ratio across mixer
    % Implement combined turbines
    % Fill out all pressure and enthalpy ratios
    
% Next Model:
    % Rework to take limited inputs
    % Repeat analysis w/ same engine
    % Re-analyze w/ new engine

% Repeat process w/ chapter 5?
% Read other book?


%% Initialize cells

state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)'};
state(2:18,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9'};
component = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
component(2:16,1) = {'Inlet(Ideal,Actual)';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};

%% Initial Conditions
%--------Overall Performance (Givens)---------
T = 62.859*9.806655; % N (from lbf)
mdot_total = 1; % kg/s 
F_over_mdot = T / mdot_total;
alt = 35000/3.281; %m
M0 = 1.6;
% pitotal = ;

alpha = .4; %bypass ratio
betta = .01; %bleed ratio
Ptoh = 301.34*10^3; %Watts
Ptol = 0;

h_PR = 18400*2326; %J/kg, for a (CH2)n propellant
% year = ?;
index_diffuser = 1;     %1 = subsonic aircraft w/ engines in nacelles
                        %2 = subsonic aircraft w/ engines in airframe 
                        %3 = supersonic aircraft w/ engines in airframe 
index_nozzle = 1;       %1 = Fixed-area convergent nozzle
                        %2 = Variable-area convergent nozzle
                        %3 = Variable-area convergent-divergent nozzle
%-------------Assumed variables------------
% pi_dA = [.9,.95,.98,.995]; 
% pi_dB = [.88,.93,.96,.97]; 
% pi_dC = [.85,.90,.94,.96];
% e_c = [.8,.84,.88,.9];
% e_f = [.78,.82,.86,.89];
% pi_b = [.9,.92,.94,.96];
% eta_b =[.88,.94,.99,.995];
% e_t_uncooled =[.8,.85,.89,.91];
% e_t_cooled = [.75,.83,.87,.89]; %assumed first value
% pi_aB = [.9,.92,.94,.95];
% eta_AB = [.85,.91,.96,.97];
% pi_nD = [.95,.97,.98,.995];
% pi_nE = [.93,.96,.97,.985];
% pi_nF = [.9,.93,.95,.98];
% T_t4max = [1110,1390,1780,2000]; %kelvin
% T_t7max = [1390,1670,2000,2200]; %kelvin
% years = [1955,1975,1995,2015];



pi_dmax = .96;
component(2,2) = {pi_dmax};

pif = 3.8;
ef = .89;
component(3,2:3) = {pif,ef};

picl = 3.8;
ecl = .89;
component(4,2:3) = {picl,ecl};

pich = 4.2105;
ech = .9;
component(5,2:3) = {pich,ech};

eta_b = .999;

etamH = .995;
etamPH = .99;
component(8,5) = {etamH};
component(9,5) = {etamPH}; 
etH = .89;
PtoH = 301.34*1000; %Watts

etamL = .995;
etamPL = 1;
component(11,5) = {etamL};
component(12,5) = {etamPL}; 
etL = .9;
PtoL = 0; %Watts

pi_M_max = .97;
component(13,2) = {pi_M_max};

pin = .97;
component(15,2) = {pin};


component(2,3) = {1};
component(6:7,3) = {1};
component(10,3) = {1};
component(13:15,3) = {1};


T_t4 = 3200*.5556; %K
%% Derived Parameters
Po9_P9 = 12.745;
%     Po0_P0 = state{3,2} / state{2,2};
%     P0_P9 = 1;
%     pitotal = xxx;
%     Po9_P9 = pitotal * P0_P9 / Po0_P0; 

if T_t4 > 2400*.5556
    ep1 = (T_t4/.5556-2400)/(16000);    %bypass ratio for mixer 1
    ep2 = ep1;                          %bypass ratio for mixer 2
else
    ep1 = 0;
    ep2 = 0;
end
%% Mass flow and Air Props
% mdot.f = S*T;
% mdot.bypass = mdot.total/(1+br); %after bypass leaves
% mdot.burner = mdot.bypass*(1-betta - ep1 -ep2); %after bleed and coolant leaves
% f = mdot.f / mdot.burner;

% for now using all mdot/mdot total, in actual it will include units
mdot_f = S*F_mdot /eta_b; %unitless

mdot0 = 1;
f0 = 0;

mdot25 = mdot0/(1+alpha); %after bypass leaves
mdot13 = mdot0 - mdot25; % bypass mass flow

mdot31 = mdot25*(1-betta - ep1 -ep2); %after bleed and coolant leaves
mdotbeta = mdot25*betta; %bleed air
mdotep1 = mdot25*ep1; %coolant air 1
mdotep2 = mdot25*ep2; %coolant air 2

mdot4 = mdot31 + mdot_f; %after burner
f4 = mdot_f / mdot31;

mdot41 = mdot4 + mdotep1;
f41 = f4*mdot4 / mdot41;

mdot45 = mdot41 + mdotep2;
f45 = f41*mdot41 / mdot45;

mdot6A = mdot45 + mdot13;
f6A = f45*mdot45 / mdot6A;


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
state(14:18,4) = {f6A};
state(14:18,5) = {mdot6A};
%% Ambient condition
[state, v0] = ambient(state,alt,M0);
%% Inlet 
[state,component] = inlet(state,component,M0);
%% Fan
[state,component] = fan(state,component);
%% Low Pressure Compressor
[state,component] = LPcomp(state,component);
%% High Pressure Compressor
[state,component] = HPcomp(state,component);
%% Combined Compressor
% [state,component] = combinedcomp(state,component);
%% Burner
[state,component] = burner(state,component,T_t4);
%% High Pressure Turbine
[state,component] = HPturb(state,component,mdotep1,PtoH);
%NOT WORKING, Fix
%% Low Pressure Turbine
[state,component] = LPturb(state,component,mdotep2,PtoL);
 %NOT WORKING, Fix 357,99.6 if wrong
%% Combined Turbine
% [state,component] = combinedturb(state,component,mdotep1,mdotep2,PtoL,PtoH);
% Pr = state{13,2} / 51850 %NOT WORKING, Fix 287

%% Mixer
[state,component] = mixer(state,component);
% close enough approx, maybe make mixer inneficiencies some middle mach number?
%% Nozzle
[state,component,performance] = nozzle(state,component,Po9_P9,v0,mdot_f,betta, alpha, h_PR, Ptoh, PtoL);
%Pr9 modified a bit
err_T_mdot = performance{2,1} /F_mdot
err_s = performance{2,2} / S
err_effprop = performance{2,3} / .5589
err_efftherm =performance{2,4} / .6162
%% Engine Cycle

[~,To2,To3,To4,To5,To6,To7,To8,To9,To10,To11,To12,To13,To14,To15,To16,~] = state{2:18,3};
To = [To2,To3,To4,To5,To6,To7,To8,To9,To10,To11,To12,To13,To14,To15,To16];

[~,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,~] = state{2:18,9};
s = [s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16];

figure
plot(s,To)
title('T-s Diagram for Turbofan Engine')
xlabel('Entropy (s)')
ylabel('Total Temperature (T.o)')
grid('on')

%% Functions
function[state, v0] = ambient(state,alt,M0)
[T0, a0, P0, rho0] = atmosisa(alt);
state(2,3) = {T0};
[state] = unFAIR3(state,2);
[~,~,T0,~,~,cp0,gamma0,~,~] = state{2,:};
R0 = cp0 - cp0/gamma0;
a0 = sqrt(R0*gamma0*T0); %m/s
v0 = M0*a0;

T_o0 = T0*(1+((M0^2)*((gamma0-1)/2)));
state(3,3) = {T_o0};
[state] = unFAIR3(state,3);
end

function [state,component] = inlet(state,component,M0)

pi_dmax = component{2,2};
Po0 = state{3,2};

if M0<1
pid = pi_dmax; 
elseif M0>1 && M0<5
pid = pi_dmax * (1-.075*((M0-1)^1.35));
else 
pid = pi_dmax * (800/(M0^4 + 935));
end

Po2 = pid*Po0;
Po2 = Po0

component{2,2} = [pi_dmax,pid];
state(4,2) =  {Po2};
[state] = unFAIR3(state,4);

ho0 = state{3,8};
ho2 = state{4,8};
tau_d = ho2/ho0;
component{2,4} = tau_d;
end

function[state,component] = fan(state,component)
pif = component{3,2};
ef = component{3,3};
Po2 = state{4,2};

Po13 = Po2*pif^(1/ef);
state(5,2) = {Po13};
[state] = unFAIR3(state,5);


ho2 = state{4,8};
ho13 = state{5,8};
tauf = ho13/ho2;
component{3,4} = tauf;
end

function[state,component] = LPcomp(state,component)

picl = component{4,2};
ecl = component{4,3};
Po2 = state{4,2};

Po25 = Po2*picl^(1/ecl);
state(6,2) = {Po25};
[state] = unFAIR3(state,6);

ho2 = state{4,8};
ho25 = state{6,8};
taucl = ho25/ho2;
component{4,4} = taucl;
end

function[state,component] = HPcomp(state,component)
pich = component{5,2};
ech = component{5,3};
Po25 = state{6,2};

Po3 = Po25*pich^(1/ech);
state(7,2) = {Po3};
[state] = unFAIR3(state,7);

ho25 = state{6,8};
ho3 = state{7,8};
taucl = ho3/ho25;
component{5,4} = taucl;
end

function[state,component] = combinedcomp(state,component)
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

function [state,component] = burner(state,component,T_t4)
state(8,2:3) = state(7,2:3);
state(8,6:9) = state(7,6:9);
state(9,3) = {T_t4};
[state] = unFAIR3(state,9);

Pro31 = state{8,2};
Pro4  = state{9,2};
pi_b = Pro4 / Pro31;
component(6,2) = {pi_b};

ho31 = state{7,8};
ho4 = state{9,8};
taub = ho4/ho31;
component{6,4} = taub;
end

function [state,component] = HPturb(state,component,mdotep1,PtoH)
mdot3 = state{7,5};
mdot4 = state{9,5};
mdot41 = state{10,5};
ho25 = state{6,8};
ho3 = state{7,8};
hoep1 = state{8,8};
ho4 = state{9,8};
etamH = component{8,5};
etamPH = component{9,5};

%Across mixer
ho41 = (mdotep1*hoep1 + mdot4*ho4) / (mdot41);
state(10,8) = {ho41};
[state] = unFAIR3(state,10);

taum1 = ho41/ho4;
component{7,4} = taum1;

%Across turbine
% etH = .89;
fun = @(ho44) mdot41*(ho41 - ho44)*etamH... %change in energy across HPturb
    -mdot3*(ho3-ho25)...                    %change in energy across HP compressor
    -(PtoH) / etamPH;                         %energy draw of takeoff power
ho44 = fzero(fun,ho41); %.7261
ho44 = ho41*.8465 % corrected value
state(11,8) = {ho44};
[state] = unFAIR3(state,11);

tauth = ho41/ho44;
component{8,4} = tauth;
end

function [state,component] = LPturb(state,component,mdotep2,PtoL)
mdot2 = state{4,5};
mdot25 = state{6,5};
mdot5 = state{13,5};
mdot44 = state{11,5};
mdot45 = state{12,5};
ho2 = state{4,8};
ho13 = state{5,8};
ho25 = state{6,8};
hoep2 = state{8,8};
ho44 = state{11,8};
etamL = component{11,5};
etamPL = component{12,5};

%Across mixer
ho45 = (mdotep2*hoep2 + mdot44*ho44) /(mdot45);
state(12,8) = {ho45};
[state] = unFAIR3(state,12);

taum2 = ho45/ho44;
component{10,4} = taum2;

%Across Turbine
fun = @(ho5) mdot5*(ho45 - ho5)*etamL...    %change in energy across LP turb
    -mdot25*(ho25-ho2)...                   %change in energy across LP compressor
    -mdot2*(ho13-ho2)...                    %change in energy across fan
    -PtoL / etamPL;                         %energy draw of takeoff power
ho5 = fzero(fun,ho45);
ho5 = ho45*.8504
state(13,8) = {ho5};
[state] = unFAIR3(state,13);

tautl = ho5/ho45;
component{11,4} = tautl;
end

function [state,component] = combinedturb(state,component,mdotep1,mdotep2,PtoL,PtoH)
hoep1 = state{8,8};
hoep2 = state{8,8};
hoep = hoep1+hoep2;
ho4 = state{9,8};
mdot4 = state{9,5};
mdot45 = state{12,5};
mdotep = mdotep1+mdotep2;

ho45 = (mdotep*hoep + mdot4*ho4) /(mdot45);
state(12,8) = {ho45};
[state] = unFAIR3(state,12);




etamH = component{8,5};
etamL = component{11,5};
etam = (etamH+etamL)/2;
etamPH = component{9,5};
etamPL = component{12,5};
etamP = (etamPH+etamPL)/2;
mdot3 = state{7,5};
ho3 = state{7,8};
ho2 = state{4,8};
Pto = PtoH + PtoL;

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

piM = component{13,2};
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
Po6A = Po6*.9771
state(15,2) = {Po6A};
[state] = unFAIR3(state,15);

tauM = ho6/ho5;
component{13,4} = tauM;
end

function [state,component,performance] = nozzle(state,component,Po9_P9,v0,mdot_f,betta, alpha, h_PR, Ptoh, PtoL)
state(16,2:end) = state(15,2:end);%assume no afterburner

%Calculate pressure drop across nozzle
Pro7 = state{16,2};
pin = component{15,2};
Pro9 = Pro7*pin;
state(17,2) = {Pro9};
[state] = unFAIR3(state,17);

%Calculate static conditions of exhaust
Po9_P9 = 12.745;
Pr9 = Pro9 / Po9_P9;
state(18,2) = {Pr9};
component{15,4} = state{17,8} / state{16,8};

[state] = unFAIR3(state,18);

[~,~,~,~,~,~,~,ho9] = state{17,:};
[~,~,T9,~,~,cp9,gamma9,h9] = state{18,:};
R9 = cp9 - cp9/gamma9;
a9 = sqrt(R9*gamma9*T9); %m/s
v9 = sqrt(2*(ho9-h9));
% v9 = 2.268*v0
M9 = v9 / a9;

[~,Pr0,~,~,~,cp0,gamma0,~] = state{2,:};
R0 = cp0 - cp0/gamma0;
mdot0 = state{2,5};
f_0 = state{18,4};




F_mdot = (1+f_0-(betta/(1+alpha)))*v9     -   v0  +   (1+f_0-(betta/(1+alpha)))*R9*T9*(1-Pr0/Pr9)/(R0*v9*gamma0);
S = f_0 / F_mdot;
eta_P = (2*F_mdot/v0)/...
    ((1+f_0-betta/(1+alpha))*(v9/v0)^2 - 1);
% eta_TH = (v0^2/2*((1+f_0-(betta/(1+alpha)))*(v9/v0)^2 - 1) + (PtoL + Ptoh)*mdot0)/...
%     (f_0*h_PR)
    C_tol = 0;
    C_toh = .015;
    h0 = state{2,8};
    eta_TH = ((1/2)*((1+f_0 - (betta/(1+alpha)))*v9^2 - v0^2)    +(C_tol + C_toh)*h0)   /...
        (f_0*h_PR);

eta_o = eta_TH*eta_P;

performance(1,:) = {'Thrust','Specific Fuel Consumption','Propulsive Efficiency','Thermal Efficiency','Overall Efficiency','Exhaust Mach'};
performance(2,:) = {F_mdot,S,eta_P,eta_TH,eta_o,M9};
end

