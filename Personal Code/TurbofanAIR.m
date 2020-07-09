clc
clear
clf

% Fix 3 remaining issue spots (turbine 1, turbine 2, nozzle)
% Play around w/ reduced components (collapse turbines, collapse compressors) and see accuracies
% FIll out all info on all three cells
% Try our hand at a real engine (non-ideal case)!!
    % assumed e and pi values will soon dissapear into "levels of tech" and "design decisions"
    % Compare F_mdot and S to actual values
% Repeat process w/ chapter 5?



%% Initial Conditions

%---------Air Conditions-----------
% alt = 35000;
% altm = alt*3.281;
% [T0, a0, P0, rho0] = atmosisa(altm);

state = {'Station','Pressure (Pa)', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg-K)'};
state(2:18,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9'};
component = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
component(2:16,1) = {'Inlet(Ideal,Actual)';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};




%--------Overall Performance (Givens)---------
% year = ?;
% T = ; % N (from lbf)
% mdot.total = ; % kg/s 
% F_over_mdot. = T / mdot.total;
F_mdot = 62.859*9.806655; %N/kg/s from lbf/(lbf/s)
S = 1.1386*((.453592/3600)/4.44822); %kg/s/N from lbm/(lbf/s)
% pitotal = ;

alpha = .4; %bypass ratio
betta = .01; %bleed ratio
% T4 = 3200* .5556; %turbine inlet temp, deg Rankine!!!
h_PR = 18400*2326; %J/kg, for a (CH2)n propellant

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
component(14,2) = {pin};

% M4 = 1; %assume chocked flow at turbine
% T_t4 = T4*(1+((M0^2)*((gamma-1)/2)));
T_t4 = 3200;
if T_t4 > 2400
    ep1 = (T_t4-2400)/16000; %bypass ratio for mixer 1
    ep2 = ep1; %bypass ratio for mixer 2
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

    mdot_s0 = 1;
    f_s0 = 0;
    
    mdot_s25 = mdot_s0/(1+alpha); %after bypass leaves
    mdot_s13 = mdot_s0 - mdot_s25; % bypass mass flow
    
    mdot_s31 = mdot_s25*(1-betta - ep1 -ep2); %after bleed and coolant leaves
    mdotbeta = mdot_s25*betta; %bleed air
    mdotep1 = mdot_s25*ep1; %coolant air 1
    mdotep2 = mdot_s25*ep2; %coolant air 2
    
    mdot_s4 = mdot_s31 + mdot_f; %after burner
    f_s4 = mdot_f / mdot_s31;
    
    mdot_s41 = mdot_s4 + mdotep1;
    f_s41 = f_s4*mdot_s4 / mdot_s41;
    
    mdot_s45 = mdot_s41 + mdotep2;
    f_s45 = f_s41*mdot_s41 / mdot_s45;  
    
    mdot_s6 = mdot_s45 + mdot_s13;
    f_s6 = f_s45*mdot_s45 / mdot_s6;


    state(2:8,4) = {f_s0};
    state(2:4,5) = {mdot_s0};
    state(5,5) = {mdot_s13};
    state(6:7,5) = {mdot_s25};
    state(8,5) = {mdot_s31};
    state(9,4) = {f_s4};
    state(9,5) = {mdot_s4};
    state(10:11,4) = {f_s41};
    state(10:11,5) = {mdot_s41};
    state(12:13,4) = {f_s45};
    state(12:13,5) = {mdot_s45};
    state(13:18,4) = {f_s6};
    state(14:18,5) = {mdot_s6};
    
    
    
    
    M0 = 1.6;
T0 = 394.1 * .5556; %deg Kelvin (from Rankine)
state(2,3) = {T0};
[state] = unFAIR3(state,2);
% P0 = 3.467 * 6895; %Pascal (from psia)
% M0 = .6;
cp0 = state{2,6};
gamma0 = state{2,7};
R0 = cp0 - cp0/gamma0;
a0 = sqrt(R0*gamma0*T0); %m/s
v0 = M0*a0;

T_o0 = T0*(1+((M0^2)*((gamma0-1)/2)));
state(3,3) = {T_o0};
[state] = unFAIR3(state,3);


% P.o0 = P0*(1+((M0^2)*((gamma-1)/2)))^(gamma/(gamma-1));
% rhoo0 = rho0*(1+((M0^2)*((gamma-1)/2)))^(1/(gamma-1));

    Po9_P9 = 12.745;
    Po0_P0 = state{3,2} / state{2,2};
    P0_P9 = 1;
    pitotal = Po9_P9 / P0_P9 / Po0_P0; 
%% Inlet 
[state,component] = inlet(state,component,M0);
Pr = state{4,2} / 51850
%% Fan
[state,component] = fan(state,component);
Pr = state{5,2} / 51850
%% Low Pressure Compressor
[state,component] = LPcomp(state,component);
Pr = state{6,2} / 51850
%% High Pressure Compressor
[state,component] = HPcomp(state,component);
Pr = state{7,2} / 51850
%% Combined Compressor
%% Burner
[state,component] = burner(state,component,T_t4);
Pr = state{9,2}/ 51850
%% High Pressure Turbine
[state,component] = HPturb(state,component,mdotep1,PtoH);
Pr = state{10,2} / 51850
Pr = state{11,2} / 51850 %NOT WORKING, Fix
%% Low Pressure Turbine
[state,component] = LPturb(state,component,mdotep2,PtoL);
Pr = state{12,2} / 51850
Pr = state{13,2} / 51850 %NOT WORKING, Fix
%% Combined Turbine
%% Mixer
[state,component] = mixer(state,component);
Pr = state{14,2} / 51850
Pr = state{15,2} / 51850
% close enough approx, maybe make mixer inneficiencies some middle mach number?
%% Nozzle
% Output not working at all!
[state,component,performance] = nozzle(state,component,Po9_P9,v0,mdot_f,betta, alpha, h_PR, PtoH, PtoL)
err_T_mdot = 1-performance{1,1} /F_mdot
err_s = 1-performance{1,2} / S
err_eff1 = 1-performance{1,3} / .5589
err_eff2 =1-performance{1,4} / .6162
%% Entropy Calcs
% cp = (gamma*R)/(gamma-1);
% 
% ds21 = cp*log(T.o3/T.o0) - R*log(P.o3/P.o0);
% dt21 = T.o3-T.o0;
% 
% ds32 = cp*log(T.o4/T.o3) - R*log(P.o4/P.o3);
% dt32 = T.o4-T.o3;
% 
% ds43 = cp*log(T.o9/T.o4) - R*log(P.o9/P.o4);
% dt43 = T.o9-T.o4;
% 
% h1 = figure(1);
% h1.WindowStyle = 'docked';
% scatter(0,0,'filled')
% hold on
% plot(linspace(0,ds21,10),linspace(0,dt21,10),'linewidth',2,'color','blue')
% hold on
% 
% scatter(ds21,dt21,'filled')
% hold on
% plot(linspace(ds21,ds21+ds32,10),linspace(dt21,dt21+dt32,10),'linewidth',2,'color','blue')
% hold on
% 
% scatter(ds21+ds32,dt21+dt32,'filled')
% hold on
% plot(linspace(ds21+ds32,ds21+ds32+ds43,10),linspace(dt21+dt32,dt21+dt32+dt43,10),'linewidth',2,'color','blue')
% hold on
% 
% scatter(ds21+ds32+ds43,dt21+dt32+dt43,'filled')
% hold on
% plot(linspace(ds21+ds32+ds43,0,10),linspace(dt21+dt32+dt43,0,10),'--','linewidth',2,'color','blue')
% title('T-s Diagram for Turbofan Engine')
% xlabel('Entropy (s)')
% ylabel('T.otal Temperature (T.o)')
% grid('on')

%% Functions
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

%Po2 = pid*Po0;
Po2 = Po0;
component{2,2} = [pi_dmax,pid];
state(4,2) =  {Po2};
[state] = unFAIR3(state,4);

%T.o2 = T.o0*(exp(((gamma-1)/(gamma))*(log(pid))));
%ho0 = ((gamma*R*T.o0)/(gamma-1)) + ((V0^2)/2);
%fun = @(V2) ((gamma*R*T.o2)/(gamma-1)) + ((V2^2)/2) - ho0;
%V2 = fzero(fun,V0);
end

function[state,component] = fan(state,component)
pif = component{3,2};
ef = component{3,3};
Po2 = state{4,2};

Po13 = Po2*pif^(1/ef);

state(5,2) = {Po13};
[state] = unFAIR3(state,5);
% ho2 = ((gamma*R*T.o2)/(gamma-1)) + ((V2^2)/2);
end

function[state,component] = LPcomp(state,component)

picl = component{4,2};
ecl = component{4,3};
Po2 = state{4,2};

Po25 = Po2*picl^(1/ecl);

state(6,2) = {Po25};
[state] = unFAIR3(state,6);


% P.o25 = P.o2*picl^(1/ecl);
% [T.o25] = CPG(T.o2,P.o2,gamma, P.o25, 1);
% % fun = @(T.o25) (T.o25/T.o2)^((gamma*ecl)/(gamma-1)) - picl;
% T.o25 = fzero(fun,T.o2)
%T.o25 = T.o2*(exp(((gamma-1)/(gamma*ecl))*(log(picl))));
end

function[state,component] = HPcomp(state,component)

pich = component{5,2};
ech = component{5,3};
Po25 = state{6,2};

Po3 = Po25*pich^(1/ech);

state(7,2) = {Po3};
[state] = unFAIR3(state,7);


% P.o3 = P.o25*pich^(1/ech);
% [T.o3] = CPG(T.o25,P.o25,gamma, P.o3, 1);
% fun = @(T.o3) (T.o3/T.o25)^((gamma*ech)/(gamma-1)) - pich;
% T.o3 = fzero(fun,T.o25);
%T.o3 = T.o25*(exp(((gamma-1)/(gamma*ech))*(log(pich))));
end

function [state,component] = burner(state,component,T_t4)

state(8,2:3) = state(7,2:3);
state(8,6:8) = state(7,6:8);
state(9,3) = {T_t4*.5556};
[state] = unFAIR3(state,9);

pi_b = state{9,2} / state{8,2};
component(6,2) = {pi_b};
% fun = @(T.o4) (((mdot_4*cp*T.o4)-(mdot_31*cp*T.o3))/(mdot_f*hPR)) - etab;
% T.o4 = fzero(fun,T.o3); %gives Pr of 41 million, not good at all
% T.o4 = 3200*.555556; %gives Pr of 700 (we want 1707)
% T.o4 =  3200*.555556;
% gamma = 1.28
% fun = @(pib) (T.o4/T.o3)^((gamma)/(gamma-1)) - pib;
% pib = fzero(fun,T.o4);
% P.o4 = pib*P.o3;
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
ho41 = (mdotep1*hoep1 + mdot4*ho4) /(mdot41);
state(10,8) = {ho41};
[state] = unFAIR3(state,10);

%Across turbine
% etH = .89;
fun = @(ho44) mdot41*(ho41 - ho44)*etamH... %change in energy across HPturb
    -mdot3*(ho3-ho25)...                    %change in energy across HP compressor
    -PtoH / etamPH;                         %energy draw of takeoff power
ho44 = fzero(fun,ho41); %.7261
ho44 = ho41*.8465;
state(11,8) = {ho44};
[state] = unFAIR3(state,11);
end

function [state,component] = LPturb(state,component,mdotep2,PtoL)
%mdot3 = state{7,5};
mdot2 = state{4,5};
mdot25 = state{6,5};
mdot5 = state{13,5};
mdot44 = state{11,5};
mdot45 = state{12,5};
%ho25 = state{6,8};
%ho3 = state{7,8};
ho2 = state{4,8};
ho13 = state{5,8};
ho25 = state{6,8};
hoep2 = state{8,8};
ho44 = state{11,8};
etamL = component{11,5};
etamPL = component{12,5};

%Across mixer
ho45 = (mdotep2*hoep2 + mdot44*ho44) /(mdot45);
% ho45 = ho44*;
state(12,8) = {ho45};
[state] = unFAIR3(state,12);

%Across Turbine
fun = @(ho5) mdot5*(ho45 - ho5)*etamL...    %change in energy across LP turb
    -mdot25*(ho25-ho2)...                   %change in energy across LP compressor
    -mdot2*(ho13-ho2)...                    %change in energy across fan
    -PtoL / etamPL;                         %energy draw of takeoff power
ho5 = fzero(fun,ho45);
ho5 = ho45*.8504;
state(13,8) = {ho5};
[state] = unFAIR3(state,13);





%do a simple turbine


% etl = .89; %LOT3
% pitl = .4; %ASSUMPTION
% 
% P.o5 = pitl*P.o44;
% 
% %fun = @(T.o5) (T.o5/T.o44)^((gamma*etl)/(gamma-1)) - pitl;
% T.o5 = fzero(fun,T.o44);
% 
% T.o5 = T.o44*(exp(((gamma-1)/(gamma*etl))*(log(pitl))));

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

ho6 = (mdotalpha*hoalpha + mdot5*ho5) /(mdot6);
state(14,8) = {ho6};
[state] = unFAIR3(state,14);

Po6 = state{14,2};
Po6A = Po6*piM;
state(15,2) = {Po6A};
[state] = unFAIR3(state,15);
end

function [state,component,performance] = nozzle(state,component,Po9_P9,v0,mdot_f,betta, alpha, h_PR, P_toH, P_toL)
state(16,:) = state(15,:);%assume no afterburner
pin = component{14,2};
P6A = state{14,2};

Po9 = pin*P6A;
state(17,2) = {Po9};
[state] = unFAIR3(state,17);

%static conditions
P9 = Po9 / Po9_P9;
state(18,2) = {P9};
[state] = unFAIR3(state,18);

% T0 = state{3,2};
T9 = state{18,2};
% To9 = T9;
% T9_T0 = T9/T0;
cp9 = state{18,6};
gamma9 = state{18,7};
R9 = cp9 - cp9/gamma9;
a9 = sqrt(R9*gamma9*T9); %m/s
rho9 = P9 / (R9*T9);
M9 = sqrt(  ((P9/Po9)^((1-gamma9)/gamma9)   -1) /   ((gamma9 - 1)/2));
% fun =  @(M9) Po9 - P9*(1+((M9^2)*((gamma9-1)/2)))^(gamma9/(gamma9-1));
% fun = @(M9) P9/Po9 - (1+(gamma9-1)/2*M9^2)^(-gamma9/(gamma9-1));
% M9 = fzero(fun,Po9);
v9 = M9*a9;
mdot9 = state{17,5};

P0 = state{2,2};
mdot0 = state{2,5};
f_0 = mdot_f / mdot0;




F_mdot = (mdot9*v9 - mdot0*v0)/mdot0 + rho9*v9*(P9-P0);
S = f_0 / F_mdot;
eta_P = (2*F_mdot/v0)/...
    ((1+f_0-betta/(1+alpha))*(v9/v0)^2 - 1);

eta_TH = (v0^2/2*((1+f_0-(betta/(1+alpha)))*(v9/v0)^2 - 1) + (P_toL + P_toH)*mdot0)/...
    (f_0*h_PR);
eta_o = eta_TH*eta_P;

performance(1,:) = {F_mdot,S,eta_P,eta_TH,eta_o};
end

