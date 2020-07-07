%% Turbofan Analysis using ONLY Air

clc
clear
clf

% need to add in power takeoff
% everything here is cpg
% why is the reference pressure such an odd number??
% assumed e and pi values will soon dissapear into "levels of tech" and
%"design decisions"
% Redo graph of accuracy with epsilons
% have JP help w/ function handles
% fix nomenclature!
% burner section 100% needs to be more than CPG if we want to even be on
% the right order of magnitude


%% Constant Setup
R = 287; %J/kg*K
gamma = 1.4;

%% Initial Conditions

%---------Air Conditions-----------
% alt = 35000;
% altm = alt*3.281;
% [T0, a0, P0, rho0] = atmosisa(altm);

M0 = 1.6;
T0 = 394.1 * .5556; %deg Kelvin (from Rankine)
P0 = 3.467 * 6895; %Pascal (from psia)
a0 = sqrt(R*gamma*T0); %m/s
% M0 = .6;
V0 = M0*a0;
T.o0 = T0*(1+((M0^2)*((gamma-1)/2)));
P.o0 = P0*(1+((M0^2)*((gamma-1)/2)))^(gamma/(gamma-1));
% rhoo0 = rho0*(1+((M0^2)*((gamma-1)/2)))^(1/(gamma-1));

%--------Overall Performance (Givens)---------
% year = ?;
% T = ; % N (from lbf)
% mdot.total = ; % kg/s 
% F_over_mdot. = T / mdot.total;
F_mdot = 62.859*9.806655; %N/kg/s from lbf/(lbf/s)
S = 1.1386*((.453592/3600)/4.44822); %kg/s/N from lbm/(lbf/s)
% pitotal = ;
    P.o9_P9 = 12.745;
    P.o0_P0 = P.o0 / P0;
    P0_P9 = 1;
    pitotal = P.o9_P9 / P0_P9 / P.o0_P0; 
alpha = .4; %bypass ratio
beta = .01; %bleed ratio
% T4 = 3200* .5556; %turbine inlet temp, deg Rankine!!!
hPR = 18400*2326; %J/kg, for a (CH2)n propellant

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



mdot.total = 1; %just for now
pi_dmax = .96;
ef = .89;
pif = 3.8;
ecl = .89;
picl = 3.8;
pich = 4.2105;
ech = .9;
etab = .999;
etamH = .995;
etamPH = .99;
etH = .89;

PtoH = 301.34*1000; %Watts

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
    mdot.f = S*F_mdot; %unitless

    mdot.s0 = 1;
    f_s0 = 0;
%     [cp_s0,gamma_s0] = unFAIR2(f_s0);
    
    mdot.s25 = mdot.s0/(1+alpha); %after bypass leaves
    mdot.s13 = mdot.s0 - mdot.s25; % bypass mass flow
    
    mdot.s31 = mdot.s25*(1-beta - ep1 -ep2); %after bleed and coolant leaves
    mdot.sbeta = mdot.s25*beta; %bleed air
    mdot.se1 = mdot.s25*ep1; %coolant air 1
    mdot.se2 = mdot.s25*ep2; %coolant air 2
    
    mdot.s4 = mdot.s31 + mdot.f; %after burner
    f.s4 = mdot.f / mdot.s31;
%     [cp_s4,gamma_s4] = unFAIR2(f.s4);
    
    mdot.s41 = mdot.s4 + mdot.se1;
    f.s41 = f.s4*mdot.s4 / mdot.s41;
%     [cp_s41,gamma_s41] = unFAIR2(f.s41);
    
    mdot.s45 = mdot.s41 + mdot.se2;
    f.s45 = f.s41*mdot.s41 / mdot.s45;
%     [cp_s45,gamma_s45] = unFAIR2(f.s45);
    
    mdot.s6 = mdot.s45 + mdot.s13;
    f.s6 = f.s45*mdot.s45 / mdot.s6;
%     [cp_s6,gamma_s6] = unFAIR2(f.s6);

%% Inlet 

[P,T] = inlet(P,T,M0,gamma,pi_dmax);
Pr = P.o2 / 51850
% good,a fraction of a percent error from somewhere (conversion rate?)
%% Fan

[P,T] = fan(P,T,gamma,ef,pif);
Pr = P.o13 / 51850
%good
%% Low Pressure Compressor

[P,T] = LPcomp(P,T,gamma,picl,ecl);
Pr = P.o25 / 51850
%good
%% High Pressure Compressor

[P,T] = HPcomp(P,T,gamma,pich,ech)
Pr = P.o3 / 51850
%good
%% Burner
% All of this doesn't work. Made massive assumptions for working function
% beta = .01; 
% ep1 = .05;
% ep2 = .05; 
% mdot.total = 100; %kg/s
% etab = .99;
% tau_r = T.o0/T0;
% tau_cl = T.o25/T.o2;
% tau_ch = T.o3/T.o25;
% hPR = 18400*2326; %J/kg
%[T.o4,P.o4] = burneribarelyknowher(gamma,R,beta,ep1,ep2,mdot.total,etab,tau_r,tau_cl,tau_ch,hPR,T0,T.o3,P.o3);

%MAKE SURE THIS IS FIXED, SHOULD OUTPUT CP, GAMMA AS WELL)
[P,T] = burner(P,T,f,gamma,R,hPR,etab);
Pr = P.o4 / 51850
%% High Pressure Turbine

%some setup we need to setup

cp25 = 1
cp3 = 1
cp4 = 1
    
[P.o44,T.o44,f41] = HPturb(T.o25,T.o3,P.o4,T.o4,mdot,f,cp25,cp3,cp4,PtoH,etamH,etamPH,etH)
Pr = P.o44 / 51850

%% Low Pressure Turbine

[P.o5,T.o5,rhoo5] = LPturb(P.o44,T.o44,R,gamma);
Pr = P.o5 / 51850
%% Nozzle

[P.o9,T.o9,rhoo9] = nozzle(P.o5,T.o5,R,gamma);

%% Entropy Calcs
cp = (gamma*R)/(gamma-1);

ds21 = cp*log(T.o3/T.o0) - R*log(P.o3/P.o0);
dt21 = T.o3-T.o0;

ds32 = cp*log(T.o4/T.o3) - R*log(P.o4/P.o3);
dt32 = T.o4-T.o3;

ds43 = cp*log(T.o9/T.o4) - R*log(P.o9/P.o4);
dt43 = T.o9-T.o4;

h1 = figure(1);
h1.WindowStyle = 'docked';
scatter(0,0,'filled')
hold on
plot(linspace(0,ds21,10),linspace(0,dt21,10),'linewidth',2,'color','blue')
hold on

scatter(ds21,dt21,'filled')
hold on
plot(linspace(ds21,ds21+ds32,10),linspace(dt21,dt21+dt32,10),'linewidth',2,'color','blue')
hold on

scatter(ds21+ds32,dt21+dt32,'filled')
hold on
plot(linspace(ds21+ds32,ds21+ds32+ds43,10),linspace(dt21+dt32,dt21+dt32+dt43,10),'linewidth',2,'color','blue')
hold on

scatter(ds21+ds32+ds43,dt21+dt32+dt43,'filled')
hold on
plot(linspace(ds21+ds32+ds43,0,10),linspace(dt21+dt32+dt43,0,10),'--','linewidth',2,'color','blue')
title('T-s Diagram for Turbofan Engine')
xlabel('Entropy (s)')
ylabel('T.otal Temperature (T.o)')
grid('on')

%% Functions
function[P,T] = inlet(P,T,M0,gamma,pi_dmax)

if M0<1
pid = pi_dmax; 
elseif M0>1 && M0<5
pid = pi_dmax * (1-.075*((M0-1)^1.35));
else 
pid = pi_dmax * (800/(M0^4 + 935));
end

P.o2 = pid*P.o0;
[T.o2] = CPG(T.o0,P.o0,gamma, P.o2, 1);

%T.o2 = T.o0*(exp(((gamma-1)/(gamma))*(log(pid))));
%ho0 = ((gamma*R*T.o0)/(gamma-1)) + ((V0^2)/2);
%fun = @(V2) ((gamma*R*T.o2)/(gamma-1)) + ((V2^2)/2) - ho0;
%V2 = fzero(fun,V0);
end

function[P,T] = fan(P,T,gamma,ef,pif)
P.o13 = P.o2*pif^(1/ef);
[T.o13] = CPG(T.o2,P.o2,gamma, P.o13, 1);
% ho2 = ((gamma*R*T.o2)/(gamma-1)) + ((V2^2)/2);
end

function[P,T] = LPcomp(P,T,gamma,picl,ecl)
P.o25 = P.o2*picl^(1/ecl);
[T.o25] = CPG(T.o2,P.o2,gamma, P.o25, 1);
% fun = @(T.o25) (T.o25/T.o2)^((gamma*ecl)/(gamma-1)) - picl;
% T.o25 = fzero(fun,T.o2)
%T.o25 = T.o2*(exp(((gamma-1)/(gamma*ecl))*(log(picl))));
end

function[P,T] = HPcomp(P,T,gamma,pich,ech)
P.o3 = P.o25*pich^(1/ech);
[T.o3] = CPG(T.o25,P.o25,gamma, P.o3, 1);
% fun = @(T.o3) (T.o3/T.o25)^((gamma*ech)/(gamma-1)) - pich;
% T.o3 = fzero(fun,T.o25);
%T.o3 = T.o25*(exp(((gamma-1)/(gamma*ech))*(log(pich))));
end

function [P,T] = burner(P,T,f,gamma,R,hPR,etab)


cp = (gamma*R)/(gamma-1);


% fun = @(T.o4) (((mdot.4*cp*T.o4)-(mdot.31*cp*T.o3))/(mdot.f*hPR)) - etab;
% T.o4 = fzero(fun,T.o3); %gives Pr of 41 million, not good at all

% T.o4 = 3200*.555556; %gives Pr of 700 (we want 1707)

T.o4 =  3200*.555556;


gamma = 1.28
fun = @(pib) (T.o4/T.o3)^((gamma)/(gamma-1)) - pib;
pib = fzero(fun,T.o4);
P.o4 = pib*P.o3;
end

function [P,T] = HPturb(P,T,mdot,f,cp25,cp3,cp4,PtoH,etamH,etamPH,etH)
%find state across mixer
T.o41 = T.o4 * (mdot.s4 + mdot.sep1*(cp3*T.o3)/(cp4*T.o4)) / ...
    (mdot.s4 + mdot.sep1);

cp41 = 1 %fix soon based on f.s41 T
gamma41 = 1 %fix soon f.s41 and T

%fun = @(P.o41) (T.o41/T.o4)^((gamma41)/(gamma41-1)) - P.o41/P.o4;
P.o41 = fzero(fun,T.o41);

%find state across turbine
%fun = @(T.o44) mdot.s41*cp41*(T.o41 - T.o44)*etamH ... %change in enthalpy across HPturb
%    - mdot.s25*cp25*(T.o3 - T.o25)... %change in enthalpy across HP compressor
%    -PtoH / etamPH; %enthalpy draw of takeoff power
T.o44 = fzero(fun,T.o41);

%fun = @(P.o44) (T.o44/T.o41)^((gamma41*etH)/(gamma41-1)) - P.o44/P.o41;
P.o44 = fzero(fun,T.o44);
end

function [P,T] = LPturb(P,T,R,gamma)

etl = .89; %LOT3
pitl = .4; %ASSUMPTION

P.o5 = pitl*P.o44;

%fun = @(T.o5) (T.o5/T.o44)^((gamma*etl)/(gamma-1)) - pitl;
T.o5 = fzero(fun,T.o44);

T.o5 = T.o44*(exp(((gamma-1)/(gamma*etl))*(log(pitl))));

end

function [P,T] = nozzle(P,T,R,gamma)

pin = .97;

P.o9 = pin*P.o5;

T.o9 = T.o5*(exp(((gamma-1)/(gamma))*(log(pin))));

end

