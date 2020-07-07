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
To0 = T0*(1+((M0^2)*((gamma-1)/2)));
Po0 = P0*(1+((M0^2)*((gamma-1)/2)))^(gamma/(gamma-1));
% rhoo0 = rho0*(1+((M0^2)*((gamma-1)/2)))^(1/(gamma-1));

%--------Overall Performance (Givens)---------
% year = ?;
% T = ; % N (from lbf)
% mdottotal = ; % kg/s 
% F_over_mdot = T / mdottotal;
F_mdot = 62.859*9.806655; %N/kg/s from lbf/(lbf/s)
S = 1.1386*((.453592/3600)/4.44822); %kg/s/N from lbm/(lbf/s)
% pitotal = ;
    Po9_P9 = 12.745;
    Po0_P0 = Po0 / P0;
    P0_P9 = 1;
    pitotal = Po9_P9 / P0_P9 / Po0_P0; 
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

% mdotf = S*T;
% mdotbypass = mdottotal/(1+br); %after bypass leaves
% mdotburner = mdotbypass*(1-betta - ep1 -ep2); %after bleed and coolant leaves
% f = mdotf / mdotburner;
    mdotf_mdot = S*F_mdot; %unitless
    mdot25_mdot = 1/(1+alpha); %after bypass leaves
    mdot31_mdot = mdot25_mdot*(1-beta - ep1 -ep2); %after bleed and coolant leaves
    f = mdotf_mdot / mdot31_mdot


%% Inlet 

[Po2,To2,rhoo2] = inlet(Po0,To0,M0,R,gamma,pi_dmax);
Pr = Po2 / 51850
% good,a fraction of a percent error from somewhere (conversion rate?)
%% Fan

[Po13,To13,rhoo13,mdot13] = fan(Po2,To2,R,gamma,mdot.total,alpha,ef,pif);
Pr = Po13 / 51850
%good
%% Low Pressure Compressor

[Po25,To25,rhoo25] = LPcomp(Po2,To2,R,gamma,picl,ecl);
Pr = Po25 / 51850
%good
%% High Pressure Compressor

[Po3,To3,rhoo3,mdot31_mdot] = HPcomp(Po25,To25,R,gamma,mdot13,beta,ep1,ep2,pich,ech);
Pr = Po3 / 51850
%good
%% Burner
% All of this doesn't work. Made massive assumptions for working function
% beta = .01; 
% ep1 = .05;
% ep2 = .05; 
% mdottotal = 100; %kg/s
% etab = .99;
% tau_r = To0/T0;
% tau_cl = To25/To2;
% tau_ch = To3/To25;
% hPR = 18400*2326; %J/kg
%[To4,Po4] = burneribarelyknowher(gamma,R,beta,ep1,ep2,mdottotal,etab,tau_r,tau_cl,tau_ch,hPR,T0,To3,Po3);



[Po4,To4] = burner(Po3,To3,mdot31_mdot,f,gamma,R,hPR,etab);
Pr = Po4 / 51850
%% High Pressure Turbine

[Po44,To44,rhoo44] = HPturb(Po4,To4,R,gamma);

%% Low Pressure Turbine

[Po5,To5,rhoo5] = LPturb(Po44,To44,R,gamma);

%% Nozzle

[Po9,To9,rhoo9] = nozzle(Po5,To5,R,gamma);

%% Entropy Calcs
cp = (gamma*R)/(gamma-1);

ds21 = cp*log(To3/To0) - R*log(Po3/Po0);
dt21 = To3-To0;

ds32 = cp*log(To4/To3) - R*log(Po4/Po3);
dt32 = To4-To3;

ds43 = cp*log(To9/To4) - R*log(Po9/Po4);
dt43 = To9-To4;

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
ylabel('Total Temperature (To)')
grid('on')

%% Functions
function[Po2,To2,rhoo2] = inlet(Po0,To0,M0,R,gamma,pi_dmax)

if M0<1
pid = pi_dmax; 
elseif M0>1 && M0<5
pid = pi_dmax * (1-.075*((M0-1)^1.35));
else 
pid = pi_dmax * (800/(M0^4 + 935));
end

Po2 = pid*Po0;

To2 = To0*(exp(((gamma-1)/(gamma))*(log(pid))));

rhoo2 = Po2/(R*To2);

%ho0 = ((gamma*R*To0)/(gamma-1)) + ((V0^2)/2);

%fun = @(V2) ((gamma*R*To2)/(gamma-1)) + ((V2^2)/2) - ho0;
%V2 = fzero(fun,V0);



end

function[Po13,To13,rhoo13,mdot13] = fan(Po2,To2,R,gamma,mdottotal,br,ef,pif)


Po13 = Po2*pif^(1/ef);

% fun = @(To13) (To13/To2)^((gamma*ef)/(gamma-1)) - pif;
% To13 = fzero(fun,To2);

To13 = To2*(exp(((gamma-1)/(gamma*ef))*(log(pif))));

rhoo13 = Po13/(R*To13);

% ho2 = ((gamma*R*To2)/(gamma-1)) + ((V2^2)/2);

mdot13 = mdottotal/(2*br);


end

function[Po25,To25,rhoo25] = LPcomp(Po2,To2,R,gamma,picl,ecl)

Po25 = Po2*picl^(1/ecl);

% fun = @(To25) (To25/To2)^((gamma*ecl)/(gamma-1)) - picl;
% To25 = fzero(fun,To2)

To25 = To2*(exp(((gamma-1)/(gamma*ecl))*(log(picl))));

rhoo25 = Po25/(R*To2);


end

function[Po3,To3,rhoo3,mdot31] = HPcomp(Po25,To25,R,gamma,mdot13,beta,ep1,ep2,pich,ech)


Po3 = Po25*pich^(1/ech);

% fun = @(To3) (To3/To25)^((gamma*ech)/(gamma-1)) - pich;
% To3 = fzero(fun,To25);

To3 = To25*(exp(((gamma-1)/(gamma*ech))*(log(pich))));

rhoo3 = Po3/(R*To3);

mdot31 = mdot13*(1-beta-ep1-ep2);

end

function [Po4,To4] = burner(Po3,To3,mdot31,f,gamma,R,hPR,etab)


cp = (gamma*R)/(gamma-1);

mdot4 = mdot31*(1+f);
mdotf = mdot31*f;

% fun = @(To4) (((mdot4*cp*To4)-(mdot31*cp*To3))/(mdotf*hPR)) - etab;
% To4 = fzero(fun,To3); %gives Pr of 41 million, not good at all

% To4 = 3200*.555556; %gives Pr of 700 (we want 1707)

To4 =  3200*.555556;


gamma = 1.28
fun = @(pib) (To4/To3)^((gamma)/(gamma-1)) - pib;
pib = fzero(fun,To4);
Po4 = pib*Po3;
end

function [Po44,To44,rhoo44] = HPturb(Po4,To4,R,gamma)
eth = .89; %LOT3
pith = .4; %ASSUMPTION

Po44 = pith*Po4;

fun = @(To44) (To44/To4)^((gamma*eth)/(gamma-1)) - pith;
To44 = fzero(fun,To4);

To44 = To4*(exp(((gamma-1)/(gamma*eth))*(log(pith))));

rhoo44 = Po44/(R*To44);

end

function [Po5,To5,rhoo5] = LPturb(Po44,To44,R,gamma)

etl = .89; %LOT3
pitl = .4; %ASSUMPTION

Po5 = pitl*Po44;

fun = @(To5) (To5/To44)^((gamma*etl)/(gamma-1)) - pitl;
To5 = fzero(fun,To44);

To5 = To44*(exp(((gamma-1)/(gamma*etl))*(log(pitl))));

rhoo5 = Po5/(R*To5);

end

function [Po9,To9,rhoo9] = nozzle(Po5,To5,R,gamma)

pin = .97;

Po9 = pin*Po5;

To9 = To5*(exp(((gamma-1)/(gamma))*(log(pin))));

rhoo9 = Po9/(R*To9);

end


