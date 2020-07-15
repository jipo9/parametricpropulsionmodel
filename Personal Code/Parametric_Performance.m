clc
clear
close all

%Assume no afterburner
%% Initialize cells

state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)'};
state(2:18,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9'};
component = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
component(2:17,1) = {'Ram Recovery';'Inlet(Ideal,Actual)';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
design = {'Parameter','Value'};
design(2:8,1) = {'alpha';'beta';'epsilon';'empty';'Pto';'Pi_total';'h_PR'};

%% Initial Conditions
%--------Overall Performance (Givens)---------
S = 3.2251e-05; %kg/s/N
T = 5.5923e+04; % N
mdot0 = 90.7185; %kg/s
alt = 35000/3.281; %m
M0 = 1.6;
pitotal = 82.9635;
P0_P9 = 1;

alpha = .4; %bypass ratio
beta = .01; %bleed ratio
Pto = 301.34*10^3; %Watts
h_PR = 18400*2326; %J/kg, for a (CH2)n propellant

design(2,2) = {alpha};
design(3,2) = {beta};
design(6,2) = {Pto};
design(7,2) = {pitotal};
design(8,2) = {h_PR};



year = 1995;

index_diffuser = 3;     %1 = subsonic aircraft w/ engines in nacelles
                        %2 = subsonic aircraft w/ engines in airframe 
                        %3 = supersonic aircraft w/ engines in airframe 
index_turbine = 2;      %1 = uncooled turbine
                        %2 = cooled turbine
index_nozzle = 2;       %1 = Fixed-area convergent nozzle
                        %2 = Variable-area convergent nozzle
                        %3 = Variable-area convergent-divergent nozzle
%-------------Level of tech------------
pi_d = [.9,.95,.98,.995;
        .88,.93,.96,.97; 
        .85,.90,.94,.96];
e_c = [.8,.84,.88,.9];
e_f = [.78,.82,.86,.89];
pi_b = [.9,.92,.94,.96];
eta_b =[.88,.94,.99,.995];
e_t =[.8,.85,.89,.91;
      .75,.83,.87,.89];
% pi_aB = [.9,.92,.94,.95];
% eta_AB = [.85,.91,.96,.97];
pi_n = [.95,.97,.98,.995;
        .93,.96,.97,.985;
        .9,.93,.95,.98];
T_t4max = [1110,1390,1780,2000]; %kelvin
T_t7max = [1390,1670,2000,2200]; %kelvin


%-------------Finding our values------------
ii = floor((year -1925)/20);

%Nozzle
pi_dmax = pi_d(index_diffuser,ii);
component(3,2) = {pi_dmax};

%Fan
pif = 3.8;
e_f = e_f(ii);
component(4,2:3) = {pif,e_f};

%Low compressor
picl = 3.8;
e_c = e_c(ii);
component(5,2:3) = {picl,e_c};

%High compressor
pich = 4.2105;
component(6,2:3) = {pich,e_c};

%Burner
pi_b = pi_b(ii);
eta_b = eta_b(ii);

%High turbine
etamH = 1; %assuming all turbine eta = 1 gives an increase of thrust of +.1% (super negligeble)
etamPH = 1; %assuming all turbine eta = 1 gives an increase of thrust of +.1% (super negligeble)
component(9,5) = {etamH};
component(10,5) = {etamPH}; 
e_t = e_t(index_turbine,ii);
T_t4 = T_t4max(ii);
state(9,3) = {T_t4};


if T_t4 > 2400*.5556 && index_turbine == 2
    ep = 2*(T_t4/.5556-2400)/(16000);    %bypass ratio for mixer 1
    design(4,2) = {ep};
else
    ep = 0;
    design(4,2) = {ep};
end

%Low turbine
etamL = 1; %assuming all turbine eta = 1 gives an increase of thrust of +.1% (super negligable)
etamPL = 1; %assuming all turbine eta = 1 gives an increase of thrust of +.1% (super negligable)
component(12,5) = {etamL};
component(13,5) = {etamPL}; 

%Mixer
pi_M_max = 1; %assuming pi = 1 gives an increase of thrust of +.5% (negligable)
component(14,2) = {pi_M_max};

%Nozzle
pi_n = pi_n(index_nozzle, ii);
component(16,2) = {pi_n};
% T_t7max = T_t7max(ii);

%Assume 100% polytropic efficiency outside of turbines
component(2:3,3) = {1}; 
component(7:8,3) = {1};
component(11,3) = {1};
component(14:16,3) = {1};


F_mdot = T / mdot0;
%% Temporary Override
alt = 35000/3.281; %altitude [m from feet]
M0 = 1.6; %freestream mach number
F_mdot = 62.859*9.806655; %thrust/mdot [N/kg/s from lbf/(lbf/s)]
S = 1.1386*((.453592/3600)/4.44822); %specific fuel consuption[kg/s/N from lbm/(lbf/s)]
T_t4 = 3200*.5556; %max burner temperature [R to K]
Po9_P9 = 12.745; % 

alpha = .4; %bypass ratio
beta = .01; %bleed ratio
PtoH = 301.34*10^3; %power takeoff high spool [watts]
PtoL = 0; %power takeoff low spool [watts]
h_PR = 18400*2326; %fuel heating value for a (CH2)n propellant [J/kg]

design(2,2) = {alpha}; %store values in design
design(3,2) = {beta};
design(7,2) = {PtoH};
design(6,2) = {PtoL};
design(8,2) = {h_PR};


pi_dmax = .96; %diffuser pressure ratio

pif = 3.8; %fan pressure ratio
ef = .89; %fan polytropic efficiency

picL = 3.8; %low pressure compressor pressure ratio
ecL = .89; %low pressure polytropic efficiency

picH = 4.2105; %high pressure compressor pressure ratio
ecH = .9; %high pressure polytropic efficiency

eta_b = .999; %burner efficiency
pi_b = .95; % burner pressure ratio

etH = .89; %high pressure turbine polytropic efficiency

etL = .9; %low pressure turbine polytropic efficiency

etamH = .995; %high pressure shaft mechanical efficiency
etamPH = .99; %high pressure shaft power takeoff mechancal efficiency
etamL = .995; %low pressure shaft mechanical efficiency
etamPL = 1; %low pressure shaft power takeoff mechancal efficiency

pi_M_max = .97; %mixer pressure ratio

pin = .97; %nozzle pressure ratio


component(3,2) = {pi_dmax}; %store values in component
component(4,2:3) = {pif,ef};
component(5,2:3) = {picL,ecL};
component(6,2:3) = {picH,ecH};
component(7,2) = {pi_b};
component(9,5) = {etamH};
component(10,5) = {etamPH}; 
component(9,3) = {etH};
component(12,5) = {etamL};
component(13,5) = {etamPL}; 
component(12,3) = {etL};
component(14,2) = {pi_M_max};
component(16,2) = {pin};
component(2:3,3) = {1};
component(7:8,3) = {1};
component(11,3) = {1};
component(14:16,3) = {1};

if T_t4 > 2400*.55556 %cooling air calculations
    ep1 = (T_t4/.5556-2400)/(16000);    
    ep2 = ep1;                         
else
    ep1 = 0;
    ep2 = 0;
end

design(4,2) = {ep1};
design(5,2) = {ep2};
%% Mass flow and Air Props
% mdot.f = S*T;
% mdot.bypass = mdot.total/(1+br); %after bypass leaves
% mdot.burner = mdot.bypass*(1-betta - ep1 -ep2); %after bleed and coolant leaves
% f = mdot.f / mdot.burner;

% for now using all mdot/mdot total, in actual it will include units
mdot_f = S*T /eta_b; %unitless

f0 = 0; %freestream fuel/air ratio 

mdot25 = mdot0/(1+alpha); %after bypass leaves
mdot13 = mdot0 - mdot25; % bypass mass flow

mdot31 = mdot25*(1-beta - ep); %after bleed and coolant leaves
mdotbeta = mdot25*beta; %bleed air
mdotep = mdot25*ep;

mdot4 = mdot31 + mdot_f; %mass flow rate post-burner
f4 = mdot_f / mdot31; %fuel/air ratio post-burner

mdot45 = mdot4 + mdotep; %mass flow rate after addtion of cooling air 2
f45 = f4*mdot4 / mdot45; %fuel/air ratio after addtion of cooling air 2

mdot6A = mdot45 + mdot13; %mass flow rate after addtion of bypass air
f6A = f45*mdot45 / mdot6A; %fuel/air ratio after addtion of bypass air


state(2:8,4) = {f0};
state(2:4,5) = {mdot0};
state(5,5) = {mdot13};
state(6:7,5) = {mdot25};
state(8,5) = {mdot31};
state(9,4) = {f4};
state(9,5) = {mdot4};
state(12:13,4) = {f45};
state(12:13,5) = {mdot45};
state(14:18,4) = {f6A};
state(14:18,5) = {mdot6A};
%% Solve for pi's

if alpha > .2 && alpha < 1 %low alpha
    pi_c = [8,30];
    pi_f = [2,4];
elseif alpha > 5 && alpha < 10 %high alpha
    pi_c = [30,40];
    pi_f = [1.4,4];
else
    error('Invalid range, model expansion or inferences required')
end

n = 20;
pif = linspace(pi_f(1),pi_f(2),n);
pic = linspace(pi_c(1),pi_c(2),n); 
picl = 0.2375*pic;
pich = pic/picl;
ep = .05;
T_t4 = linspace(1600,1900,n);
for ii = 1:n
    for jj = 1:n
        for kk = 1:n
        component(4:6,2) = {pif(ii);picl(jj);picl(jj)};
        state(9,3) = {T_t4(kk)};
        [statei, componenti, performancei, designi, errori] = Parametric_Function(state, design, component, M0, alt, mdotep,F_mdot,S);
        T_error(ii,jj,kk) = errori(1,1);
        S_error(ii,jj,kk) = errori(1,2);
        eff_T_error(ii,jj,kk) = errori(1,3);
        eff_P_error(ii,jj,kk) = errori(1,4);
        err_total = T_error(ii,jj,kk)*eff_T_error(ii,jj)*eff_P_error(ii,jj);
        end
    end
    ii
end



err_total = (1-S_error).*(1-eff_T_error).*(1-eff_P_error);
err_total = abs(err_total);
%%
[M,I] = min(err_total,[],[1 2, 3],'linear');

m = 0;
ii = 0;
jj = 0;
kk = 0;
%%
m = 0;
ii = 0;
jj = 0;
kk = 0;
while m ~= I
  ii = ii+1;
  while jj>n
    jj = jj+1;
    while kk>n
          kk = kk+1;
          m = m+1
    end
    kk = 0;
  end
  jj = 0;
  if ii >20
      break
  end
end

%%
figure
hold on
contour(pif,pic,err_total(:,:,1))
plot(pif(90),pic(37),'x')
plot(pif(50),pic(50),'x')
title('Thrust')

figure
hold on
contour(pif,pic,eff_T_error)
plot(pif(90),pic(37),'x')
plot(pif(50),pic(50),'x')
title('Prop eff')

figure
hold on
contour(pif,pic,eff_P_error)
plot(pif(90),pic(37),'x')
plot(pif(50),pic(50),'x')
title('Thermal eff')


figure
hold on
err_total = T_error.*eff_T_error.*eff_P_error;
contour(pif,pic,err_total)
plot(pif(90),pic(37),'x')
plot(pif(50),pic(50),'x')
title('Total eff')


%define range of compressor and fan pressure ratios
% [state, component, performance, design, error] = Turbofan_Barebones_iteration(state, design, component, M0, alt, mdotep,F_mdot,S)
%Find minimum error for T at each fan pressure
% %Fan
% pif = 3.8;
% component(4,2) = {pif,e};
% 
% %Low compressor
% picl = 3.8;
% component(5,2) = {picl};
% 
% %High compressor
% pich = 4.2105;
% component(6,2) = {pich};

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
err_T_mdot = performance{2,1} /F_mdot %T/mdot error compared to book
err_s = performance{2,2} / S %S error compared to book
err_efftherm = performance{2,3} / .5589 %thermal efficiency error compared to book
err_effprop =performance{2,4} / .6162 %propulsive efficiency compared to book
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



F_mdot = (1+f_0-(beta/(1+alpha)))*v9     -   v0  +   (1+f_0-(beta/(1+alpha)))*R9*T9*(1-Pr0/Pr9)/(R0*v9*gamma0);
S = f_0 / F_mdot;
eta_TH = ((v0^2/2*((1+f_0-(beta/(1+alpha)))*(v9/v0)^2 - 1) + (Pto)/mdot0))/...
    (f_0*h_PR);
%eta_P = 2/(1+v9/v0); Simplified case, neglected for now
eta_P = (2*F_mdot/v0)/...
    ((1+f_0-beta/(1+alpha))*(v9/v0)^2 - 1);
eta_o = eta_TH*eta_P;
% eta_o = v0/(h_PR*S); Alternate eqn
performance(1,:) = {'Thrust','Specific Fuel Consumption','Propulsive Efficiency','Thermal Efficiency','Overall Efficiency','Exhaust Mach'};
performance(2,:) = {F_mdot,S,eta_TH,eta_P,eta_o,M9};
end

