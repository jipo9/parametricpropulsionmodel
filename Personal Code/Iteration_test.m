clear
clc
close all

%% Define inputs
%input number corresponds to what the analysis is most sensitive too
alt = 35000/3.281; %5
M0 = 1.6; %2
mdot0 = 200*0.45359237;
F_mdot = 62.859*9.806655; %13
S = 1.1386*((.453592/3600)/4.44822); %14
alpha = .4; %9
beta = .01;
PtoH = 301.34*10^3; 
PtoL = 0;
h_PR = 18400*2326; %8
pi_dmax = .96;
pif = 3.8;
ef = .89;
picl = 3.8; %11
ecl = .89; %7
pich = 4.2105; %12
ech = .9; %6
eta_b = .999; %15
pi_b = .95;
etamH = .995;
etamPH = .99;
etH = .89;
etamL = .995; %10
etamPL = 1;
etL = .9;
pi_M_max = .97;
pin = .97; %4
T_t4 = 3200*.5556; %1
Po9_P9 = 12.745; %3

%% Define inputs
in = [alt M0 mdot0 F_mdot S alpha beta PtoH PtoL h_PR pi_dmax pif ef picl ecl pich ech eta_b pi_b etamH etamPH etH etamL etamPL etL pi_M_max pin T_t4 Po9_P9 ];
[state,component,performance] = turbofan_iteration(in);
T = performance{2,1};
S = performance{2,2};

for ii = 1:29
    A = in(ii);
    in(ii) = A*1.005;
    [state,component,performance] = turbofan_iteration(in);
    T_high(ii) = performance{2,1};
    S_high(ii) = performance{2,2};
    in(ii) = A*.995;
    [state,component,performance] = turbofan_iteration(in);
    T_low(ii) = performance{2,1};
    S_low(ii) = performance{2,2};
    sensitivity(ii,1) = (T_high(ii)-T)/T;
    sensitivity(ii,2) = (T_low(ii)-T)/T;
    sensitivity(ii,3) = (S_high(ii)-S)/S;
    sensitivity(ii,4) = (S_low(ii)-S)/S;
    in(ii) = A;
    ii
end



