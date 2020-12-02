clear
clc
close all


alpha = 7.5;
pi_total = 35;
pi_f = 3;


n = 30;
M_SL = linspace(.1,1,n);
M_10 = linspace(.1,1,n);
M_20 = linspace(.1,1,n);
M_30 = linspace(.1,1,n);
M_36 = linspace(.1,1,n);
M_40 = linspace(.1,1,n);
M_50 = linspace(.1,1,n);

M = [M_SL;M_10;M_20;M_30;M_36;M_40;M_50];
altitude = [0,10000,20000,30000,36000,40000,50000] ./ 3.281;

figure
hold on
F_A = zeros(size(M));
for ii = 1:size(altitude,2)
    ii
    alt = altitude(ii);
    for jj = 1:size(M,2)
        jj
        M0 = M(ii,jj);
        [F_A0] = simplistic_model_(alpha,pi_total,pi_f,M0,alt);
        F_A(ii,jj) = F_A0;
    end
    plot(M(ii,:),F_A(ii,:),'linewidth',1.5)
    legend('SL','10k ft','20k ft','30k ft','36k ft','40k ft','50k ft','location','best')
    title('Mach Number vs. Thrust/Inlet Area at Various Altitudes')
    xlabel('Mach Number')
    ylabel('F/A0')
    grid('on')
    warning('off','all')
end


function [F_A] = simplistic_model_(alpha,pi_total,pi_f,M0,alt)
state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
state(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};

gamma = 1.4;
R = 287;


%% Calculate Velocities
[T0, ~, P0, rho0] = atmosisa(alt);
state(2,3) = {T0};
state(2,4) = {0};
[state] = unFAIR3(state,2);
[~,~,T0,~,~,cp0,gamma0,~,~] = state{2,:};
a0 = sqrt(R*gamma*T0); %m/s
v0 = M0*a0;

To0 = T0*(1+((M0^2)*((gamma0-1)/2))); %find total temperature using isentropic
state(3,3) = {To0};
state(3,4) = {0};
[state] = unFAIR3(state,3);
Pro0 = state{3,2};


Pro19 = Pro0*pi_f;
P_Pt19 = Pro19/state{2,2};
state(4,2) = {Pro19};
state(4,4) = {0};
[state] = unFAIR3(state,4);
M19 = sqrt((( P_Pt19^(-(gamma-1)/gamma)  -  1 ) / (- (gamma - 1)/2)));
To19 = state{4,3};
T19 = To19 * (1 - (gamma - 1)/2*M19^2)^-1;
a19 = sqrt(R*gamma*T19); %m/s
v19 = M19*a19;


Pro9 = Pro0*pi_total;
P_Pt9 = Pro9/state{2,2};
state(5,2) = {Pro9};
state(5,4) = {0};
[state] = unFAIR3(state,5);
M9 = sqrt((( P_Pt9^(-(gamma-1)/gamma)  -  1 ) / (- (gamma - 1)/2)));
To9 = state{5,3};
T9 = To9 * (1 - (gamma - 1)/2*M19^2)^-1;
a9 = sqrt(R*gamma*T9); %m/s
v9 = M9*a9;
%% Calculate Mass Flows
mdot0_A = rho0*v0;
mdot9_A = mdot0_A/(1+alpha);
mdot19_A = mdot0_A - mdot9_A;

F_A = mdot19_A*v19 +mdot9_A*v9 - mdot0_A*v0;
end



% P_Pt = (1 - (gamma - 1)/2*M^2)^(-gamma/(gamma-1));
% M = sqrt((( P_Pt^(-(gamma-1)/gamma)  -  1 ) / (- (gamma - 1)/2)));





% F = mdot19*v19 +mdot9*v9 - mdot0*v0;
