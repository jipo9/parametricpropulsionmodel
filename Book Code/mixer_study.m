clear
clc
close all

%% Analyze M6 vs M13 at each P6/P13
pi = [60, 50, 40, 30, 20, 10, 5, 1, 0]; %P6/P13 whats a realistic range here?
M_mixer = linspace(0,1);

Pro13 = 8; %Assumption: Relative pressure at the fan
f6 = .2; %Assumption: fuel to air ratio at the mixer entry
alpha = .4; %assumption
beta = .01; %assumption

state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
state(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};

state(5,2) = {Pro13};
state(5,4) = {0};
[state] = unFAIR3(state,5);
[~,~,To13,~,~,cp13,gamma13,~,~,R13,~] = state{5,:};

state(14,4) = {f6};
mdot13 = alpha/(1+alpha); % bypass mass flow
mdot6 = (1-mdot13) * (1-beta); %mixer flow

for ii = 1:size(pi,2)
    pi_local = pi(ii);
    Pro6 = Pro13*pi_local;
    state(14,2) = {Pro6};
    state(14,3) = {[]};
    state(14,8) = {[]};
    [state] = unFAIR3(state,14);
    [~,~,To6,~,~,cp6,gamma6,~,~,R6,~] = state{14,:};

    for jj = 1:100
        M6 = M_mixer(jj);
        [M13] = Kutta_mach(gamma13,M6,gamma6,pi_local);
        [M13a] = Kutta_mach_alt(gamma13,M6,gamma6,pi_local);
        [A13_6] = area(gamma13, M13, To13, Pro13, mdot13, R13, gamma6, M6, To6, Pro6, mdot6, R6);        
        M_bypass(ii,jj) = M13;
        M_bypass_a(ii,jj) = M13a;
        A_ratio(ii,jj) = A13_6;
    end
end

figure
plot(M_mixer,M_bypass(1,:),M_mixer,M_bypass(2,:),M_mixer,M_bypass(3,:),M_mixer,M_bypass(4,:),M_mixer,M_bypass(5,:),M_mixer,M_bypass(6,:),M_mixer,M_bypass(7,:),M_mixer,M_bypass(8,:),M_mixer,M_bypass(9,:))
xlabel('Mixer Mach Number')
ylabel('Bypass Mach Number')
legend('Pi = 200','Pi = 100','Pi = 50','Pi = 25')

figure
plot(M_mixer,M_bypass_a(1,:),M_mixer,M_bypass_a(2,:),M_mixer,M_bypass_a(3,:),M_mixer,M_bypass_a(4,:),M_mixer,M_bypass_a(5,:),M_mixer,M_bypass_a(6,:),M_mixer,M_bypass_a(7,:),M_mixer,M_bypass_a(8,:),M_mixer,M_bypass_a(9,:))
xlabel('Mixer Mach Number')
ylabel('Bypass Mach Number')
legend('Pi = 200','Pi = 100','Pi = 50','Pi = 25')

figure
plot(M_mixer,A_ratio(1,:),M_mixer,A_ratio(2,:),M_mixer,A_ratio(3,:),M_mixer,A_ratio(4,:))
xlabel('Mixer Mach Number')
ylabel('Area Ratio 13/6')
legend('Pi = 200','Pi = 100','Pi = 50','Pi = 25')
%% Analyze M6 vs A6/A13 at each P6/P13






P_Pt = (1 - (gamma - 1)/2*M^2)^(-gamma/(gamma-1));
T_Tt = (1 - (gamma - 1)/2*M^2)^(-gamma/(gamma-1));
A = mdot*sqrt(Tt)/Pt    *   1/M * sqrt(R/gamma) * 1/P_Pt * sqrt(T_Tt);




for ii = 1:size(pi,2)
    pi_local = pi(ii);
    for jj = 1:100
        M6 = M_mixer(jj);
        M13 = sqrt(...
            (2/(gamma13 - 1))  *  ...
            ((pi_local * (1 + (gamma6-1)*M6^2 /2)^(-gamma6/(gamma6-1)))^ ((gamma13 - 1)/gamma13)  - 1));    
%         fun = @(M13) pi_local - ...
%             ((1+(gamma13 - 1)/2 * M13^2)^(-gamma13/(gamma13-1))) /...
%             ((1+(gamma6  - 1)/2 * M6 ^2 )^(-gamma6 /(gamma6 -1)));
%         M13 = fzero(fun,pi_local);
        M_bypass(ii,jj) = M13;
    end
end

plot(M_mixer,M_bypass(1,:),M_mixer,M_bypass(2,:),M_mixer,M_bypass(3,:),M_mixer,M_bypass(4,:))
xlabel('Mixer Mach Number')
ylabel('Bypass Mach Number')
legend('Pi = 200','Pi = 100','Pi = 50','Pi = 25')



%% Analyze M6 vs tau, pi M at each P6/P13

%% Analyze M6 vs tau, pi m1 locally

%% Analyze M6 vs tau, pi tH locally

%% Analyze M6 vs tau, pi m2 locally

%% Analyze M6 vs tau, pi tC locally

%% Analyze T5/T13 vs alpha at each M6

%% Look into bypass design?

function [A1_2] = area(gamma1, M1, Tt1, Pt1, mdot1,R1,gamma2, M2, Tt2, Pt2, mdot2,R2)
P_Pt1 = (1 - (gamma1 - 1)/2*M1^2)^(-gamma1/(gamma1-1));
T_Tt1 = (1 - (gamma1 - 1)/2*M1^2)^(-gamma1/(gamma1-1));
Ar1 = mdot1*sqrt(Tt1)/Pt1    *   1/M1 * sqrt(R1/gamma1) * 1/P_Pt1 * sqrt(T_Tt1); %Relative pressure means relative area

P_Pt2 = (1 - (gamma2 - 1)/2*M2^2)^(-gamma2/(gamma2-1));
T_Tt2 = (1 - (gamma2 - 1)/2*M2^2)^(-gamma2/(gamma2-1));
Ar2 = mdot2*sqrt(Tt2)/Pt2    *   1/M2 * sqrt(R2/gamma2) * 1/P_Pt2 * sqrt(T_Tt2);

A1_2 =Ar1/Ar2;
end





function [M1] = Kutta_mach(gamma1,M2,gamma2,pi2_1)
%calculates the mach if kutta condition is satisfied between the two states
M1 = sqrt(...
    (2/(gamma1 - 1))  *  ...
    ((pi2_1 * (1 + (gamma2-1)*M2^2 /2)^(-gamma2/(gamma2-1)))^ ((gamma1 - 1)/gamma1)  - 1));
end


function [M1] = Kutta_mach_alt(gamma1,M2,gamma2,pi2_1)
%calculates the mach if kutta condition is satisfied between the two states
P_Pt2 = (1 - (gamma2 - 1)/2*M2^2)^(-gamma2/(gamma2-1));
P_Pt1 = P_Pt2*pi2_1;
M1 = sqrt((2/(gamma1 - 1)) * (P_Pt1 ^ ((gamma1 - 1)/gamma1)  - 1));
end
