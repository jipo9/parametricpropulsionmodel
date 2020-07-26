clear
clc
close all

%idk if were using the right pi
%alpha seems very off
%does gamma vary over the same station at different mach?

%% Analyze M6 vs M16 at each P6/P16
pi_b = [.9,.92,.94,.96,1]; %pi_burner
pi = 1./pi_b; %P6/P16 whats a realistic range here?
M_mixer = linspace(0,1);

Pro16 = 8; %Assumption: Relative pressure at the fan. Has almost no effect on M16 vs M6, small noticeable (max of 5%) effect on Area Ratio
f6 = .2; %Assumption: fuel to air ratio at the mixer entry
alpha = .4; %assumption
beta = .01; %assumption

state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
state(2:22,1) = {'0';'o0';'o2';'o16';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};

state(5,2) = {Pro16};
state(5,4) = {0};
[state] = unFAIR3(state,5);
[~,~,To16,~,~,cp16,gamma16,~,~,R16,~] = state{5,:};

state(14,4) = {f6};
mdot16 = alpha/(1+alpha); % bypass mass flow
mdot6 = (1-mdot16) * (1-beta); %mixer flow

for ii = 1:size(pi,2)
    pi_local = pi(ii);
    Pro6 = Pro16*pi_local;
    state(14,2) = {Pro6};
    state(14,3) = {[]};
    state(14,8) = {[]};
    [state] = unFAIR3(state,14);
    [~,~,To6,~,~,cp6,gamma6,~,~,R6,~] = state{14,:};

    for jj = 1:100
        M6 = M_mixer(jj);
        [M16] = Kutta_mach(gamma16,M6,gamma6,pi_local);
        MFP16 = MFP2(M16, gamma16, R16);
        MFP6 = MFP2(M6, gamma6, R6);
        [A16_6] = area(gamma16, M16, To16, Pro16, mdot16,MFP16, gamma6, M6, To6, Pro6, mdot6,MFP6);
        [alpha_prime] = bypass(R16,gamma16,M16,Pro16,To16,R6,gamma6,M6,Pro6,To6,A16_6);
        [alpha_prime_approx] = M16/M6*A16_6/pi_local^2;
        M_bypass(ii,jj) = M16;
        A_ratio(ii,jj) = A16_6;
        alpha_p(ii,jj) = alpha_prime;
        alpha_p_approx(ii,jj) = alpha_prime_approx;

    end
end

figure
plot(M_mixer,M_bypass(1,:),M_mixer,M_bypass(2,:),M_mixer,M_bypass(3,:),M_mixer,M_bypass(4,:),M_mixer,M_bypass(5,:))
xlabel('Mixer Mach Number')
ylabel('Bypass Mach Number')
legend('1945-1965','1965-1985','1985-2005','2005-2025','pi = 1')


figure
plot(M_mixer(30:70),A_ratio(1,30:70),M_mixer(30:70),A_ratio(2,30:70),M_mixer(30:70),A_ratio(3,30:70),M_mixer(30:70),A_ratio(4,30:70),M_mixer(30:70),A_ratio(5,30:70))
xlabel('Mixer Mach Number')
ylabel('Area Ratio 16/6')
legend('1945-1965','1965-1985','1985-2005','2005-2025','pi = 1')

figure
plot(M_mixer(30:70),alpha_p(1,30:70),M_mixer(30:70),alpha_p(2,30:70),M_mixer(30:70),alpha_p(3,30:70),M_mixer(30:70),alpha_p(4,30:70),M_mixer(30:70),alpha_p(5,30:70))
hold on
plot(M_mixer(30:70),alpha_p_approx(1,30:70),M_mixer(30:70),alpha_p_approx(2,30:70),M_mixer(30:70),alpha_p_approx(3,30:70),M_mixer(30:70),alpha_p_approx(4,30:70),M_mixer(30:70),alpha_p_approx(5,30:70))

xlabel('Mixer Mach Number')
ylabel('Alpha Prime')
% legend('1945-1965','1965-1985','1985-2005','2005-2025','pi = 1')

%% Analyze M6 vs tau, pi M at each P6/P16

%% Analyze M6 vs tau, pi m1 locally

%% Analyze M6 vs tau, pi tH locally

%% Analyze M6 vs tau, pi m2 locally

%% Analyze M6 vs tau, pi tC locally

%% Analyze T5/T16 vs alpha at each M6

%% Look into bypass design?

function [A1_2] = area(gamma1, M1, Tt1, Pt1, mdot1,MFP1,gamma2, M2, Tt2, Pt2, mdot2,MFP2)
P_Pt1 = pressure(M1,gamma1);
T_Tt1 = temperature(M1,gamma1);
Ar1 = mdot1*sqrt(Tt1)/Pt1    *   MFP1; %Relative pressure means relative area

P_Pt2 = pressure(M2,gamma2);
T_Tt2 = temperature(M2,gamma2);
Ar2 = mdot2*sqrt(Tt2)/Pt2    *   MFP2;

A1_2 =Ar1/Ar2;
end


function [M1] = Kutta_mach(gamma1,M2,gamma2,pi2_1)
%calculates the mach if kutta condition is satisfied between the two states
P_Pt2 = pressure(M2,gamma2);
P_Pt1 = P_Pt2*pi2_1;
M1 = sqrt((2/(gamma1 - 1)) * (P_Pt1 ^ ((gamma1 - 1)/gamma1)  - 1));
end

function [alpha_prime] = bypass(R1,gamma1,M1,Pt1,Tt1,R2,gamma2,M2,Pt2,Tt2,A1_2)
    P_Pt1 = (1 - (gamma1 - 1)/2*M1^2)^(-gamma1/(gamma1-1));
    T_Tt1 = (1 - (gamma1 - 1)/2*M1^2)^-1;
    MFP1 = M1*sqrt(gamma1/R1)*sqrt(T_Tt1)/P_Pt1;
%     MFP = MFP2(M1, gamma1, R1)
    mdot1 = MFP1*Pt1/sqrt(Tt1);
    
    
    P_Pt2 = (1 - (gamma2 - 1)/2*M2^2)^(-gamma2/(gamma2-1));
    T_Tt2 = (1 - (gamma2 - 1)/2*M2^2)^-1;
    MFP2 = M2*sqrt(gamma2/R2)*sqrt(T_Tt2)/P_Pt2;
    mdot2 = MFP2*Pt2/sqrt(Tt2);
    
    alpha_prime = mdot1/mdot2*A1_2;
end

function [MFP] = MFP1(mdot, Tt, Pt, A)
    MFP = mdot*sqrt(Tt) / (Pt*A);
end

function [MFP] = MFP2(M, gamma, R)
    [P_Pt] = pressure(M,gamma);
    [T_Tt] = temperature(M,gamma);
    MFP = M*sqrt(gamma/R)/sqrt(T_Tt)*P_Pt;
end

function [P_Pt] = pressure(M,gamma)
P_Pt = (1 - (gamma - 1)/2*M^2)^(-gamma/(gamma-1));
end

function [T_Tt] = temperature(M,gamma)
T_Tt = (1 - (gamma - 1)/2*M^2)^-1;
end