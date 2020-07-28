clear
clc
close all

%idk if were using the right pi
%alpha seems very off
%does gamma vary over the same station at different mach?
MFP6 = MFP2(.3835, 1.3041, 287)
MFP16 = MFP2(.4559, 1.3871, 287)

%% Analyze M6 vs M16 at each P6/P16 (fixed alpha)
pi_b = [.9,.92,.94,.96,1]; %pi_burner
pi = 1./pi_b; %P16/P6 whats a realistic range here?
M_mixer = linspace(0,1);

Pro16 = 8; %Assumption: Relative pressure at the fan. Has almost no effect on M16 vs M6, small noticeable (max of 5%) effect on Area Ratio
f = .03; %Assumption: fuel to air ratio at the mixer entry
alpha = .449; %Start Assumption
beta = .01; %Assumption
ep1 = .05; %Assumption
ep2 = .05; %Assumption


state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
state(2:22,1) = {'0';'o0';'o2';'o16';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
[state] = derived_parameters(state,alpha,beta,ep1,ep2,f);

state(5,2) = {Pro16};
state(5,4) = {0};
[state] = unFAIR3(state,5);
[~,~,To16,~,~,cp16,gamma16,~,~,R16,~] = state{5,:};

mdot6 = state{13,5};
mdot16 = state{5,5};
alpha_prime = mdot16/mdot6;

for ii = 1:size(pi,2)
    pi_local = pi(ii);
    Pro6 = Pro16/pi_local;
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
        [alpha_prime] = bypass(MFP16,Pro16,To16,MFP6,Pro6,To6,A16_6);
        [alpha_prime_approx] = M16/M6*A16_6/pi_local^2;
        M_bypass(ii,jj) = M16;
        A_ratio(ii,jj) = A16_6;
        alpha_p(ii,jj) = alpha_prime;
        alpha_p_approx(ii,jj) = alpha_prime_approx;

    end
end

figure
plot(M_mixer,M_bypass(1,:),M_mixer,M_bypass(2,:),M_mixer,M_bypass(3,:),M_mixer,M_bypass(4,:),M_mixer,M_bypass(5,:),[.4,.3835],[.394,.4559],'x')
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




%% 
clear
clc
close all

%% Analyze M6 vs M16 at each P6/P16 (fixed A16_6)
A16_6 = 0.2485;
pi = [1.0042,1.0492];
M_mixer = linspace(.01,1);

Pro16 = 8; %Assumption: Relative pressure at the fan. Has almost no effect on M16 vs M6, small noticeable (max of 5%) effect on Area Ratio
f = .03; %Assumption: fuel to air ratio at the mixer entry
beta = .01; %Assumption
ep1 = .05; %Assumption
ep2 = .05; %Assumption

alpha_ref = .5;

state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
state(2:22,1) = {'0';'o0';'o2';'o16';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
[state] = derived_parameters(state,alpha_ref,beta,ep1,ep2,f);

state(5,2) = {Pro16};
state(5,4) = {0};
[state] = unFAIR3(state,5);
[~,~,To16,~,~,cp16,gamma16,~,~,R16,~] = state{5,:};

mdot6 = state{13,5};
mdot16 = state{5,5};
alpha_prime = mdot16/mdot6;

for ii = 1:size(pi,2)
    pi_local = pi(ii);
    Pro6 = Pro16*pi_local;
    state(14,2) = {Pro6};
    state(14,3) = {[]};
    state(14,8) = {[]};
    [state] = unFAIR3(state,14);
    [~,~,To6,~,~,cp6,gamma6,~,~,R6,~] = state{14,:};
    
    for jj = 1:100
        error =1;
        alpha = alpha_ref;
        while error > .0001
            M6 = M_mixer(jj);
            [M16] = Kutta_mach(gamma16,M6,gamma6,pi_local);
            
            MFP16 = MFP2(M16, gamma16, R16);
%             MFP = mdot*sqrt(Tt) / (Pt*A);
            MFP6 = MFP2(M6, gamma6, R6);
            
            [alpha_prime] = bypass(MFP16,Pro16,To16,MFP6,Pro6,To6,A16_6);
            
            fun = @(alpha) alpha_prime - alpha ./ ((1-beta-ep1-ep2)*(1+f) + ep1 +ep2);
            alpha_i = fzero(fun,alpha_prime);
            error = norm((alpha-alpha_i)/alpha);
            alpha = alpha_i;
            
            [state] = derived_parameters(state,alpha,beta,ep1,ep2,f);
            state(5,2) = {[]};
            state(5,8) = {[]};
            [state] = unFAIR3(state,5);
            [~,~,To16,~,~,cp16,gamma16,~,~,R16,~] = state{5,:};
            pi_local = pi(ii);
            Pro6 = Pro16*pi_local;
            state(14,2) = {[]};
            state(14,8) = {[]};
            [state] = unFAIR3(state,14);
            [~,~,To6,~,~,cp6,gamma6,~,~,R6,~] = state{14,:};
        end
        [pi_M] = mixer_ideal(state,alpha_prime,A16_6,M16,M6);
        M_bypass(ii,jj) = M16;
        alpha_p(ii,jj) = alpha_prime;
        alpha_(ii,jj) = alpha;
        pi_M_ideal(ii,jj) = pi_M;
    end
end



figure
plot(M_mixer,M_bypass(1,:),M_mixer,M_bypass(2,:),[.4,.3835],[.394,.4559],'x')
xlabel('Mixer Mach Number')
ylabel('Bypass Mach Number')
legend('reference','test')


figure
plot(M_mixer(30:70),alpha_p(1,30:70),M_mixer(30:70),alpha_p(2,30:70))
xlabel('Mixer Mach Number')
ylabel('Alpha Prime')
legend('reference','test')

figure
plot(M_mixer(30:70),alpha_(1,30:70),M_mixer(30:70),alpha_(2,30:70),[.4,.3835],[.449,.520],'x')
xlabel('Mixer Mach Number')
ylabel('Alpha')
legend('reference','test')

figure
plot(M_mixer(30:70),pi_M_ideal(1,30:70),M_mixer(30:70),pi_M_ideal(2,30:70),[.4,.3835],[.9935,1],'x')
xlabel('Mixer Mach Number')
ylabel('Pi_mixer_ideal')
legend('reference','test')
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

function [alpha_prime] = bypass(MFP1,Pt1,Tt1,MFP2,Pt2,Tt2,A1_2)
alpha_prime = Pt1/Pt2 * MFP1/MFP2 * sqrt(Tt2/Tt1) * A1_2;
end

function [pi] = mixer_ideal(state,alpha_prime,A16_6,M16,M6)
[~,~,To16,~,~,cp16,gamma16,ho16,~,R16,~] = state{5,:};
[~,~,To6,~,mdot6,cp6,gamma6,ho6,~,R6,~] = state{14,:};

A6_6A = 1/(1 + A16_6);

ho6A = (ho6 + alpha_prime*ho16)/(1+alpha_prime);
state(15,2) = {[]};
state(15,3) = {[]};
state(15,8) = {ho6A};
[state] = unFAIR3(state,15);
[~,~,To6A,~,mdot6A,cp6A,gamma6A,ho6A,~,R6A,~] = state{15,:};

[T_Tt16] = temperature(M16,gamma16);
[T_Tt6] = temperature(M6,gamma6);
T16 = To16 * T_Tt16;
T6 = To6 * T_Tt6;

const = 1/(1+alpha_prime) * ...
    (sqrt(R6*T6/gamma6)*(1+gamma6*M6^2)/(M6)   +   alpha_prime*sqrt(R16*T16/gamma16)*(1+gamma16*M16^2)/(M16));

M6Ai = M6;
error = 1;

while error > .0001
    [T_Tt6A] = temperature(M6Ai,gamma6A);
    T6A = T_Tt6A * To6A;
    M6A = sqrt(R6A*T6A/gamma6A)*(1+gamma6A*M6Ai^2)/(const);
    error = norm((M6Ai - M6A)/M6Ai);
    M6Ai = M6A;
end

[MFP6] = MFP2(M6, gamma6, R6);
[MFP6A] = MFP2(M6A, gamma6A, R6A);

pi = (1+alpha_prime) * sqrt(To6A/To6) * A6_6A * MFP6 / MFP6A;
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

function [state] = derived_parameters(state,alpha,beta,ep1,ep2,f)
%% Derived Parameters
% COOLING AIR CALCULATIONS
% Found in Aircraft Engine Design - Mattingly but
% unable to find rationale...

mdot = 1;



% MASS FLOW AND FUEL TO AIR CALCULATIONS
f0 = 0; %freestream fuel/air ratio

mdot25 = mdot/(1+alpha); %after bypass leaves
mdot13 = mdot - mdot25; % bypass mass flow

mdot31 = mdot25*(1-beta - ep1 -ep2); %after bleed and coolant leaves
mdotbeta = mdot25*beta; %bleed air
mdotep1 = mdot25*ep1; %coolant air 1
mdotep2 = mdot25*ep2; %coolant air 2
mdotep = mdotep1+mdotep2; %combined coolant air

mdot4 = mdot31 * (1+f); %mass flow rate post-burner
mdot_f = mdot31*f;
f4 = mdot_f / mdot31; %fuel/air ratio post-burner

mdot41 = mdot4 + mdotep1; %mass flow rate after addtion of cooling air 1
f41 = f4*mdot4 / mdot41; %fuel/air ratio after addtion of cooling air 1

mdot45 = mdot41 + mdotep2; %mass flow rate after addtion of cooling air 2
f45 = f41*mdot41 / mdot45; %fuel/air ratio after addtion of cooling air 2

mdot6A = mdot45 + mdot13; %mass flow rate after addtion of bypass air
f6A = f45*mdot45 / mdot6A; %fuel/air ratio after addtion of bypass air

% Store all values
state(2:8,4) = {f0};
state(2:4,5) = {mdot};
state(5,5) = {mdot13};
state(6:7,5) = {mdot25};
state(8,5) = {mdot31};
state(9,4) = {f4};
state(9,5) = {mdot4};
state(10:11,4) = {f41};
state(10:11,5) = {mdot41};
state(12:14,4) = {f45};
state(12:14,5) = {mdot45};
state(15:18,4) = {f6A};
state(15:18,5) = {mdot6A};
state(19,5) = {mdotbeta};
state(20,5) = {mdotep};
state(21,5) = {mdotep1};
state(22,5) = {mdotep2};

end