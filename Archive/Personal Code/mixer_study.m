clear
clc
close all

%% Read Me
% This study looks at the mixer of a turbofan engine and aims to develop
% an understanding of mach ratios, area ratios, and pressure ratios between
% stations.

% With an endgoal of approximating an area ratio and mach number, known
% values are plotted against our approximations and used as a reference to
% check accuracy. Inputs come from the text, and equations are a mixture of
% self derived and found in the text
%% Notes
% Current state
    % M6/M16 seems accurate
    % Bypass ratios and ideal pressure ratio seems off
 
% To-do
    % Hardcode station 6 to what it should be
    % iron out descrepency of Pt16/Pt6 and Pro16/Pro6
    % does gamma vary over the same station at different mach?
    % Fix
    % Quantify effect on thrust
    % Rederive everything without the book and implement?
    % Look into realistic ranges of Mach numbers
    % Talk w/ propulsion reference?
    % Cry a bit then move on to a more rudimentary model?

    
M = linspace(0,1.5);
gamma = 1.4;
for ii = 1:100
    [pi_14] = pressure(M(ii),1.4);
    Pt_P_14(ii) = 1/pi_14;
    [pi_13] = pressure(M(ii),1.3);
    Pt_P_13(ii) = 1/pi_13;
end
plot(M,Pt_P_14,M,Pt_P_13)
xlabel('Mach Number')
ylabel('P / Pt')
legend('gamma = 1.4','gamma = 1.3')

%% Analyze M6 vs M16 at each P6/P16 (fixed A16_6)
A16_6 = 0.2715;

pi = [1.0042,1.0492,1.0706];
pi = [1.0492];

M_mixer = linspace(.01,1);

Pro16 = 8; %Assumption: Relative pressure at the fan. Has almost no effect on M16 vs M6, small noticeable (max of 5%) effect on Area Ratio
Pro6 = 300;
f = .03; %Assumption: fuel to air ratio at the mixer entry, very close to both cases
beta = .01; %Assumption, used in book results
ep1 = .05; %Assumption, used in book results
ep2 = .05; %Assumption, used in book results
R = 287; %J/kg k 

alpha_ref = .5; %For reset on iteration

%Setup state cell
state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
state(2:22,1) = {'0';'o0';'o2';'o16';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
[state] = derived_parameters(state,alpha_ref,beta,ep1,ep2,f);

%Setup state of bypass air
state(5,2) = {Pro16};
[state] = unFAIR3(state,5);
[~,~,To16,~,~,~,gamma16,~,~,R16,~] = state{5,:};

state(14,2) = {Pro6};
[state] = unFAIR3(state,14);
[~,~,To6,~,~,~,gamma6,~,~,R6,~] = state{14,:};

for ii = 1:size(pi,2)
    pi_bypass = pi(ii);
    Pro6 = Pro16*pi_bypass;
    %Setup state of core air
%     state(14,2) = {Pro6};
%     state(14,3) = {[]};
%     state(14,8) = {[]};
%     [state] = unFAIR3(state,14);
%     [~,~,~,~,~,cp6,gamma6,~,~,R6,~] = state{14,:};
    for jj = 1:100
        error =1; %error on alpha value, iterated due to small changes in 
        alpha = alpha_ref;
        while error > .00016
            gamma16 = 1.4;
            gamma6 = 1.3;
            
            %Iterate M6 and find MFP
            M6 = M_mixer(jj);
            MFP6 = MFP2(M6, gamma6, R);
            %Calculate M16 and MFP16
            [M16] = Kutta_mach(gamma16,M6,gamma6,pi_bypass);
            MFP16 = MFP2(M16, gamma16, R);
            %Find bypass ratios and quantify error
            [alpha_prime,alpha_i] = bypass_ratio(MFP16,Pro16,To16,MFP6,Pro6,To6,A16_6,beta,ep1,ep2,f);
            error = norm((alpha-alpha_i)/alpha);
            alpha = alpha_i;
            
            %Cell housekeeping
%             [state] = derived_parameters(state,alpha,beta,ep1,ep2,f);
%             state(5,2) = {[]};
%             state(5,8) = {[]};
%             [state] = unFAIR3(state,5);
%             [~,~,To16,~,~,cp16,gamma16,~,~,~,~] = state{5,:};
%             state(14,2) = {[]};
%             state(14,8) = {[]};
%             [state] = unFAIR3(state,14);
%             [~,~,~,~,~,cp6,gamma6,~,~,~,~] = state{14,:};
        end
        %Find pressure ratio
        state(15,2) = {[]};
        state(15,3) = {[]};
        state(15,8) = {[]};
        [pi_M_ideal,M6A,state] = mixer_ideal(state,alpha_prime,A16_6,M16,M6);
        %Store values
        M_bypass(ii,jj) = M16;
        M_post_mixer(ii,jj) = M6A;
        alpha_p(ii,jj) = alpha_prime;
        alpha_(ii,jj) = alpha;
        pi_M(ii,jj) = pi_M_ideal;
    end
end



figure  
plot(M_mixer,M_bypass(1,:),M_mixer,M_bypass(2,:),M_mixer,M_bypass(3,:),[.4,.3835, .3748],[.394,.4559,.4814],'x')
xlabel('Mixer Mach Number (M_6)')
ylabel('Bypass Mach Number (M_1_6)')
legend('reference','test1','test2')

figure
plot(M_mixer,M_post_mixer(1,:),M_mixer,M_post_mixer(2,:),M_mixer,M_post_mixer(3,:),[.4,.3835, .3748],[.4188,.4187,.4185],'x')
xlabel('Mixer Mach Number (M_6)')
ylabel('Post Mixer Mach Number (M_6_A)')
legend('reference','test1','test2')

% figure
% plot(M_mixer(30:70),alpha_p(1,30:70),M_mixer(30:70),alpha_p(2,30:70))
% xlabel('Mixer Mach Number')
% ylabel('Alpha Prime (\alpha')')
% legend('reference','test')

figure
plot(M_mixer(30:70),alpha_(1,30:70),M_mixer(30:70),alpha_(2,30:70),M_mixer(30:70),alpha_(3,30:70),[.4,.3835,.3748],[.449,.520,.576],'x')
xlabel('Mixer Mach Number (M_6)')
ylabel('Bypass Ratio (\alpha)')
legend('reference','test1','test2')

figure
plot(M_mixer(30:70),pi_M(1,30:70),M_mixer(30:70),pi_M(2,30:70))
xlabel('Mixer Mach Number (M_6)')
ylabel('\pi_i_d_e_a_l')
legend('reference','test')





A = M_bypass(1,:)./M_mixer;
B = M_bypass(2,:)./M_mixer;
C = M_bypass(3,:)./M_mixer;


figure  
plot(M_mixer,A,M_mixer,B,M_mixer,C)
xlabel('Mixer Mach Number (M_6)')
ylabel('Bypass Mach Number (M_1_6 / M_6)')
legend('reference','test1','test2')


%% Area
clear
clc
close all

A = [.1,.2,.4,.8,1.6];
pi = [1/.90,1/.92,1/.94,1/.96];

M_mixer = linspace(.01,1);

Pro16 = 8; %Assumption: Relative pressure at the fan. Has almost no effect on M16 vs M6, small noticeable (max of 5%) effect on Area Ratio
Pro6 = 300;
f = .03; %Assumption: fuel to air ratio at the mixer entry, very close to both cases
beta = .01; %Assumption, used in book results
ep1 = .05; %Assumption, used in book results
ep2 = .05; %Assumption, used in book results
R = 287; %J/kg k 

alpha_ref = .5; %For reset on iteration

%Setup state cell
state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
state(2:22,1) = {'0';'o0';'o2';'o16';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
[state] = derived_parameters(state,alpha_ref,beta,ep1,ep2,f);

% Setup state of bypass air
state(5,2) = {Pro16};
[state] = unFAIR3(state,5);
[~,~,To16,~,~,cp16,gamma16,~,~,R16,~] = state{5,:};

state(14,2) = {Pro6};
[state] = unFAIR3(state,14);
[~,~,To6,~,~,~,~,~,~,~,~] = state{14,:};

figure
hold on
xlabel('Mixer Mach Number (M_6)')
ylabel('Bypass Ratio (\alpha)')

for hh = 1:size(A,2)
A16_6 = A(hh);

for ii = 1:size(pi,2)
    pi_bypass = pi(ii);
    Pro6 = Pro16*pi_bypass;
    for jj = 1:100
        error =1; %error on alpha value, iterated due to small changes in 
        alpha = alpha_ref;
        while error > .0001
            gamma16 = 1.4;
            gamma6 = 1.3;
            
            %Iterate M6 and find MFP
            M6 = M_mixer(jj);
            MFP6 = MFP2(M6, gamma6, R);
            %Calculate M16 and MFP16
            [M16] = Kutta_mach(gamma16,M6,gamma6,pi_bypass);
            MFP16 = MFP2(M16, gamma16, R);
            %Find bypass ratios and quantify error
            [alpha_prime,alpha_i] = bypass_ratio(MFP16,Pro16,To16,MFP6,Pro6,To6,A16_6,beta,ep1,ep2,f);
            error = norm((alpha-alpha_i)/alpha);
            alpha = alpha_i;
        end
        %Find pressure ratio
        state(15,2) = {[]};
        state(15,3) = {[]};
        state(15,8) = {[]};
        [pi_M_ideal,M6A,state] = mixer_ideal(state,alpha_prime,A16_6,M16,M6);
        %Store values
        M_bypass(ii,jj) = M16;
        M_post_mixer(ii,jj) = M6A;
        alpha_p(ii,jj) = alpha_prime;
        alpha_(ii,jj) = alpha;
        pi_M(ii,jj) = pi_M_ideal;
    end
    plot(M_mixer(20:80),alpha_(ii,20:80))
end

% legend('1945-1965','1965-1985','1985-2005','2005-2025')


end

%% Edits:
% Fix pi mixer?

% Find estimate for alpha
    % Estimate P16/P6 as pi_burner
    % Estimate M6 as .5?

% Turn into function

%% Archived Code

% %% Analyze M6 vs M16 at each P6/P16 (fixed alpha)
% pi_b = [.9,.92,.94,.96,1]; %pi_burner
% pi = 1./pi_b; %P16/P6 whats a realistic range here?
% M_mixer = linspace(0,1);
% 
% Pro16 = 8; %Assumption: Relative pressure at the fan. Has almost no effect on M16 vs M6, small noticeable (max of 5%) effect on Area Ratio
% f = .03; %Assumption: fuel to air ratio at the mixer entry
% alpha = .449; %Start Assumption
% beta = .01; %Assumption
% ep1 = .05; %Assumption
% ep2 = .05; %Assumption
% 
% 
% state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
% state(2:22,1) = {'0';'o0';'o2';'o16';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
% [state] = derived_parameters(state,alpha,beta,ep1,ep2,f);
% 
% state(5,2) = {Pro16};
% state(5,4) = {0};
% [state] = unFAIR3(state,5);
% [~,~,To16,~,~,cp16,gamma16,~,~,R16,~] = state{5,:};
% 
% mdot6 = state{13,5};
% mdot16 = state{5,5};
% alpha_prime = mdot16/mdot6;
% 
% for ii = 1:size(pi,2)
%     pi_local = pi(ii);
%     Pro6 = Pro16/pi_local;
%     state(14,2) = {Pro6};
%     state(14,3) = {[]};
%     state(14,8) = {[]};
%     [state] = unFAIR3(state,14);
%     [~,~,To6,~,~,cp6,gamma6,~,~,R6,~] = state{14,:};
% 
%     for jj = 1:100
%         M6 = M_mixer(jj);
%         [M16] = Kutta_mach(gamma16,M6,gamma6,pi_local);
%         MFP16 = MFP2(M16, gamma16, R16);
%         MFP6 = MFP2(M6, gamma6, R6);
%         [A16_6] = area(To16, Pro16, mdot16,MFP16, To6, Pro6, mdot6,MFP6);
%         [alpha_prime] = bypass(MFP16,Pro16,To16,MFP6,Pro6,To6,A16_6);
%         [alpha_prime_approx] = M16/M6*A16_6/pi_local^2;
%         M_bypass(ii,jj) = M16;
%         A_ratio(ii,jj) = A16_6;
%         alpha_p(ii,jj) = alpha_prime;
%         alpha_p_approx(ii,jj) = alpha_prime_approx;
% 
%     end
% end
% 
% figure
% plot(M_mixer,M_bypass(1,:),M_mixer,M_bypass(2,:),M_mixer,M_bypass(3,:),M_mixer,M_bypass(4,:),M_mixer,M_bypass(5,:),[.4,.3835],[.394,.4559],'x')
% xlabel('Mixer Mach Number')
% ylabel('Bypass Mach Number')
% legend('1945-1965','1965-1985','1985-2005','2005-2025','pi = 1')
% 
% 
% figure
% plot(M_mixer(30:70),A_ratio(1,30:70),M_mixer(30:70),A_ratio(2,30:70),M_mixer(30:70),A_ratio(3,30:70),M_mixer(30:70),A_ratio(4,30:70),M_mixer(30:70),A_ratio(5,30:70))
% xlabel('Mixer Mach Number')
% ylabel('Area Ratio 16/6')
% legend('1945-1965','1965-1985','1985-2005','2005-2025','pi = 1')
% 
% figure
% plot(M_mixer(30:70),alpha_p(1,30:70),M_mixer(30:70),alpha_p(2,30:70),M_mixer(30:70),alpha_p(3,30:70),M_mixer(30:70),alpha_p(4,30:70),M_mixer(30:70),alpha_p(5,30:70))
% hold on
% plot(M_mixer(30:70),alpha_p_approx(1,30:70),M_mixer(30:70),alpha_p_approx(2,30:70),M_mixer(30:70),alpha_p_approx(3,30:70),M_mixer(30:70),alpha_p_approx(4,30:70),M_mixer(30:70),alpha_p_approx(5,30:70))
% 
% xlabel('Mixer Mach Number')
% ylabel('Alpha Prime')
% % legend('1945-1965','1965-1985','1985-2005','2005-2025','pi = 1')













%% Functions

function [A16_6] = area(Tt16, Pt16, mdot16,MFP16, Tt6, Pt6, mdot6,MFP6)
%Calculates area ratio for bypass and core of mixer
Ar16 = mdot16*sqrt(Tt16)/Pt16    *   MFP16; %Relative pressure means relative area
Ar6 = mdot6*sqrt(Tt6)/Pt6    *   MFP6;
A16_6 =Ar16/Ar6;
end

function [M16] = Kutta_mach(gamma16,M6,gamma6,pi6_16)
%calculates the mach if kutta condition is satisfied between the two states
P_Pt6 = pressure(M6,gamma6);
P_Pt16 = P_Pt6*pi6_16;
M16 = sqrt((2/(gamma16 - 1)) * (P_Pt16 ^ ((gamma16 - 1)/gamma16)  - 1));
end

function [alpha_prime,alpha] = bypass_ratio(MFP16,Pt16,Tt16,MFP6,Pt6,Tt6,A16_6,beta,ep1,ep2,f)
%Calculates bypass ratio at mixer
alpha_prime = Pt16/Pt6 * MFP16/MFP6 * sqrt(Tt6/Tt16) * A16_6;
% alpha_prime = MFP16/MFP6 * A16_6 * 1.65;
alpha = alpha_prime * ((1-beta-ep1-ep2)*(1+f) + ep1 +ep2);
end

function [pi,M6A,state] = mixer_ideal(state,alpha_prime,A16_6,M16,M6)
%Calculates state and ideal pressure ratio of the mixer
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
%Calculates mass flow, stagnation temp and pressure, and Area dependent MFP
MFP = mdot*sqrt(Tt) / (Pt*A); %kg/s*K^.5 / Pa*m^2 = s*K^.5 / m
end

function [MFP] = MFP2(M, gamma, R)
%Calculates Mach, gamma, and R dependent MFP
[P_Pt] = pressure(M,gamma);
[T_Tt] = temperature(M,gamma);
MFP = M*sqrt(gamma/R)/sqrt(T_Tt)*sqrt(P_Pt); %sqrt(s^2*K / m^2) = s*K^.5 / m
end

function [P_Pt] = pressure(M,gamma)
%Calculates static over stagnation pressure for isentropic compressible
%gas
Pt_P = (1 + (gamma - 1)/2*M^2)^(gamma/(gamma-1));
P_Pt = 1/Pt_P;
end

function [T_Tt] = temperature(M,gamma)
%Calculates static over stagnation tempurature for isentropic compressible
%gas
T_Tt = (1 - (gamma - 1)/2*M^2)^-1;
end

function [state] = derived_parameters(state,alpha,beta,ep1,ep2,f)
% Calculates mass flow and fuel to air ratio at every station given inputs
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