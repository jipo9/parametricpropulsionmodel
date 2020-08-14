clear
clc
clf
close all

%% To do:

% Establish Driving Equaitons/ Functions for new turbofan model
% Adjust on design inputs
% Put on design and off design into funcitons

% Assume 1 spool?
% How do we use static/takeoff thrust?
%% Read Me
% The following function creates performance charts of a subsonic turbofan
% w/ no afterburner or mixer

%% On design
% Input data from "On-Design Case". This case should be the ideal case that
% the engine is designed to perform under


%Engine being utilized is the DGEN-390 Price Induction
alt = 10000 / 3.281; %altitude [m from feet]
M0 = .338;
pi_f = 1.2; %Mid approx for turbofan %2
pi_cL = 6; %Mid approx for turbofan, pi_cL = pi_cH %5.9
pi_cH = 6; %Mid approx for turbofan, pi_cL = pi_cH
alpha = 6.9;
beta = 0;
PtoH = 0;
PtoL = 0;
A0 = pi*(.469/2)^2;
year = 2011;

%Engine JT9D
% To4 == 1970F

% alt = 30000 / 3.281;
% M0 = .8;
% pi_f = 22.6/14.7;
% pi_cL = 32.1/14.7;
% pi_cH = 316/32.1;
% alpha = 5;
% beta = 0;
% PtoH = 0;
% PtoL = 0;
% A0 = pi*(2.34/2)^2;
% year = 1966;

[state,component,design,inputs,performance] = on_design(alt,M0,pi_f,pi_cL,pi_cH,alpha,beta,PtoH,PtoL,A0,year);

M9 = performance{2,7};
[A45_9] = A_turb(component, M9);

F = performance{2,1} * .224809
mdot0 = state{2,5};
F_mdot = F/mdot0 * (1/9.806655);
S = performance{2,2} / ((.453592/3600)/4.44822)
disp('We want F of 291   Range for S is .45 to ,8')


%% Off design



[state2,component2,design2,inputs2,performance2] = off_design(state,component,design,inputs,component,M0,alt,A0,A45_9);



componentR = component;

alt = [0,10000,20000,30000]./3.281;
M0 = linspace(.1,.45,20);


[T_std, ~, P_std, ~] = atmosisa(0); %obtain standard atmospheric conditions
for i = 1:length(alt)
    for j = 1:length(M0)
        %[state,component,design,inputs,performance] = on_design(alt(i),M0(j),pi_f,pi_cL,pi_cH,alpha,beta,PtoH,PtoL,A0,year);
        [state,component,design,inputs,performance] = off_design(state,component,design,inputs,componentR,M0(j),alt(i),A0,A45_9);
        [~,pi_r,pi_d,pi_f,pi_cL,pi_cH,pi_b,~,pi_tH,~,~,pi_tL,~,~,~,pi_n,~] = component{:,2};
        [~, ~, P0, ~] = atmosisa(alt(i)); %obtain standard atmospheric conditions
        
        Po25_Std = pi_r*pi_d(2)*pi_cL*P0/P_std;
        To25 = state{6,3};
        mdot25 = state{6,5};
        mdot25_cor(i,j) = mdot25 * sqrt(To25/T_std) / (Po25_Std);
%         pif(i,j) = pi_f;
        picL(i,j) = pi_cL;
       
        S(i,j) = performance{2,2} / ((.453592/3600)/4.44822);
        T = performance{2,1};
        F(i,j) = T;
    end
end

h1 = figure(1);
h1.WindowStyle = 'docked'; 
for i = 1:length(alt)
plot(M0,F(i,:)*0.224809,'linewidth',1.5)
hold on
title('Method 1 - Holding Turbine Values Constant')
legend('SL','10k','20k','30k')
xlabel('Mach Number')
ylabel('Thrust (lb)')
grid('on')
end

h2 = figure(2);
h2.WindowStyle = 'docked';
for i = 1:length(alt)
plot(M0,S(i,:),'linewidth',1.5)
title('Method 1 - Holding Turbine Values Constant')
hold on
legend('SL','10k','20k','30k')
xlabel('Mach Number')
ylabel('SFC')
grid('on')
end


% h3 = figure(3);
% h3.WindowStyle = 'docked';
% for i = 1:length(alt)
% plot(M0,pif(i,:),'linewidth',1.5)
% hold on
% legend('SL','10k','20k','30k','40k')
% xlabel('Mach Number')
% ylabel('Fan Pressure Ratio')
% grid('on')
% end


h4 = figure(4);
h4.WindowStyle = 'docked';
for i = 1:length(alt)
plot(mdot25_cor(i,:),picL(i,:),'linewidth',1.5)
hold on
title('Method 1 - Holding Turbine Values Constant')
legend('SL','10k','20k','30k')
xlabel('Corrected Mass Flow')
ylabel('LP Compressor Ratio')
grid('on')
end

% Loop for each case
    % Use modified off design analysis for 2 spool, no afterburner, no mixer
        % No iterative scheme, performance is funciton of flight condition
        % Add theta break in later (or now since its easy?

%% Eqns
% F = mdot16*v16 + mdot9*v9 - mdot0*v0; %Assume perfectly expanded
% S = f_0 * mdot0 / F;
% Overall_efficieny = v0/(hPR*S);
% Propulsive_efficiency = (F*v0) / (.5*mdot16*v16^2 + .5*mdot9*v9^2 - .5*mdot0*v0^2 + PTOH + PTOL)
% Thermal_efficiency = Overall_efficieny/Propulsive_efficiency;
function [A45_6] = A_turb(component, M6)
    gamma = 1.3;
    [pi_tL,~,tau_tL] = component{12,2:4};
    A45_6 = pi_tL/sqrt(tau_tL) * M6 * 1/(2/(gamma+1)*(1+(gamma-1)/2*M6^2))^((gamma+1)/(2*(gamma-1)));
end