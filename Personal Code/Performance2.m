clear
clc
clf
close all

%% Read Me
% The following function creates performance charts of a subsonic turbofan
% no afterburner and no mixer

%% On design

% Engine : DGEN-390 Price Induction
% Design point is 10kft and M = .338
alt = 10000 / 3.281; %altitude [m from feet]
M0 = .338;
pi_f = 1.4; %Mid approx for turbofan %2
pi_cL = 5.9; %Mid approx for turbofan, pi_cL = pi_cH %5.9
pi_cH = 5.9; %Mid approx for turbofan, pi_cL = pi_cH
alpha = 6.9;
beta = 0;
PtoH = 0;
PtoL = 0;
A0 = pi*(.469/2)^2;
year = 2011;

% Engine : JT9D
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

% Run on-design
[state,component,design,inputs,performance] = on_design(alt,M0,pi_f,pi_cL,pi_cH,alpha,beta,PtoH,PtoL,A0,year);

% Calculate thrust and SFC
F = performance{2,1} * .224809;
mdot0 = state{2,5};
F_mdot = F/mdot0 * (1/9.806655);
S = performance{2,2} / ((.453592/3600)/4.44822);

fprintf('%s\n','--------On Design---------')
fprintf('%s%.2f\n','Calculated Thrust = ',F)
fprintf('%s%.2f\n\n','Calculated SFC = ',S)
fprintf('%s\n\n','For the DGEN390, We want F of 354 and S of .72 at design point')
fprintf('%s\n','Range for F ~ 300-570')
fprintf('%s\n','Tange for S ~ .45-.83')
fprintf('%s\n','--------------------------')

M9 = performance{2,7};
[A45_9] = A_turb(component, M9);

%% Off design


alt = [0,10000,20000,30000]./3.281;
M0 = linspace(.1,.45,20);


[T_std, ~, P_std, ~] = atmosisa(0); %obtain standard atmospheric conditions at SL
for i = 1:length(alt)
    for j = 1:length(M0)
        [state,component,design,inputs,performance] = off_design(state,component,design,inputs,M0(j),alt(i),A0,A45_9);      
        [~,pi_r,pi_d,pi_f,pi_cL,pi_cH,pi_b,~,pi_tH,~,~,pi_tL,~,~,~,pi_n,~] = component{:,2};
        [~, ~, P0, ~] = atmosisa(alt(i)); %obtain standard atmospheric conditions
        Po25_Std = pi_r*pi_d*pi_cL*P0/P_std;
        [state,component,design,inputs,performance] = off_design(state,component,design,inputs,M0(j),alt(i),A0,A45_9);
        
        %finding component pressure ratios and corrected mass flow
        [~,pi_r,pi_d,pi_f,pi_cL,pi_cH,pi_b,~,pi_tH,~,~,pi_tL,~,~,~,pi_n,~] = component{:,2};
        [~, ~, P0, ~] = atmosisa(alt(i)); 
        Po25_Std = pi_r*pi_d*pi_cL*P0/P_std;
        To25 = state{6,3};
        mdot25 = state{6,5};
        mdot25_cor(i,j) = mdot25 * sqrt(To25/T_std) / (Po25_Std);
        pif(i,j) = pi_f;
        picL(i,j) = pi_cL;
       
        S(i,j) = performance{2,2} / ((.453592/3600)/4.44822);
        F(i,j) = performance{2,1} * 0.224809;
    end
end

%% Plots

h1 = figure(1);
h1.WindowStyle = 'docked'; 
for i = 1:length(alt)
plot(M0,F(i,:),'linewidth',1.5)
hold on
title('Thrust vs. Mach Number')
legend('SL','10k','20k','30k')
xlabel('Mach Number')
ylabel('Thrust (lb)')
grid('on')
end


h2 = figure(2);
h2.WindowStyle = 'docked';
for i = 1:length(alt)
plot(M0,S(i,:),'linewidth',1.5)
title('SFC vs. Mach Number')
hold on
legend('SL','10k','20k','30k')
xlabel('Mach Number')
ylabel('SFC')
grid('on')
end


h3 = figure(3);
h3.WindowStyle = 'docked';
for i = 1:length(alt)
plot(M0,pif(i,:),'linewidth',1.5)
title('Fan Pressure Ratio vs. Mach Number')
hold on
legend('SL','10k','20k','30k','40k')
xlabel('Mach Number')
ylabel('Fan Pressure Ratio')
grid('on')
end


h4 = figure(4);
h4.WindowStyle = 'docked';
for i = 1:length(alt)
plot(mdot25_cor(i,:),picL(i,:),'linewidth',1.5)
hold on
title('LP Comp Pressure Ratio vs. Corr Mass Flow')
legend('SL','10k','20k','30k')
xlabel('Corrected Mass Flow')
ylabel('LP Compressor Ratio')
grid('on')
end


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