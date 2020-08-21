clear
clc
clf
close all

%% Read Me
% The following function creates performance charts of a subsonic turbofan
% no afterburner and seperate exhausts

% Summary:
    % Inputs - Enter 9 engine parameters at engine on design case (usually either takeoff of max cruise)
    % On-Design - Calculates component efficiencies and performance at specified design case
    % Off-Design - Iterates through operational envelope to find max thrust at any feasible condition
    % Plots - Returns figures detailing engine performance including thrust, specific fuel consumption, inlet and exhaust momentum, etc.

% All relavent information is stored in 1 of 4 cells:
    % Component parameters, both given and found, are stored in
    % "component". These parameters include pressure ratio, polytropic
    % efficiency, enthalpy ratio, and mechanical efficiency

    % Inputs general beyond components is stored in "design"

    % Thermodynamic states of gas in the engine are modeled at each state
    % through the engine and live in "state"

    % Finally, overall performance parameters are given in "performance"
    % for comparisson

%% To do
% Get sensitivity study back?
% Get a better way of estimating inputs? Limit on pi_f from input design
% point and alpha???
%% Inputs

% Engine : DGEN-390 Price Induction
% Design point is 10kft and M = .338
alt = 10000 / 3.281; %altitude [m from feet]
M0 = .338;
pi_f = 1.5; %Mid approx for turbofan %2
pi_cL = 1.5; %Mid approx for turbofan, pi_cL = pi_cH %5.9
pi_cH = 10; %Mid approx for turbofan, pi_cL = pi_cH
alpha = 6.9;
beta = 0;
PtoH = 0;
PtoL = 0;
A0 = pi*(.469/2)^2;
year = 2011;

% Engine : JT9D
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

%% On design
% Run on-design
[state,component,design,inputs,performance] = on_design(alt,M0,pi_f,pi_cL,pi_cH,alpha,beta,PtoH,PtoL,A0,year);

% Display thrust and SFC
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
%% Off design
% Setup off design
altR = alt;
mdotc_R = CorrectedMassFlow(state,altR,altR,component);
alt = [0,10000,20000,30000]./3.281; %Operational envelop in altitude
M0 = linspace(.1,.8,10); %Operational envelop in Mach
[T_std, ~, P_std, ~] = atmosisa(0); %obtain standard atmospheric conditions at SL

% Perform off design analysis over operational scope
for i = 1:length(alt)
    for j = 1:length(M0)
        j
        %Perform off-design analysis
        [state_i,component_i,design_i,inputs_i,performance_i] = off_design2(state,component,design,inputs,M0(j),alt(i),A0);            
        
        %Find component pressure ratios and corrected mass flow
        [~,pi_r,pi_d,pi_f,pi_cL,pi_cH,pi_b,~,pi_tH,~,~,pi_tL,~,~,~,pi_n,~] = component_i{:,2};
        [~, ~, P0, ~] = atmosisa(alt(i)); 
        Po25_Std = pi_r*pi_d*pi_cL*P0/P_std;
        To25 = state_i{6,3};
        mdot25 = state_i{6,5};
        mdot25_cor(i,j) = mdot25 * sqrt(To25/T_std) / (Po25_Std);
        picL(i,j) = pi_cL;
        
        Po44_Std = pi_r*pi_d*pi_cL*pi_cH*pi_b*pi_tH*P0/P_std;
        To44 = state_i{11,3};
        mdot44 = state_i{11,5};
        mdot44_cor(i,j) = mdot44 * sqrt(To44/T_std) / (Po44_Std);
        pitH(i,j) = pi_tH;
       
        % Store performance parameters
        S(i,j) = performance_i{2,2} / ((.453592/3600)/4.44822);
        F(i,j) = performance_i{2,1} * 0.224809;
        pinlet(i,j) = performance_i{2,8};
        pcore(i,j) = performance_i{2,9};
        pbypass(i,j) = performance_i{2,10};
        mdot0(i,j) = performance_i{2,11};
        mdot9(i,j) = performance_i{2,12};
        mdot19(i,j) = performance_i{2,13};
        v0(i,j) = performance_i{2,14};
        v9(i,j) = performance_i{2,15};
        v19(i,j) = performance_i{2,16};
        
        %Store effective bypass ratio
        alpha(i,j) = design_i{2,2};
    end
    a= 1;
end

%% Plots
% Execution of plots. See figures for more details

for plots = 1
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
plot(mdot44_cor(i,:),pitH(i,:),'linewidth',1.5)
hold on
title('HP Turb Pressure Ratio vs. Corr Mass Flow')
legend('SL','10k','20k','30k')
xlabel('Corrected Mass Flow')
ylabel('HP Turbine Pressure Ratio')
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


h5 = figure(5);
h5.WindowStyle = 'docked'; 
for i = 1:length(alt)
subplot(3,1,1)
plot(M0,pinlet(i,:),'linewidth',1.5)
hold on
title('Mach Number vs. Inlet Air Momentum')
legend('SL','10k','20k','30k')
xlabel('Mach Number')
ylabel('Inlet Air Momentum')
grid('on')
end
for i = 1:length(alt)
subplot(3,1,2)
plot(M0,pcore(i,:),'linewidth',1.5)
hold on
title('Mach Number vs. Core Air Momentum')
legend('SL','10k','20k','30k')
xlabel('Mach Number')
ylabel('Core Air Momentum')
grid('on')
end
for i = 1:length(alt)
subplot(3,1,3)
plot(M0,pbypass(i,:),'linewidth',1.5)
hold on
title('Mach Number vs. Bypass Air Momentum')
legend('SL','10k','20k','30k')
xlabel('Mach Number')
ylabel('Bypass Air Momentum')
grid('on')
end

h6 = figure(6);
h6.WindowStyle = 'docked'; 
for i = 1:length(alt)
subplot(3,1,1)
plot(M0,mdot0(i,:),'linewidth',1.5)
hold on
title('Mach Number vs. Inlet Mass Flow')
legend('SL','10k','20k','30k')
xlabel('Mach Number')
ylabel('Inlet Air Mass Flow')
grid('on')
end
for i = 1:length(alt)
subplot(3,1,2)
plot(M0,mdot9(i,:),'linewidth',1.5)
hold on
title('Mach Number vs. Core Mass Flow')
legend('SL','10k','20k','30k')
xlabel('Mach Number')
ylabel('Core Air Mass Flow')
grid('on')
end
for i = 1:length(alt)
subplot(3,1,3)
plot(M0,mdot19(i,:),'linewidth',1.5)
hold on
title('Mach Number vs. Bypass Mass Flow')
legend('SL','10k','20k','30k')
xlabel('Mach Number')
ylabel('Bypass Air Mass Flow')
grid('on')
end

h7 = figure(7);
h7.WindowStyle = 'docked'; 
for i = 1:length(alt)
subplot(3,1,1)
plot(M0,v0(i,:),'linewidth',1.5)
hold on
title('Mach Number vs. Inlet Velocity')
legend('SL','10k','20k','30k')
xlabel('Mach Number')
ylabel('Inlet Velocity')
grid('on')
end
for i = 1:length(alt)
subplot(3,1,2)
plot(M0,v9(i,:),'linewidth',1.5)
hold on
title('Mach Number vs. Core Velocity')
legend('SL','10k','20k','30k')
xlabel('Mach Number')
ylabel('Core Velocity')
grid('on')
end
for i = 1:length(alt)
subplot(3,1,3)
plot(M0,v19(i,:),'linewidth',1.5)
hold on
title('Mach Number vs. Bypass Velocity')
legend('SL','10k','20k','30k')
xlabel('Mach Number')
ylabel('Bypass Velocity')
grid('on')
end

h8 = figure(8);
h8.WindowStyle = 'docked';
for i = 1:length(alt)
plot(M0(:),alpha(i,:),'linewidth',1.5)
hold on
title('Bypass Ratio vs. Mach Number')
legend('SL','10k','20k','30k')
xlabel('Mach Number')
ylabel('\alpha')
grid('on')
end
end