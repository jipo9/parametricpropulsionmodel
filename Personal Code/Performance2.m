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
alt = 0 / 3.281; %altitude [m from feet]
M0 = .1;
pi_cL = 5; %Mid approx for turbofan, pi_cL = pi_cH %5.9
pi_cH = 6; %Mid approx for turbofan, pi_cL = pi_cH
alpha = 6.9;
[pi_f] = optimalfan(alpha); %approximation for turbofan based on bypass ratio
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





% alt = 10000 / 3.281; %altitude [m from feet]
% M0 = .338;
% pi_f = 1.5; %Mid approx for turbofan %2
% pi_cL = 1.5; %Mid approx for turbofan, pi_cL = pi_cH %5.9
% pi_cH = 10; %Mid approx for turbofan, pi_cL = pi_cH
% alpha = 6.9;
% beta = 0;
% PtoH = 0;
% PtoL = 0;
% A0 = pi*(.469/2)^2;
% year = 2011;



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
alt = [0,12500,25000]./3.281; %Operational envelop in altitude
M0 = linspace(.1,.45,15); %Operational envelop in Mach
[T_std, ~, P_std, ~] = atmosisa(0); %obtain standard atmospheric conditions at SL

% Perform off design analysis over operational scope
for i = 1:length(alt)
    for j = 1:length(M0)
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
    display(['Analysis is ',num2str(100*i/length(alt)),'% complete'])
end

%% DGEN DATA
dgen390.M(:,1) = linspace(.1,.45,15);
dgen390.F(:,1) = [638.87;616.06;594.09;572.88;552.39;532.56;513.31;494.59;476.35;458.52;441.05;423.89;406.97;390.25;373.67];
dgen390.F(:,2) = [460.82;444.91;429.71;415.16;401.2;387.78;374.83;362.3;350.13;338.27;326.67;315.26;303.99;292.77;279.78];
dgen390.F(:,3) = [309.99;299.64;289.89;280.27;271.07;262.25;253.76;245.55;237.56;229.74;222.03;214.39;206.75;199.08;191.3];

dgen390.S(:,1) = [.53;.55;.57;.59;.61;.63;.66;.68;.7;.73;.75;.78;.81;.84;.87];
dgen390.S(:,2) = [.52;.54;.56;.58;.6;.62;.64;.66;.68;.7;.73;.75;.77;.8;.83];
dgen390.S(:,3) = [.52;.54;.55;.57;.59;.61;.62;.64;.66;.68;.7;.72;.74;.76;.78];

dgen390.Fnondim(:,1) = dgen390.F(:,1)./dgen390.F(1,1);
dgen390.Fnondim(:,2) = dgen390.F(:,2)./dgen390.F(1,1);
dgen390.Fnondim(:,3) = dgen390.F(:,3)./dgen390.F(1,1);

dgen390.Snondim(:,1) = dgen390.S(:,1)./dgen390.S(1,1);
dgen390.Snondim(:,2) = dgen390.S(:,2)./dgen390.S(1,1);
dgen390.Snondim(:,3) = dgen390.S(:,3)./dgen390.S(1,1);

Fnondim(:,1) = F(1,:)./F(1,1);
Fnondim(:,2) = F(2,:)./F(1,1);
Fnondim(:,3) = F(3,:)./F(1,1);

Snondim(:,1) = S(1,:)./S(1,1);
Snondim(:,2) = S(2,:)./S(1,1);
Snondim(:,3) = S(3,:)./S(1,1);




%% Plots
% Execution of plots. See figures for more details

for plots = 1
h1 = figure(1);
h1.WindowStyle = 'docked'; 
for i = 1:length(alt)
plot(M0,F(i,:),'linewidth',1.5)
hold on
title('Thrust vs. Mach Number')
xlabel('Mach Number')
ylabel('Thrust (lb)')
grid('on')
end
legend('SL','12.5k','25k')


h2 = figure(2);
h2.WindowStyle = 'docked';
for i = 1:length(alt)
plot(M0,S(i,:),'linewidth',1.5)
title('SFC vs. Mach Number')
hold on
xlabel('Mach Number')
ylabel('SFC')
grid('on')
end
legend('SL','12.5k','25k')


% h3 = figure(3);
% h3.WindowStyle = 'docked';
% for i = 1:length(alt)
% plot(mdot44_cor(i,:),pitH(i,:),'linewidth',1.5)
% hold on
% title('HP Turb Pressure Ratio vs. Corr Mass Flow')
% xlabel('Corrected Mass Flow')
% ylabel('HP Turbine Pressure Ratio')
% grid('on')
% end
% legend('SL','12.5k','25k')
% 
% 
% h4 = figure(4);
% h4.WindowStyle = 'docked';
% for i = 1:length(alt)
% plot(mdot25_cor(i,:),picL(i,:),'linewidth',1.5)
% hold on
% title('LP Comp Pressure Ratio vs. Corr Mass Flow')
% xlabel('Corrected Mass Flow')
% ylabel('LP Compressor Ratio')
% grid('on')
% end
% legend('SL','12.5k','25k')


h5 = figure(5);
h5.WindowStyle = 'docked'; 
for i = 1:length(alt)
subplot(3,1,1)
plot(M0,pinlet(i,:),'linewidth',1.5)
hold on
title('Mach Number vs. Inlet Air Momentum')
xlabel('Mach Number')
ylabel('Inlet Air Momentum')
grid('on')
end
legend('SL','12.5k','25k')
for i = 1:length(alt)
subplot(3,1,2)
plot(M0,pcore(i,:),'linewidth',1.5)
hold on
title('Mach Number vs. Core Air Momentum')
xlabel('Mach Number')
ylabel('Core Air Momentum')
grid('on')
end
legend('SL','12.5k','25k')
for i = 1:length(alt)
subplot(3,1,3)
plot(M0,pbypass(i,:),'linewidth',1.5)
hold on
title('Mach Number vs. Bypass Air Momentum')
xlabel('Mach Number')
ylabel('Bypass Air Momentum')
grid('on')
end
legend('SL','12.5k','25k')

% h6 = figure(6);
% h6.WindowStyle = 'docked'; 
% for i = 1:length(alt)
% subplot(3,1,1)
% plot(M0,mdot0(i,:),'linewidth',1.5)
% hold on
% title('Mach Number vs. Inlet Mass Flow')
% xlabel('Mach Number')
% ylabel('Inlet Air Mass Flow')
% grid('on')
% end
% legend('SL','12.5k','25k')
% for i = 1:length(alt)
% subplot(3,1,2)
% plot(M0,mdot9(i,:),'linewidth',1.5)
% hold on
% title('Mach Number vs. Core Mass Flow')
% xlabel('Mach Number')
% ylabel('Core Air Mass Flow')
% grid('on')
% end
% legend('SL','12.5k','25k')
% for i = 1:length(alt)
% subplot(3,1,3)
% plot(M0,mdot19(i,:),'linewidth',1.5)
% hold on
% title('Mach Number vs. Bypass Mass Flow')
% xlabel('Mach Number')
% ylabel('Bypass Air Mass Flow')
% grid('on')
% end
% legend('SL','12.5k','25k')
% 
% h7 = figure(7);
% h7.WindowStyle = 'docked'; 
% for i = 1:length(alt)
% subplot(3,1,1)
% plot(M0,v0(i,:),'linewidth',1.5)
% hold on
% title('Mach Number vs. Inlet Velocity')
% xlabel('Mach Number')
% ylabel('Inlet Velocity')
% grid('on')
% end
% legend('SL','12.5k','25k')
% for i = 1:length(alt)
% subplot(3,1,2)
% plot(M0,v9(i,:),'linewidth',1.5)
% hold on
% title('Mach Number vs. Core Velocity')
% xlabel('Mach Number')
% ylabel('Core Velocity')
% grid('on')
% end
% legend('SL','12.5k','25k')
% for i = 1:length(alt)
% subplot(3,1,3)
% plot(M0,v19(i,:),'linewidth',1.5)
% hold on
% title('Mach Number vs. Bypass Velocity')
% xlabel('Mach Number')
% ylabel('Bypass Velocity')
% grid('on')
% end
% legend('SL','12.5k','25k')

h8 = figure(8);
h8.WindowStyle = 'docked';
for i = 1:length(alt)
plot(M0(:),alpha(i,:),'linewidth',1.5)
hold on
title('Bypass Ratio vs. Mach Number')
xlabel('Mach Number')
ylabel('\alpha')
grid('on')
end
legend('SL','12.5k','25k')
end

h9 = figure(9);
h9.WindowStyle = 'docked';
jj = size(dgen390.F);
for i = 1:jj(2)
plot(dgen390.M,dgen390.Fnondim(:,i),'linewidth',1.5)
hold on
end
plot(dgen390.M,Fnondim(:,1),'--','linewidth',1.5,'color','#0072BD')
hold on
plot(dgen390.M,Fnondim(:,2),'--','linewidth',1.5,'color','#D95319')
hold on
plot(dgen390.M,Fnondim(:,3),'--','linewidth',1.5,'color','#EDB120')
legend('DGEN SL','DGEN 12.5k Feet','DGEN 25k Feet','Off Design SL','Off Design 12.5k Feet','Off Design 25k Feet')
title('Normalized Thrust vs. Mach Number')
xlabel('Mach Number')
ylabel('Normalized Thrust')
grid on

h10 = figure(10);
h10.WindowStyle = 'docked';
jj = size(dgen390.S);
for i = 1:jj(2)
plot(dgen390.M,dgen390.Snondim(:,i),'linewidth',1.5)
hold on
end
plot(dgen390.M,Snondim(:,1),'--','linewidth',1.5,'color','#0072BD')
hold on
plot(dgen390.M,Snondim(:,2),'--','linewidth',1.5,'color','#D95319')
hold on
plot(dgen390.M,Snondim(:,3),'--','linewidth',1.5,'color','#EDB120')
legend('DGEN SL','DGEN 12.5k Feet','DGEN 25k Feet','Off Design SL','Off Design 12.5k Feet','Off Design 25k Feet')
title('Normalized SFC vs. Mach Number')
xlabel('Mach Number')
ylabel('Normalized SFC')
grid on


