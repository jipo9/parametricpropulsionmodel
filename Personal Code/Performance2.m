clear
clc
% close all

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
alt = 20000 / 3.281; %altitude [m from feet]
M0 = .35;
pi_f = 2; %Mid approx for turbofan
pi_cL = 5.9; %Mid approx for turbofan, pi_cL = pi_cH
pi_cH = 5.9; %Mid approx for turbofan, pi_cL = pi_cH
alpha = 7.5;
beta = 0;
PtoH = 0;
PtoL = 0;
A0 = pi*(.469/2)^2;
year = 2011;

%Engine JT9D
% To4 == 1970F

alt = 30000;
M0 = .8;
pi_f = 22.6/14.7;
pi_cL = 32.1/14.7;
pi_cH = 316/32.1;
alpha = 5;
beta = 0;
PtoH = 0;
PtoL = 0;
A0 = pi*(2.34/2)^2;
year = 1966;

[state,component,design,inputs,performance] = on_design(alt,M0,pi_f,pi_cL,pi_cH,alpha,beta,PtoH,PtoL,A0,year);


F = performance{2,1}
mdot0 = state{2,5}
F_mdot = F/mdot0 * (1/9.806655)
S = performance{2,2} / ((.453592/3600)/4.44822)
disp('Range for F/mdot is 13-27 and range for S is .67 to 1.03')


%% Off design

componentR = component;

alt = [0,10000,20000,30000,40000]./3.281;
M0 = linspace(.1,.85,10);

for i = 1:length(alt)
    for j = 1:length(M0)
        [state,component,design,inputs,performance] = off_design(state,component,design,inputs,componentR,M0(j),alt(i),A0);
        S(i,j) = performance{2,2};
        T = performance{2,1};
        F(i,j) = T;
    end
end

h1 = figure(1);
h1.WindowStyle = 'docked'; 
plot(M0,F(3,:))
xlabel('Mach Number')
ylabel('Thrust (N)')
grid('on')

h2 = figure(2);
h2.WindowStyle = 'docked'; 
plot(M0,S(1,:))
xlabel('Mach Number')
ylabel('SFC')
grid('on')


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