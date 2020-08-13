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

%Engine being ustilized is the DGEN-390 Price Induction
alt = 5000 / 3.281; % Max altitude [m from feet]
M0 = .35; %Max Mach
pi_f = 2.7; %Mid approx for fan is 2.7
pi_cL = 5.9; %Mid approx for low pressure compressor is 5.9, pi_cL = pi_cH
pi_cH = 5.9; %Mid approx for high pressure compressor is 5.9, pi_cL = pi_cH
alpha = 7.5; % Bypass Ratio Given in engine data
beta = 0; % Assume no bleed air
PtoH = 0; % Assume no power takeoff on high pressure shaft
PtoL = 0; % Assume no power takeoff on high pressure shaft
A0 = pi*(.469/2)^2; % Assume inlet area (only really matters wjen power takeoff is given)
year = 2011;

[state,component,design, inputs, performance] = on_design(alt,M0,pi_f,pi_cL,pi_cH,alpha,beta,PtoH,PtoL,A0,year);


F = performance{2,1};
mdot0 = state{2,5};
F_mdot = F/mdot0 * (1/9.806655)
S = performance{2,2} / ((.453592/3600)/4.44822)
disp('Range for F/mdot is 13-27 and range for S is .67 to 1.03')


%% Off design plots
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