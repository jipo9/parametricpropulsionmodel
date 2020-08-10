clear
clc
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
S = 3.2251e-05; %kg/s/N
T = 5.5923e+04; % N Cruise
mdot0 = 90.7185; %kg/s
alt = 35000/3.281; %m
M0 = 1.6;
pitotal = 82.9635;
P0_P9 = 1;

alpha = .4; %bypass ratio
beta = .01; %bleed ratio
Pto = 301.34*10^3; %Watts
h_PR = 18400*2326; %J/kg, for a (CH2)n propellant



design(2,2) = {alpha};
design(3,2) = {beta};
design(6,2) = {Pto};
design(7,2) = {pitotal};
design(8,2) = {h_PR};



year = 1995;



% Use Parametric Performance function to compute off design conditions
        % Change the way that burner fuel is calculated/ way max
        % temperature is estimated
        
        % Streamline inputs
        
% Find efficiency parameters

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