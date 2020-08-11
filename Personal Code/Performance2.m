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

T = 1; % N Max or Takeoff
A = 1; % Area of inlet 
pitotal = 1; %Assume constant?
alpha = 1; %bypass ratio
beta = 1; %bleed ratio
Pto = 0; %Watts
h_PR = 18400*2326; %J/kg, for a (CH2)n propellant
year = 1995;

% Add spot to hardcode in any given information!!!!
on_design(T,A,pitotal,alpha,beta,Pto,h_PR,year)

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