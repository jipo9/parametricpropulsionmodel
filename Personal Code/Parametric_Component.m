%% Read Me
% This code is a parametric analysis of a turbofan engine w/out an afterburner.

% Error outputs of 4 major performance parameters output in the command
% window, next to a P-V diagram and a T-S diagram

% All relavent information is stored in 1 of 4 cells:
% Component parameters, both given and found, are stored in
% "component". These parameters include pressure ratio, polytropic
% efficiency, enthalpy ratio, and mechanical efficiency

% Inputs general beyond components is stored in "design"

%Thermodynamic states of gas in the engine are modeled at each state
% through the engine and live in "state"

% Finally, overall performance parameters are given in "performance"
% for comparisson

% Sensitivity study is currently unfinished and returning incorrect values
%% Assumptions
%----------------------Book Assumptions----------------------%

% Assumptions used by both the book analysis and (for most) our analysis

% FLOW IS STEADY ON AVERAGE

% FLOW IS 1D THROUGH EACH COMPONENT

% FLOW BEHAVES LIKE PERFECT GAS (OUTSIDE OF BURNER?)

% fair UTILIZES NASA GLENN THERMOCHEMICAL DATA TO OBTAIN PROPERTIES OF AIR
% AND (CH2)n propellant COMBUSTION

% TOTAL PRESSURE RATIO OF DIFFUSER OR INLET IS:
% PI_TOTAL = PI_MAX*ETA
% WHERE PI_MAX IS PRESSURE RATIO CAUSED BY WALL FRICTION AND ETA IS RAM
% RECOVERY OF MIL-E-5008B

% FAN, LP COMPRESSOR, AND POWER TAKEOFF(LP) POWERED BY LP TURBINE (FURTHER
% ASSUMPTION IN COMBINED STATIONS)

% HP COMPRESSOR AND POWER TAKEOFF(HP) POWERED BY HP TURBINE (FURTHER
% ASSUMPTION IN COMBINED STATIONS)

% HIGH PRESSURE BLEED AIR AND TURBINE COOLING AIR REMOVED BETWEEN HP
% COMPRESSOR AND BURNER

% INNEFICIENCIES IN TURBINE AND COMPRESSOR PRESSURE RATIOS ACOUNTED W/
% POLYTROPIC INNEFICIENCIES, AND ARE LOWER AT TURBINES DUE TO COOLING AIR

% FAN AND CORE STREAMS MIX COMPLETELY IN THE MIXER W/ A PRESSURE RATIO:
% PI_TOTAL = PI_IDEAL*PI_MAX
% WHERE PI_IDEAL IS TOTAL PRESSURE RATIO ACROSS IDEAL CONSTANT AREA MIXER
% AND PI_MAX IS PRESSURE RATIO CAUSED BY WALL FRICTION

% COOLANT AIR DEFINED BY BURNER EXIT TEMP (SEE

%-------------------Additional Assumptions-------------------%

% ALL ANALYSIS DONE IN METRIC

% THRUST AND SPECIFIC FUEL CONSUMPTION ASSUMED AND USED TO FIND THE FUEL TO
% AIR RATIO

% NO AFTERBURNER AND NO PRESSURE RATIO ACROSS

% SIMPLIFIED MIXER (and stagnation)
% Mixer is assumed to be independent of mass flow through the bypass.
% (PI_IDEAL = 1)
% Final Pr -.5%
% T -.1%
% Note: This data is on a low bypass engine, will most likely matter more
% for high bypass engines Also, value was close to 1 anyways, so older
% engines will probably be much less accurate

% NOZZLE IS PERFECTLY EXPANDED (PO = P9)

% COMBINED STATION
% Using combined functions on turbine (HP turbine, Coolant 1, LP turbine,
% Coolant 2) and compressors (LP compressor, HP compressor) yields (from the book):
% Final Pr +22%
% +4.5% from compressor, +17% from turbines
% T +4.5%
% S -4.3%
% eta_thermal +6%
% eta_prop -1.8%
% With the added benefit of:
% Ptol and Ptoh simplify to just 1 variable
% pi_cH and pi_cL simplify to just 1 variable
% Assumptions
% LP and HP shafts are one in the same
% All values are averaged between shafts
% Power takeoff is combined
% All coolant air is mixed prior to the first turbine
% I tried to vary mixing ratios a bit, seems to affect
% results negatively to switch things around too much...

%--------------------Future Assumptions--------------------%

% Assumptions that will be used in comparison against actual engines

% MOST ENGINE PARAMETERS "LEVEL OF TECH" based

% TURBINE MECHANICAL EFFICIENCIES AND MIXER PRESSURE RATIO  = 1
% Assuming mechanical efficiency of both spools and takeoff power is at 1
% differs under this analysis by T = +.1%

% Mixer is assumed to be 1
% Final Pr +2.3%
% T +.5%

% Note: This data is on a low bypass engine, mixer will most likely matter more
% for high bypass engines Also, these values were all very close to 1 anyways, so older
% engines will probably be much less accurate

%% FUTURE WORK (ON COMPONENT MODEL):
% Fix combined turbine function
%   -> added fan contribution. Pressure ratios are closer but overall...
%   ...performance values are way off
% Fix sensitivity study
% no pi_b or pi_AB
% Tt/T values are off, need to find out why


%% Initialize cells
clc
clear
close all
state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
state(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
component = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
component(2:17,1) = {'Ram Recovery';'Inlet(Ideal,Actual)';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
design = {'Parameter','Value'};
design(2:8,1) = {'alpha';'beta';'epsilon1';'epsilon2';'PtoL';'PtoH';'h_PR'};
inputs = {'Parameter','Value'};
inputs(2:8,1) = {'Altitude (m)','Mach Number','F/mdot','Mass Flow Rate (kg/s)','SFC (kg/s/N)','Max Burner Temp (K)','Po9/P9'};

%% Initial Conditions
%--------Overall Performance (Givens)---------
alt = 35000/3.281; %altitude [m from feet]
M0 = 1.6; %freestream mach number
F_mdot = 62.859*9.806655; %thrust/mdot [N/kg/s from lbf/(lbf/s)]
mdot = 200*0.45359237; %freestream mass flow rate [kg/s from lbm/s]
S = 1.1386*((.453592/3600)/4.44822); %specific fuel consuption[kg/s/N from lbm/(lbf/s)]
T_t4 = 3200*.555556; %max burner temperature [R to K]
Po9_P9 = 12.745; %stagnation to static pressure ratio at the nozzle (note: in future analysis calculated by total pressure ratio)

alpha = .4; %bypass ratio
beta = .01; %bleed ratio
PtoH = 301.34*10^3; %power takeoff high spool [W]
PtoL = 0.000000001; %power takeoff low spool [W]
h_PR = 18400*2326; %fuel heating value for a (CH2)n propellant [J/kg]

pi_dmax = .96; %diffuser pressure ratio
pif = 3.8; %fan pressure ratio
ef = .89; %fan polytropic efficiency
picL = 3.8; %low pressure compressor pressure ratio
ecL = .89; %low pressure polytropic efficiency
picH = 4.2105; %high pressure compressor pressure ratio
ecH = .9; %high pressure polytropic efficiency
eta_b = .999; %burner efficiency
pi_b = .95; % burner pressure ratio
etH = .89; %high pressure turbine polytropic efficiency
etL = .9; %low pressure turbine polytropic efficiency
etamH = .995; %high pressure shaft mechanical efficiency
etamPH = .99; %high pressure shaft power takeoff mechancal efficiency
etamL = .995; %low pressure shaft mechanical efficiency
etamPL = 1; %low pressure shaft power takeoff mechancal efficiency
pi_M_max = .97; %mixer pressure ratio
pin = .97; %nozzle pressure ratio

% Store all values
design(2,2) = {alpha}; %store values in design
design(3,2) = {beta};
design(7,2) = {PtoH};
design(6,2) = {PtoL};
design(8,2) = {h_PR};
state(9,3) = {T_t4};
component(3,2) = {pi_dmax}; %store values in component
component(4,2:3) = {pif,ef};
component(5,2:3) = {picL,ecL};
component(6,2:3) = {picH,ecH};
component(7,2) = {pi_b};
component(7,5) = {eta_b};
component(9,5) = {etamH};
component(10,5) = {etamPH};
component(9,3) = {etH};
component(12,5) = {etamL};
component(13,5) = {etamPL};
component(12,3) = {etL};
component(14,2) = {pi_M_max};
component(16,2) = {pin};
component(2:3,3) = {1};
component(7:8,3) = {1};
component(11,3) = {1};
component(14:16,3) = {1};
inputs(2,2) = {alt};
inputs(3,2) = {M0};
inputs(4,2) = {F_mdot};
inputs(5,2) = {mdot};
inputs(6,2) = {S};
inputs(7,2) = {T_t4};
inputs(8,2) = {Po9_P9};

%% Derived Parameters
[state,design] = derived_parameters(state,inputs,design,component);
state_joey = state;
design_joey = design;
inputs_joey = inputs;
component_joey = component;

%% Analysis

prompt = 'Run analysis with combined compressors/turbines? Y/N : ';
str = input(prompt,'s');
if str == ('Y') || str == ('y')
    clc
    [state,component,performance] = component_combined(state,component,design,inputs);
elseif str == ('N') || str == ('n')
    clc
    [state,component,performance] = component_seperate(state,component,design,inputs);
else
    error('Invalid Input. Try again and only type Y or N you dummy')
end

%% Sensitivity Study
% state_i = state;
% component_i = component;
% alt_i = inputs{2,2};
% M0_i = inputs{3,2};
% mdotep1_i = state{21,5};
% Po9_P9_i = inputs{8,2};
% design_i = design;


% Define inputs
% Free: alt M0 Po9_P9
% Design: alpha beta PtoH h_PR
% Component: pi_dmax pif ef picl ecl pich ech pi_b etamH etamPH etH etamL etamPL etL pi_M_max pin T_t4];
% State: T_t4 (and thus mdot)
%
% design(2,2) = {alpha}; %store values in design
% design(3,2) = {beta};
% design(7,2) = {PtoH};
% design(6,2) = {PtoL};
% design(8,2) = {h_PR};
%
% component(3,2) = {pi_dmax}; %store values in component
% component(4,2:3) = {pif,ef};
% component(5,2:3) = {picL,ecL};
% component(6,2:3) = {picH,ecH};
% component(9,5) = {etamH};
% component(10,5) = {etamPH};
% component(9,3) = {etH};
% component(12,5) = {etamL};
% component(13,5) = {etamPL};
% component(12,3) = {etL};
% component(14,2) = {pi_M_max};
% component(16,2) = {pin};
%
%
% state(9,3) = {T_t4};


% A = [3,2;%pi_dmax
%     4,2;%pif
%     4,3;%ef
%     5,2;%picL
%     5,3;%ecL
%     6,2;%picH
%     6,3;%ecH
%     9,5;%etamH
%     10,5;%etamPH
%     9,3;%etH
%     12,5;%etamL
%     13,5;%etamPL
%     12,3;%etL
%     14,2;%pi_M_max
%     16,2];%pin
% %F_mdot = performance_s{2,1};
% %S = performance_s{2,2};
% 
% for n = 1:15
%     [ii] = A(n,1);
%     [jj] = A(n,2);
%     component_i(ii,jj) = {n*component_i{ii,jj}};
%     delta_val(n) = .01;
%     [error(n,:)] = sensitivity(S,F_mdot,state_i,component_i,alt_i,M0_i,mdotep1_i,Po9_P9_i,design_i,delta_val(n),inputs);
%     component_i(ii,jj) = {component_i{ii,jj}/n};
% end
% 
% alt_i = alt;
% M0_i = M0;
%     alt_i = n*alt_i;
%     delta_val(16) = .01;
%     [error(16,:)] = sensitivity(S,F_mdot,state_i,component_i,alt_i,M0_i,mdotep1_i,Po9_P9_i,design_i,delta_val(16),inputs);
%     alt_i = alt_i/n;
% 
%     M0_i = n*M0_i;
%     delta_val(17) = .01;
%     [error(17,:)] = sensitivity(S,F_mdot,state_i,component_i,alt_i,M0_i,mdotep1_i,Po9_P9_i,design_i,delta_val(17),inputs);
%     M0_i = M0_i/n;

%% Results
err_T_mdot = performance{2,1} /F_mdot; %T/mdot error compared to book
err_s = performance{2,2} / S; %S error compared to book
err_efftherm = performance{2,3} / .5589; %thermal efficiency error compared to book
err_effprop =performance{2,4} / .6162; %propulsive efficiency compared to book

fprintf('%s%.3f%s\n','Thrust                    of this analysis is ',abs(100*(1-err_T_mdot)),'% off book solution.')
fprintf('%s%.3f%s\n','Specific Fuel Consumption of this analysis is ',abs(100*(1-err_s)),'% off book solution.')
fprintf('%s%.3f%s\n','Thermal Efficiency        of this analysis is ',abs(100*(1-err_efftherm)),'% off book solution.')
fprintf('%s%.3f%s\n','Propulsive Efficiency     of this analysis is ',abs(100*(1-err_effprop)),'% off book solution.')

%% Visuals

if str == ('N') || str == ('n') %only run pv and ts diagrams if using 'seperate' function
    
    %T-s diagram
    for i = 1:15
        To(i) = state{i+2,3};
        s(i) = state{i+2,9};
        
        Ravg = (state{i+2,10}+state{3,10})/2;
        pratio = (state{i+2,2}/state{3,2});
        ds(i) = state{i+2,9} - state{3,9} - (Ravg*log(pratio));
    end
    
    h1 = figure(1);
    h1.WindowStyle = 'docked';
    %plot(ds,To,'linewidth',2)
    plot([ds(2),ds(6),ds(7),ds(15)],[To(2),To(6),To(7),To(15)],'linewidth',2)
    title('T-s Diagram for Turbofan Engine')
    xlabel('Entropy (s)')
    ylabel('Total Temperature (To)')
    grid('on')
    
    %P-v diagram
    for i = 1:14
        Vr(i) = state{i+2,12};
        Por(i) = state{i+2,2};
    end
    
    h2 = figure(2);
    h2.WindowStyle = 'docked';
    plot(Vr,Por,'linewidth',2)
    title('P-V Diagram for Turbofan Engine')
    xlabel('Relative Volume')
    ylabel('Relative Pressure')
    grid('on')
    
end

% Error Bar chart
bookthrme = .5589;
bookprope = .6162;

h3 = figure(3);
h3.WindowStyle = 'docked';
sgtitle('Error Analysis')
x = categorical({'Calculated','Book'});
x = reordercats(x,{'Calculated','Book'});

subplot(2,2,1)
bc = [performance{2,1}; F_mdot];
bar(x,bc,.4)
title('Thrust per unit mass flow')

subplot(2,2,2)
bc = [performance{2,2}; S];
bar(x,bc,.4)
title('Specific Fuel Consumption')

subplot(2,2,3)
bc = [performance{2,3}; bookprope];
bar(x,bc,.4)
title('Propulsive Efficiency')

subplot(2,2,4)
bc = [performance{2,4}; bookthrme];
bar(x,bc,.4)
title('Thermal Efficiency')

h4 = figure(4);
h4.WindowStyle = 'docked';
x = categorical({'Thrust','SFC','Thermal','Propulsive'});
x = reordercats(x,{'Thrust','SFC','Thermal','Propulsive'});
bc = [err_T_mdot; err_s; err_efftherm; err_effprop];
bar(x,bc*100,.4,'FaceColor','flat')
hold on
yline(100,'--r','linewidth',1.5);
title('Percent Accuracy')
ylabel('% Accuracy Compared to Book')

%% Sens2
sens = {'Parameter','Value','Effect on F','Effect on S'};
sens(2:30,1) = {'alt';'M0';'mdot';'F_mdot';'s';'alpha';'beta';'Ptoh';'PtoL';'h_PR';'pi_dmax';'pif';'ef';'picL';'ecL';'picH';'ecH';'eta_b';'pi_b';'etamH';'etamPH';'etH';'etamL';'etamPL';'etL';'pi_M_max';'pin';'T_t4';'Po9_P9'};
sens(2:30,2) = {alt M0 mdot F_mdot S alpha beta PtoH PtoL h_PR pi_dmax pif ef picL ecL picH ecH eta_b pi_b etamH etamPH etH etamL etamPL etL pi_M_max pin T_t4 Po9_P9};

performance_f = performance;


for i = 2:length(sens)
    clear state;clear component;clear design;clear inputs
    state_f = state_joey;
    component_f = component_joey;
    design_f = design_joey;
    inputs_f = inputs_joey;
    [inputs,design,component,state,changedval] = sensedit(state_f,component_f,design_f,inputs_f,i); %sw1 == change, sw2==revert
    [state,design] = derived_parameters(state,inputs,design,component);
    [~,~,performance_i] = component_seperate(state,component,design,inputs);
    errT = ((performance_i{2,1} - performance_f{2,1})/(performance_f{2,1}))/((changedval(1,1)-sens{i,2})/(sens{i,2}));
    errS = ((performance_i{2,2} - performance_f{2,2})/(performance_f{2,2}))/((changedval(1,1)-sens{i,2})/(sens{i,2}));
    sens{i,3} = errT; %take abs out eventually
    sens{i,4} = errS;
    
    %     sens{i,3} = abs(errT); %take abs out eventually
%     sens{i,4} = abs(errS);
end

h5 = figure(5);
h5.WindowStyle = 'docked';
pie(cell2mat(sens(2:end,3)),{'Alt','M0','','F/mdot','SFC','\alpha','\beta','','PtoL','h_PR','','\pi_f','','\pi_cL','ecL','\pi_cH','ecH','\eta_b','\pi_b','\eta_mH','','etH','\eta_mL','\eta_mPL','etL','\pi_Mmax','\pi_n','To4','Po9/P9'})
warning off
title('Input Sensitivity to Overall Thrust')

h6 = figure(6);
h6.WindowStyle = 'docked';
pie(cell2mat(sens(2:end,4)),{'Alt','M0','','F/mdot','SFC','\alpha','\beta','','PtoL','h_PR','','\pi_f','','\pi_cL','ecL','\pi_cH','ecH','\eta_b','\pi_b','\eta_mH','','etH','\eta_mL','\eta_mPL','etL','\pi_Mmax','\pi_n','To4','Po9/P9'})
warning off
title('Input Sensitivity to Specific Fuel Consumption')

%% Functions
% High level analysis functions

function [state,component,performance] = component_seperate(state,component,design,inputs)
% Runs an engine analysis w/ a seperated LP and HP spools
[state, component,v0] = ambient(state,component,inputs);
[state,component] = inlet(state,component,inputs);
[state,component] = fan(state,component);
[state,component] = LPcomp(state,component);
[state,component] = HPcomp(state,component);
[state,component] = burner(state,component);
[state,component] = HPturb(state,component,design);
[state,component] = LPturb(state,component,design);
[state,component] = mixer(state,component);
[state,component,performance] = nozzle(state,component,inputs,v0,design);
fprintf('%s\n\n','This analysis was completed using SEPERATE high and low spools.')
end

function [state,component,performance] = component_combined(state,component,design,inputs)
% Runs an engine analysis w/ a combined LP and HP spools
[state, component,v0] = ambient(state,component,inputs);
[state,component] = inlet(state,component,inputs);
[state,component] = fan(state,component);
[state,component] = combinedcomp(state,component);
[state,component] = burner(state,component);
[state,component] = combinedturb(state,component,design);
[state,component] = mixer(state,component);
[state,component,performance] = nozzle(state,component,inputs,v0,design);
fprintf('%s\n\n','This analysis was completed using COMBINED high and low spools.')
end

function [inputs,design,component,state,changedval] = sensedit(state,component,design,inputs,index)

n = .95;

    if index == 2
        inputs{2,2} = inputs{2,2}*n;
        changedval = inputs{2,2};
    elseif index == 3
        inputs{3,2} = inputs{3,2}*n;
        changedval = inputs{3,2};
    elseif index == 4
        inputs{5,2} = inputs{5,2}*n;
        changedval = inputs{5,2};
    elseif index == 5
        inputs{4,2} = inputs{4,2}*n;
        changedval = inputs{4,2};
    elseif index == 6
        inputs{6,2} = inputs{6,2}*n;
        changedval = inputs{6,2};
    elseif index == 7
        design{2,2} = design{2,2}*n;
        changedval = design{2,2};
    elseif index == 8
        design{3,2} = design{3,2}*n;
        changedval = design{3,2};
    elseif index == 9
        design{7,2} = design{7,2}*n;
        changedval = design{7,2};
    elseif index == 10
        design{6,2} = design{6,2}*n;
        changedval = design{6,2};
    elseif index == 11
        design{8,2} = design{8,2}*n;
        changedval = design{8,2};
    elseif index == 12
        component{3,2} = component{3,2}*n;
        changedval = component{3,2};
    elseif index == 13
        component{4,2} = component{4,2}*n;
        changedval = component{4,2};
    elseif index == 14
        component{4,3} = component{4,3}*n;
        changedval = component{4,3};
    elseif index == 15
        component{5,2} = component{5,2}*n;
        changedval = component{5,2};
    elseif index == 16
        component{5,3} = component{5,3}*n;
        changedval = component{5,3};
    elseif index == 17
        component{6,2} = component{6,2}*n;
        changedval = component{6,2};
    elseif index == 18
        component{6,3} = component{6,3}*n;
        changedval = component{6,3};
    elseif index == 19
        component{7,5} = component{7,5}*n;
        changedval = component{7,5};
    elseif index == 20
        component{7,2} = component{7,2}*n;
        changedval = component{7,2};
    elseif index == 21
        component{9,5} = component{9,5}*n;
        changedval = component{9,5};
    elseif index == 22
        component{10,5} = component{10,5}*n;
        changedval = component{10,5};
    elseif index == 23
        component{9,3} = component{9,3}*n;
        changedval = component{9,3};
    elseif index == 24
        component{12,5} = component{12,5}*n;
        changedval = component{12,5};
    elseif index == 25
        component{13,5} = component{13,5}*n;
        changedval = component{13,5};
    elseif index == 26
        component{12,3} = component{12,3}*n;
        changedval = component{12,3};
    elseif index == 27
        component{14,2} = component{14,2}*n;
        changedval = component{14,2};
    elseif index == 28
        component{16,2} = component{16,2}*n;
        changedval = component{16,2};
    elseif index == 29
        inputs{7,2} = inputs{7,2}*n;
        state{9,3} = state{9,3}*n;
        changedval = inputs{7,2};
    elseif index == 30
        inputs{8,2} = inputs{8,2}*n;
        changedval = inputs{8,2};
    end
   
    
end

function [error] = sensitivity(S,F_mdot,state_i,component_i,alt_i,M0_i,mdotep1_i,Po9_P9_i,design_i,delta_val,inputs);
[~,~,performance_i] = component_seperate(state_i,component_i,design_i,inputs);
F_mdot_i = performance_i{2,1};
S_i = performance_i{2,2};

delta_F = (F_mdot_i - F_mdot)/9.806;
delta_S = (S_i - S)/S;
error(1,:) = [delta_F/delta_val,delta_S/delta_val];
end


function [state,design] = derived_parameters(state,inputs,design,component)
%% Derived Parameters
% COOLING AIR CALCULATIONS
% Found in Aircraft Engine Design - Mattingly but
% unable to find rationale...

alt = inputs{2,2};
M0 = inputs{3,2};
F_mdot = inputs{4,2};
mdot = inputs{5,2};
S = inputs{6,2};
T_t4 = inputs{7,2};
alpha = design{2,2};
beta = design{3,2};
eta_b = component{7,5};


if T_t4 > 2400*.55556
    ep1 = (T_t4/.5556-2400)/(16000);
    ep2 = ep1;
else
    ep1 = 0;
    ep2 = 0;
end

design(4,2) = {ep1};
design(5,2) = {ep2};


% MASS FLOW AND FUEL TO AIR CALCULATIONS
mdot_f = S*F_mdot*mdot /eta_b; %mass flow per fuel/air ratio

f0 = 0; %freestream fuel/air ratio

mdot25 = mdot/(1+alpha); %after bypass leaves
mdot13 = mdot - mdot25; % bypass mass flow

mdot31 = mdot25*(1-beta - ep1 -ep2); %after bleed and coolant leaves
mdotbeta = mdot25*beta; %bleed air
mdotep1 = mdot25*ep1; %coolant air 1
mdotep2 = mdot25*ep2; %coolant air 2
mdotep = mdotep1+mdotep2; %combined coolant air

mdot4 = mdot31 + mdot_f; %mass flow rate post-burner
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
state(12:13,4) = {f45};
state(12:13,5) = {mdot45};
state(14:18,4) = {f6A};
state(14:18,5) = {mdot6A};
state(19,5) = {mdotbeta};
state(20,5) = {mdotep};
state(21,5) = {mdotep1};
state(22,5) = {mdotep2};

end


% Component analysis functions
function [state, component,v0] = ambient(state,component,inputs)
alt = inputs{2,2};
M0 = inputs{3,2};
[T0, ~, ~, ~] = atmosisa(alt); %obtain standard atmospheric conditions
state(2,3) = {T0};
[state] = unFAIR3(state,2);
[~,~,T0,~,~,cp0,gamma0,~,~] = state{2,:};
R0 = cp0 - cp0/gamma0;
a0 = sqrt(R0*gamma0*T0); %[m/s]
v0 = M0*a0; %[m/s]

T_o0 = T0*(1+((M0^2)*((gamma0-1)/2))); %find total temperature using isentropic
state(3,3) = {T_o0};
[state] = unFAIR3(state,3);


P0 = state{2,2};
Po0 = state{3,2};
pi_r = Po0/P0;
component{2,2} = pi_r;

h0 = state{2,8};
ho0 = state{3,8};
tau_r = ho0/h0;
component{2,4} = tau_r;
end

function [state,component] = inlet(state,component,inputs)
M0 = inputs{3,2};
pi_dmax = component{3,2};
Po0 = state{3,2};

if M0<1
    pid = pi_dmax;
elseif M0>1 && M0<5
    pid = pi_dmax * (1-.075*((M0-1)^1.35));
else
    pid = pi_dmax * (800/(M0^4 + 935));
end

Po2 = pid*Po0;

component{3,2} = [pi_dmax,pid];
state(4,2) =  {Po2};
[state] = unFAIR3(state,4);

ho0 = state{3,8};
ho2 = state{4,8};
tau_d = ho2/ho0;
component{3,4} = tau_d;
end

function [state,component] = fan(state,component)
pif = component{4,2};
ef = component{4,3};
Po2 = state{4,2};

Po13 = Po2*pif^(1/ef);
state(5,2) = {Po13};
[state] = unFAIR3(state,5);


ho2 = state{4,8};
ho13 = state{5,8};
tauf = ho13/ho2;
component{4,4} = tauf;
end

function [state,component] = LPcomp(state,component)

picl = component{5,2};
ecl = component{5,3};
Po2 = state{4,2};

Po25 = Po2*picl^(1/ecl);
state(6,2) = {Po25};
[state] = unFAIR3(state,6);

ho2 = state{4,8};
ho25 = state{6,8};
taucl = ho25/ho2;
component{5,4} = taucl;
end

function [state,component] = HPcomp(state,component)
pich = component{6,2};
ech = component{6,3};
Po25 = state{6,2};

Po3 = Po25*pich^(1/ech);
state(7,2) = {Po3};
[state] = unFAIR3(state,7);

ho25 = state{6,8};
ho3 = state{7,8};
taucl = ho3/ho25;
component{6,4} = taucl;
end

function [state,component] = combinedcomp(state,component)
picl = component{4,2};
pich = component{5,2};
pic = picl*pich;

ecl = component{4,3};
ech = component{5,3};
ec = (ecl+ech)/2;

Po2 = state{4,2};

Po3 = Po2*pic^(1/ec);
state(7,2) = {Po3};
[state] = unFAIR3(state,7);
end

function [state,component] = burner(state,component)
state(8,2:3) = state(7,2:3);
state(8,6:12) = state(7,6:12);
[state] = unFAIR3(state,9);

Pro31 = state{8,2};
Pro4  = state{9,2};
pi_b_total = Pro4 / Pro31;
component(7,2) = { pi_b_total};

ho31 = state{7,8};
ho4 = state{9,8};
taub = ho4/ho31;
component{7,4} = taub;
end

function [state,component] = HPturb(state,component,design)
mdot3 = state{7,5};
mdot4 = state{9,5};
mdot41 = state{10,5};
ho25 = state{6,8};
ho3 = state{7,8};
hoep1 = state{8,8};
ho4 = state{9,8};
etamH = component{9,5};
etamPH = component{10,5};
PtoH = design{7,2};
mdotep1 = state{21,5};

%Across mixer
ho41 = (mdotep1*hoep1 + mdot4*ho4) / (mdot41);
state(10,8) = {ho41};
[state] = unFAIR3(state,10);

taum1 = ho41/ho4;
component{8,4} = taum1;
pim1 = state{10,2} / state{9,2};
component{8,2} = pim1;

%Across turbine
% etH = .89;
fun = @(ho44) mdot41*(ho41-ho44)*etamH... %change in energy across HPturb
    -mdot3*(ho3-ho25)...                    %change in energy across HP compressor
    -(PtoH) / etamPH;                         %energy draw of takeoff power
ho44 = fzero(fun,ho41);

state(11,8) = {ho44};
[state] = unFAIR3(state,11);

tauth = ho44/ho41;
component{9,4} = tauth;
pitH = state{11,2} / state{10,2};
component{9,2} = pitH;
end

function [state,component] = LPturb(state,component,design)
mdot13 = state{5,5};
mdot25 = state{6,5};
mdot5 = state{13,5};
mdot44 = state{11,5};
mdot45 = state{12,5};
ho2 = state{4,8};
ho13 = state{5,8};
ho25 = state{6,8};
hoep2 = state{8,8};
ho44 = state{11,8};
etamL = component{12,5};
etamPL = component{13,5};
PtoL = design{6,2};
mdotep2 = state{22,5};

%Across mixer
ho45 = (mdotep2*hoep2 + mdot44*ho44) /(mdot45);
state(12,8) = {ho45};
[state] = unFAIR3(state,12);

taum2 = ho45/ho44;
component{11,4} = taum2;
pim2 = state{12,2} / state{11,2};
component{11,2} = pim2;

%Across Turbine
fun = @(ho5) mdot5*(ho45 - ho5)*etamL...    %change in energy across LP turb
    -mdot25*(ho25-ho2)...                   %change in energy across LP compressor
    -mdot13*(ho13-ho2)...                    %change in energy across fan
    -PtoL / etamPL;                         %energy draw of takeoff power
ho5 = fzero(fun,ho45);

state(13,8) = {ho5};
[state] = unFAIR3(state,13);

tautl = ho5/ho45;
component{12,4} = tautl;
pitL = state{13,2} / state{12,2};
component{12,2} = pitL;
end

function [state,component] = combinedturb(state,component,design)
hoep = state{8,8};
ho4 = state{9,8};
ho25 = state{6,8};
ho13 = state{5,8};
ho2 = state{4,8};
mdot13 = state{5,5};
mdot25 = state{6,5};
mdot4 = state{9,5};
mdot45 = state{12,5};
mdotep = state{20,5};

ho45 = (mdotep*hoep + mdot4*ho4) /(mdot45);
state(12,8) = {ho45};
[state] = unFAIR3(state,12);


etamL = component{12,5};
etamPL = component{13,5};
PtoL = design{6,2};
etamH = component{9,5};
etamPH = component{10,5};
PtoH = design{7,2};


etam = (etamH+etamL)/2;
etamP = (etamPH+etamPL)/2;

mdot3 = state{7,5};
ho3 = state{7,8};
ho2 = state{4,8};
Pto = PtoH + PtoL;

fun = @(ho5) mdot45*(ho4-ho5)*etam... %change in energy across turbine
    -mdot3*(ho3-ho13)...              %change in energy across comp
    -mdot13*(ho13-ho2)...             %change in energy across fan
    -(PtoL) / etamPL...               %energy draw off LP spool
    -(PtoH) / etamPH;                 %energy draw off HP spool

ho5 = fzero(fun,ho45);
state(13,8) = {ho5};
[state] = unFAIR3(state,13);

taut = ho5/ho45;
component{12,4} = taut;
component{9,4} = taut;
pit = state{13,2} / state{12,2};
component{12,2} = pit;
component{9,2} = pit;
end

function [state,component] = mixer(state,component)
%Assume ideal pressure ratio (mach independent)
%Assume perfect polytropic efficiency

piM = component{14,2};
hoalpha = state{5,8};
ho5 = state{13,8};
mdotalpha = state{5,5};
mdot5 = state{13,5};
mdot6 = state{14,5};

ho6 = (mdotalpha*hoalpha + mdot5*ho5) / (mdot6);
state(14,8) = {ho6};
[state] = unFAIR3(state,14);

Po6 = state{14,2};
Po6A = Po6*piM;
state(15,2) = {Po6A};
[state] = unFAIR3(state,15);

tauM = ho6/ho5;
component{14,4} = tauM;
end

function [state,component] = afterburner(state,component)
[state] = unFAIR3(state,16);

Pro31 = state{8,2};
Pro4  = state{9,2};
pi_ab_total = Pro4 / Pro31;
component(7,2) = { pi_ab_total};

ho31 = state{7,8};
ho4 = state{9,8};
tauab = ho4/ho31;
component{7,4} = tauab;
end

function [state,component,performance] = nozzle(state,component,inputs,v0,design)
alpha = design{2,2};
beta = design{3,2};
h_PR = design{8,2};
PtoL = design{6,2};
PtoH = design{7,2};
Po9_P9 = inputs{8,2};

state(16,2:end) = state(15,2:end); %assume no afterburner

%Calculate pressure drop across nozzle
Pro7 = state{16,2};
pin = component{16,2};
Pro9 = Pro7*pin;
state(17,2) = {Pro9};
[state] = unFAIR3(state,17);

%Calculate static conditions of exhaust
Pr9 = Pro9 / Po9_P9;
state(18,2) = {Pr9};
component{16,4} = state{17,8} / state{16,8};

[state] = unFAIR3(state,18);

[~,~,~,~,~,~,~,ho9] = state{17,:};
[~,~,T9,~,~,cp9,gamma9,h9] = state{18,:};
R9 = cp9 - cp9/gamma9;
a9 = sqrt(R9*gamma9*T9); %m/s
v9 = sqrt(2*(ho9-h9));
M9 = v9 / a9;

[~,Pr0,~,~,~,cp0,gamma0,~] = state{2,:};
R0 = cp0 - cp0/gamma0;
mdot0 = state{2,5};
f_0 = state{18,4};




F_mdot = (1+f_0-(beta/(1+alpha)))*v9     -   v0  +   (1+f_0-(beta/(1+alpha)))*R9*T9*(1-Pr0/Pr9)/(R0*v9*gamma0);
S = f_0 / F_mdot;
eta_TH = ((v0^2/2*((1+f_0-(beta/(1+alpha)))*(v9/v0)^2 - 1) + (PtoL + PtoH)/mdot0))/...
    (f_0*h_PR);
%eta_P = 2/(1+v9/v0); Simplified case, neglected for now
eta_P = (2*F_mdot/v0)/...
    ((1+f_0-beta/(1+alpha))*(v9/v0)^2 - 1);
eta_o = eta_TH*eta_P;
% eta_o = v0/(h_PR*S); Alternate eqn
performance(1,:) = {'Thrust','Specific Fuel Consumption','Propulsive Efficiency','Thermal Efficiency','Overall Efficiency','Exhaust Mach'};
performance(2,:) = {F_mdot,S,eta_TH,eta_P,eta_o,M9};
end

