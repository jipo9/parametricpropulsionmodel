clear
clc
% close all

% Error 1: Mixer area (hardcodded in)
        % Comes from error in turbine pressure ratios
        % Maybe somewhere else?
% M6 initiial guess is correct, most of the time itll off (were gonna say
        % M6i = .5
% Inputs for off design hardcoded
% Hardcoded in mdot
        % Use alt, M0, and user input diameter?
        % Add in limiters later
% Hardcoded in F
        % Use change in enthalpy across burner to find f
        % Eventually get rid of afterburner
% Hardcode in M6A reference, find it via A16_6 estimator
% Fix compressor and fan slightly lower pressures
% Fix efficiency in nozzle function, and specific fuel consumption
% Fix  iterative plots


% Results are hardcoded into on design error readouts rn to support afterburner!!!
% Hardcoded value and fuel to air ratio change in off design afterburner
% function!!!
        

% On design is off???? fix afterburner maybe

% Need to fix fan and compressor functions
% Need to fix A estimator
% Need to fix turbine functions
% Need to fix mdot estimate
% Need to corraborate afterburner

%% To do
% Check right on design parameters and plot actual F, mdot and S


% Fix area ratio calculator!!!
% Fix burner and AB functions
% Throw into M vs T chart

%% On Design Analysis

for setup = 1:1
    stateR = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
    stateR(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
    componentR = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
    componentR(2:17,1) = {'Ram Recovery';'Inlet(Ideal,Actual)';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
    designR = {'Parameter','Value'};
    designR(2:8,1) = {'alpha';'beta';'epsilon1';'epsilon2';'PtoL';'PtoH';'h_PR'};
    inputsR = {'Parameter','Value'};
    inputsR(2:8,1) = {'Altitude (m)','Mach Number','F/mdot','Mass Flow Rate (kg/s)','SFC (kg/s/N)','Max Burner Temp (K)','Po9/P9'};
end
for input = 1:1
    % gamma_c = 1.4;
    % gamma_t = 1.3;
    % gamma_AB = 1.3;
    % cp_c = .24* 4186.8; %J/kg K
    % cp_t = .295* 4186.8; %J/kg K
    % cp_AB = .295* 4186.8; %J/kg K
    
    alt = 36000/3.281; %altitude [m from feet]
    M0 = 1.451; %freestream mach number
    F_mdot = 62.493*9.806655; %thrust/mdot [N/kg/s from lbf/(lbm/s)]
    mdot = 200*0.45359237; %freestream mass flow rate [kg/s from lbm/s]
    S = 1.0862*((.453592/3600)/4.44822); %specific fuel consuption[kg/s/N from lbm/(lbf/s)]
    T_t4 = 3200*.555556; %max burner temperature [R to K]
    T_t7 = 3600*.555556; %max burner temperature [R to K]
    Po9_P9 = 11.327; %stagnation to static pressure ratio at the nozzle (note: in future analysis calculated by total pressure ratio)
    
    alpha = .449; %bypass ratio (-1?)
    beta = .01; %bleed ratio
    PtoH = 300.61*10^3; %power takeoff high spool [W]
    PtoL = 0.000000001; %power takeoff low spool [W]
    h_PR = 18400*2326; %fuel heating value for a (CH2)n propellant [J/kg]
    
    pi_dmax = .96; %diffuser pressure ratio
    pif = 3.9; %fan pressure ratio
    ef = .89; %fan polytropic efficiency
    picL = 3.9; %low pressure compressor pressure ratio
    ecL = .89; %low pressure polytropic efficiency
    picH = 5.1282; %high pressure compressor pressure ratio
    ecH = .9; %high pressure polytropic efficiency
    eta_b = .999; %burner efficiency
    pi_b = .95; % burner pressure ratio
    etH = .89; %high pressure turbine polytropic efficiency
    etL = .9; %low pressure turbine polytropic efficiency
    etamH = .995; %high pressure shaft mechanical efficiency
    etamPH = 1; %high pressure shaft power takeoff mechancal efficiency
    etamL = .995; %low pressure shaft mechanical efficiency
    etamPL = 1; %low pressure shaft power takeoff mechancal efficiency
    pi_M_max = .97; %mixer pressure ratio
    pin = .97; %nozzle pressure ratio
    pi_AB = .95;
    % %Control limits
    % pi_c = 20;
end
for storing_values = 1:1 
    % Store all values
    designR(2,2) = {alpha}; %store values in design
    designR(3,2) = {beta};
    designR(7,2) = {PtoH};
    designR(6,2) = {PtoL};
    designR(8,2) = {h_PR};
    stateR(9,3) = {T_t4};
    stateR(16,3) = {T_t7};
    componentR(3,2) = {pi_dmax}; %store values in component
    componentR(4,2:3) = {pif,ef};
    componentR(5,2:3) = {picL,ecL};
    componentR(6,2:3) = {picH,ecH};
    componentR(7,2) = {pi_b};
    componentR(7,5) = {eta_b};
    componentR(9,5) = {etamH};
    componentR(10,5) = {etamPH};
    componentR(9,3) = {etH};
    componentR(12,5) = {etamL};
    componentR(13,5) = {etamPL};
    componentR(12,3) = {etL};
    componentR(14,2) = {pi_M_max};
    componentR(15,2) = {pi_AB};
    componentR(16,2) = {pin};
    componentR(2:3,3) = {1};
    componentR(7:8,3) = {1};
    componentR(11,3) = {1};
    componentR(14:16,3) = {1};
    inputsR(2,2) = {alt};
    inputsR(3,2) = {M0};
    inputsR(4,2) = {F_mdot};
    inputsR(5,2) = {mdot};
    inputsR(6,2) = {S};
    inputsR(7,2) = {T_t4};
    inputsR(8,2) = {Po9_P9};
end
for run_offdesign = 1:1
    [stateR,designR] = on_derived_parameters(stateR,inputsR,designR,componentR);
%     prompt = 'Run analysis with combined compressors/turbines? Y/N : ';
%     str = input(prompt,'s');
%     if str == ('Y') || str == ('y')
%         clc
%         [stateR,componentR,performanceR] = component_combined(stateR,componentR,designR,inputsR);
%     elseif str == ('N') || str == ('n')
%         clc
        [stateR,componentR,performanceR] = component_seperate(stateR,componentR,designR,inputsR);
%     else
%         error('Invalid Input. Try again and only type Y or N you dummy')
%     end
    
    err_T_mdot = performanceR{2,1} /1.0869e+03; %T/mdot error compared to book
    err_s = performanceR{2,2} / 4.7985e-05; %S error compared to book
    err_efftherm = performanceR{2,3} / .4525; %thermal efficiency error compared to book
    err_effprop =performanceR{2,4} / .4625; %propulsive efficiency compared to book
    
    fprintf('%s%.3f%s\n','Thrust                    of this analysis is ',abs(100*(1-err_T_mdot)),'% off book solution.')
    fprintf('%s%.3f%s\n','Specific Fuel Consumption of this analysis is ',abs(100*(1-err_s)),'% off book solution.')
    fprintf('%s%.3f%s\n','Thermal Efficiency        of this analysis is ',abs(100*(1-err_efftherm)),'% off book solution.')
    fprintf('%s%.3f%s\n','Propulsive Efficiency     of this analysis is ',abs(100*(1-err_effprop)),'% off book solution.')
end
for estimating_mixer_props = 1:1
    M6_est = .4;
    [A16_6,P6A_ref] = A_mixer(stateR,componentR,designR,inputsR,M6_est) %all error comes from errors in turbines
    A16_6 = .2715
    alt = inputsR{2,2};
    [~, ~, P0, ~] = atmosisa(alt);
    Po6A_P0 = 3.421 * .935 * 20 * .95 * .4231 * .4831 *.9637;
    P6A_ref = Po6A_P0 * P0
    [A45_6] = A_turb(componentR, M6_est)
end

performanceR
alpha
M6_est


% M6A_ref = .4188;
% M6A_ref = .3221;
% M6A_ref = .3746;
% [state,design,~,M6,~] = mixer(stateR,componentR,designR,A16_6,M6_est,M6A_ref,1);
% 
% alpha = design{2,2}
% M6


%% Off Design Analysis

for setup = 1:1
    state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
    state(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
    component = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
    component(2:17,1) = {'Ram Recovery';'Inlet(Ideal,Actual)';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
    design = {'Parameter','Value'};
    design(2:8,1) = {'alpha';'beta';'epsilon1';'epsilon2';'PtoL';'PtoH';'h_PR'};
    inputs = {'Parameter','Value'};
    inputs(2:8,1) = {'Altitude (m)','Mach Number','F/mdot','Mass Flow Rate (kg/s)','SFC (kg/s/N)','Max Burner Temp (K)','Po9/P9'};
end
for input = 1:1
    component = componentR;
    
        %should be automatically input, not hardcoded in
        eta_f = .8674;
        eta_cL = .8674;
        eta_cH = .8751;
        eta_tH = .8995;
        eta_tL = .9074;
%         component(:,5) = componentR(:,5); %Actually how this will look, code this in better?
%         eta_b = .999;
%         eta_PL = 1;
%         eta_PH = 1;
%         eta_AB = .99;
        component(4,5) = {eta_f};
        component(5,5) = {eta_cL};
        component(6,5) = {eta_cH};
        component(9,6) = {eta_tH};
        component(12,6) = {eta_tL};

%         pi_b = .95;
%         pi_AB = .95;
%         component(3,2) = {componentR{3,2}(1)};
%         component(7,2) = {pi_b};
%         component(15,2) = {pi_AB};
%         component(16,2) = componentR(16,2);
    
    alt = 36000/3.281; %altitude [m from feet]
    M0 = 1.451; %freestream mach number
    
%     alt = 40000/3.281; %altitude [m from feet]
%     M0 = 1.8; %freestream mach number
    inputs(2,2) = {alt};
    inputs(3,2) = {M0};
    
    design = designR;
    PtoH = 281.9*10^3; %power takeoff high spool [W]
    PtoL = 0; %power takeoff high spool [W]
%     design(2,2) = {[]};
    design(6,2) = {PtoL};
    design(7,2) = {PtoH};
    
state = stateR;
%mdot = 188.72*0.45359237;
%state{2,5} = mdot;
A0 = 5.836/10.764; %area of inlet [sqft => sqm]



    f = .02975;
    state{9,4} = f;
    fAB = .03371;
    
    M6 = .4; %inital guess
    
    
    
    
    
    
    
    
    
    
    
    M6A_ref = .4188;
end


[performance_results,inputs_results,state_results,design_results,component_results,M6] = offdesign(inputs,state,design,component,componentR,A16_6,A45_6,M6,P6A_ref,fAB,A0);
performance_results
alpha = design_results{2,2}
M6

%% Off Design Plots

n = 10;
M_SL = linspace(.1,1,n);
M_10 = linspace(.2,1.2,n);
M_20 = linspace(.3,1.5,n);
M_30 = linspace(.4,1.7,n);
M_36 = linspace(.5,1.8,n);
M_40 = linspace(.6,1.9,n);
M_50 = linspace(.7,1.9,n);

M = [M_SL;M_10;M_20;M_30;M_36;M_40;M_50];
altitude = [0,10000,20000,30000,36000,40000,50000] ./ 3.281;
fAB = 0;





% Sea level, 20 kft, 40 kft
M = [.1,.2,.3,.4,.5,.6,.7,.8,.9,1;
      .4,.5,.6,.7,.8,.9,1,1.1,1.2,1.3;
      .6,.7,.8,.9,1,1.1,1.2,1.3,1.4,1.5];
altitude = [0,20000,40000] ./ 3.281;
mdot = [142.5,145,150,155,162.5,170,182.5,190,197.5,205;
        77.5,80,85,92.5,100,110,120,130,137.5,147.5;
        37.5,40,45,50,55,60, 65,72.5,82.5,87.5] .* 0.45359237;
fAB = 0;

% f1 = figure;
% f2 = figure;
% hold on
F = zeros(size(M));
S = zeros(size(M));
for ii = 1:size(altitude,2)
    ii
    for jj = 1:size(M,2)
        jj
        state{2,5} = mdot(ii,jj);
        M0 = M(ii,jj);
        alt = altitude(ii);
        inputs(2,2) = {alt};
        inputs(3,2) = {M0};
        [performance_i,inputs_i,state_i,design_i,component_i] = offdesign(inputs,state,design,component,componentR,A16_6,A45_6,M6,P6A_ref,fAB,A0);
        % Find thurst based on F/ mdot * mdot
        T = performance_i{2,1} * state_i{2,5};
        mdot_ = state_i{2,5};
        F(ii,jj) = T;
        S(ii,jj) = performance_i{2,2};
        P6A_ref(ii,jj) = state_i{23,1};
    end
%     plot(M(ii,:),F(ii,:),'linewidth',1.5)
    %plot(M(ii,:),S(ii,:),'linewidth',1.5)

end

figure
    legend('SL','20k ft','40kft')
    title('Mach Number vs. Thrust at Various Altitudes')
    xlabel('Mach Number')
    ylabel('Thrust (N)')
    grid('on')
hold on
for ii = 1:size(altitude,2)
    plot(M(ii,:),F(ii,:),'linewidth',1.5)
end 

figure
    legend('SL','20k ft','40kft')
    title('Mach Number vs. Specific Fuel Consumption')
    xlabel('Mach Number')
    ylabel('S (kg/N)')
    grid('on')
hold on
for ii = 1:size(altitude,2)
    plot(M(ii,:),S(ii,:),'linewidth',1.5)
end 

figure
    legend('SL','20k ft','40kft')
    title('Mach Number vs. Mass Flow at Various Altitudes')
    xlabel('Mach Number')
    ylabel('Mass Flow (kg/s)')
    grid('on')
hold on
for ii = 1:size(altitude,2)
    plot(M(ii,:),mdot(ii,:))
end 

figure
    legend('SL','20k ft','40kft')
    title('Mach Number vs. Reference Pressure at Various Altitudes')
    xlabel('Mach Number')
    ylabel('Pa')
    grid('on')
hold on
for ii = 1:size(altitude,2)
    plot(M(ii,:),P6A_ref(ii,:),'linewidth',1.5)
end 


%% Printout Checks

for ref = 1:1
    state_check = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
    state_check(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
    component_check = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
    component_check(2:17,1) = {'Ram Recovery';'Inlet Actual';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
    design_check = {'Parameter','Value'};
    design_check(2:8,1) = {'alpha';'beta';'epsilon1';'epsilon2';'PtoL';'PtoH';'h_PR'};
    inputs_check = {'Parameter','Value'};
    inputs_check(2:8,1) = {'Altitude (m)','Mach Number','F/mdot','Mass Flow Rate (kg/s)','SFC (kg/s/N)','Max Burner Temp (K)','Po9/P9'};
    
    M0 = 1.451;
    alt = 36000/3.281; %altitude [m from feet]
    To4 = 3200*.555556;
    To7 = 3600*.555556;
    pi_n = .97;
    
    component_check(2:16,2) = {3.4211; .9354; 3.9000; 3.9000; 5.1282; .95; [];    .4231; []; [];    .4831; []; .9779; .95; pi_n};
    component_check(2:16,4) = {1.4211; [];    1.5479; 1.5479; 1.6803;  []; .9673; .8381; []; .9731; .8598; []; .8404;  []; []};
    
    mdot = 188.72*0.45359237;
    alpha = .449;
    beta = .01;
    ep1 = .05;
    ep2 = .05;
    f = .03070;
    fAB = .03353;
    
    state_check{2,5} = mdot;
    state_check{9,3} = To4;
    state_check{9,4} = f;
    state_check{16,3} = To7;
    
    inputs_check{2,2} = alt;
    inputs_check{3,2} = M0;
    
    design_check{2,2} = alpha;
    design_check{3,2} = beta;
    design_check{4,2} = ep1;
    design_check{5,2} = ep2;
    
    [state_check,design_check] = off_derived_parameters(state_check,design_check,fAB);
    [state_check,component_check] = a(state_check,component_check,alt,To4);
    
    F_mdot = 110.83*9.806655;
    S = 1.6941*2.8325e-05;
    eta_TH = .4525;
    eta_P = .4626;
    M9 = M0*1.531;
    
    eta_o = eta_TH*eta_P;
    
performance_check(1,:) = {'Thrust','Specific Fuel Consumption','Propulsive Efficiency','Thermal Efficiency','Overall Efficiency','Exhaust Mach'};
performance_check(2,:) = {F_mdot,S,eta_TH,eta_P,eta_o,M9};
end
for test = 1:1
%     state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
%     state(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
%     component = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
%     component(2:17,1) = {'Ram Recovery';'Inlet Actual';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
%     design = {'Parameter','Value'};
%     design(2:8,1) = {'alpha';'beta';'epsilon1';'epsilon2';'PtoL';'PtoH';'h_PR'};
%     inputs = {'Parameter','Value'};
%     inputs(2:8,1) = {'Altitude (m)','Mach Number','F/mdot','Mass Flow Rate (kg/s)','SFC (kg/s/N)','Max Burner Temp (K)','Po9/P9'};
%     
%     M0 = 1.8000;
%     alt = 40000/3.281; %altitude [m from feet]
%     To4 = 3200*.555556;
%     To7 = 3600*.555556;
%     pi_n = .97;
%     
%     component(2:16,2) = {5.7458; .9067; 3.0054; 3.0054; 4.7208; .95;    [];    .4231; []; [];    .5023; []; .9735; .95; pi_n};
%     component(2:16,4) = {1.6480; [];    1.4259; 1.4259; 1.6377;  []; .9673;    .8381; []; .9731; .8667; []; .8268;  []; []};
%     
%     mdot = 188.72*0.45359237;
%     alpha = .530;
%     beta = .01;
%     ep1 = .05;
%     ep2 = .05;
%     f = .02875;
%     fAB = .03371;
%     
%     state{2,5} = mdot;
%     state{9,3} = To4;
%     state_check{9,4} = f;
%     state{16,3} = To7;
%     
%     inputs{2,2} = alt;
%     inputs{3,2} = M0;
%     
%     design{2,2} = alpha;
%     design{3,2} = beta;
%     design{4,2} = ep1;
%     design{5,2} = ep2;
%     
%     [state_check,design_check] = off_derived_parameters(state_check,design_check,fAB);
%     [state_check,component_check] = a(state_check,component_check,alt,To4);
%     
%     F_mdot = 104.69*9.806655;
%     S = 1.7468*2.8325e-05;
%     eta_TH = .471;
%     eta_P = .5342;
%     M9 = 2.3377;
%     
%     eta_o = eta_TH*eta_P;
%     
% performance_check(1,:) = {'Thrust','Specific Fuel Consumption','Propulsive Efficiency','Thermal Efficiency','Overall Efficiency','Exhaust Mach'};
% performance_check(2,:) = {F_mdot,S,eta_TH,eta_P,eta_o,M9};
end

performance_check;









%% On Design Functions
function [state,component,performance] = component_seperate(state,component,design,inputs)
% Runs an engine analysis w/ a seperated LP and HP spools
[state, component,v0] = on_ambient(state,component,inputs);
[state,component] = on_inlet(state,component,inputs);
[state,component] = on_fan(state,component);
[state,component] = on_LPcomp(state,component);
[state,component] = on_HPcomp(state,component);
[state,component] = on_burner(state,component);
[state,component] = on_HPturb(state,component,design);
[state,component] = on_LPturb(state,component,design);
[state,component] = on_mixer(state,component);
[state,component] = on_afterburner(state,component);
[state,component,performance] = on_nozzle(state,component,inputs,v0,design);
fprintf('%s\n\n','This analysis was completed using SEPERATE high and low spools.')
end
function [state,component,performance] = component_combined(state,component,design,inputs)
% Runs an engine analysis w/ a combined LP and HP spools
[state, component,v0] = on_ambient(state,component,inputs);
[state,component] = on_inlet(state,component,inputs);
[state,component] = on_fan(state,component);
[state,component] = on_combinedcomp(state,component);
[state,component] = on_burner(state,component);
[state,component] = on_combinedturb(state,component,design);
[state,component] = on_mixer(state,component);
[state,component,performance] = on_nozzle(state,component,inputs,v0,design);
fprintf('%s\n\n','This analysis was completed using COMBINED high and low spools.')
end
function [state,design] = on_derived_parameters(state,inputs,design,component)
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
function [state, component,v0] = on_ambient(state,component,inputs)
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
function [state,component] = on_inlet(state,component,inputs)
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
function [state,component] = on_fan(state,component)
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

statet = state;
Po13i = Po2*pif; %mechanical efficiency
state(5,2) = {Po13i};
state(5,3) = {[]};
state(5,8) = {[]};
[state] = unFAIR3(state,5);
ho13i = state{5,8};
state = statet;
etaf = (ho13i-ho2)/(ho13-ho2);
component{4,5} = etaf;
end
function [state,component] = on_LPcomp(state,component)

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

statet = state;
Po25i = Po2*picl; %mechanical efficiency
state(6,2) = {Po25i};
state(6,3) = {[]};
state(6,8) = {[]};
[state] = unFAIR3(state,6);
ho25i = state{6,8};
state = statet;
etacL = (ho25i-ho2)/(ho25-ho2);
component{5,5} = etacL;
end
function [state,component] = on_HPcomp(state,component)
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

statet = state;
Po3i = Po25*pich; %mechanical efficiency
state(7,2) = {Po3i};
state(7,3) = {[]};
state(7,8) = {[]};
[state] = unFAIR3(state,7);
ho3i = state{7,8};
state = statet;
etach = (ho3i-ho25)/(ho3-ho25);
component{6,5} = etach;
end
function [state,component] = on_combinedcomp(state,component)
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
function [state,component] = on_burner(state,component)
state(8,2:3) = state(7,2:3);
state(8,6:12) = state(7,6:12);
[state] = unFAIR3(state,9);

Pro31 = state{8,2};
Pro4  = state{9,2};
% pi_b_total = Pro4 / Pro31;
% component(7,2) = { pi_b_total};

ho31 = state{7,8};
ho4 = state{9,8};
taub = ho4/ho31;
component{7,4} = taub;
end
function [state,component] = on_HPturb(state,component,design)
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
Po41 = state{10,2};

%Across turbine
% etH = .89;
fun = @(ho44) mdot41*(ho41-ho44)*etamH... %change in energy across HPturb
    -mdot3*(ho3-ho25)...                    %change in energy across HP compressor
    -(PtoH) / etamPH;                         %energy draw of takeoff power
ho44 = fzero(fun,ho41);

state(11,8) = {ho44};
[state] = unFAIR3(state,11);
Po44 = state{11,2};

tauth = ho44/ho41;
component{9,4} = tauth;
eth = component{9,3};
pitH = (state{11,2} / state{10,2})^(1/eth);
component{9,2} = pitH;

statet = state;
Po44i = Po41*pitH; %mechanical efficiency
state(11,2) = {Po44i};
state(11,3) = {[]};
state(11,8) = {[]};
[state] = unFAIR3(state,11);
ho44i = state{11,8};
state = statet;
etatH = (ho41-ho44)/(ho41-ho44i);
component{9,6} = etatH;

end
function [state,component] = on_LPturb(state,component,design)
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
Po45 = state{12,2};

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
Po5 = state{13,2};

tautl = ho5/ho45;
component{12,4} = tautl;
etl = component{12,3};
pitL = (state{13,2} / state{12,2})^(1/etl);
component{12,2} = pitL;

Po5i = Po45*pitL; %mechanical efficiency
state(13,2) = {Po5i};
state(13,3) = {[]};
state(13,8) = {[]};
[state] = unFAIR3(state,13);
ho5i = state{13,8};
state(13,2) = {Po5};
state(13,3) = {[]};
state(13,8) = {[]};
[state] = unFAIR3(state,13);
etatH = (ho45-ho5)/(ho45-ho5i);
component{12,6} = etatH;
end
function [state,component] = on_combinedturb(state,component,design)
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
function [state,component] = on_mixer(state,component)
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
function [state,component] = on_afterburner(state,component)
[state] = unFAIR3(state,16);

ho31 = state{7,8};
ho4 = state{9,8};
tauab = ho4/ho31;
component{7,4} = tauab;

        fAB = .03352;
        state{16,4} = state{16,4} + fAB;
        state{17,4} = state{16,4};
        state{18,4} = state{16,4};
end
function [state,component,performance] = on_nozzle(state,component,inputs,v0,design)
alpha = design{2,2};
beta = design{3,2};
h_PR = design{8,2};
PtoL = design{6,2};
PtoH = design{7,2};
Po9_P9 = inputs{8,2};

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
%% Off Design Functions
function [performance,inputs,state,design,component,M6] = offdesign(inputs,state,design,component,componentR,A16_6,A45_6,M6,P6A_ref,fAB,A0)
     [state, component,v0] = off_ambient(state,component,inputs,A0);
    [state,design] = off_derived_parameters(state,design,fAB);    
    [state,component] = off_inlet(state,component,inputs);
    check = 0;
    alpha_track = 0;
    while check == 0
        [state,design] = off_derived_parameters(state,design,fAB);
        [state,component] = off_fan(state,component,componentR,design); %flag
        [state,component] = off_comp(state,component,design, componentR); %flag
        [state,component] = off_burner(state,component,design,fAB);
        [state,component] = off_turb(state,component,M6,A45_6); %flag
        [state,design,check,M6,alpha_track,M6A_i] = mixer(state,component,design,inputs,A16_6,M6,P6A_ref,alpha_track);
%         M6
        alpha_track = alpha_track+1;
    end  
    state = off_afterburner(state,fAB);
    [state,component,performance] = off_nozzle(state,component,v0,design);
end
function [state,design] = off_derived_parameters(state,design,fAB)
%Derived parameters for performance model, w/ changing values

mdot = state{2,5};
alpha = design{2,2};
beta = design{3,2};
ep1 = design{4,2};
ep2 = design{5,2};
f = state{9,4};
f0 = 0;

% MASS FLOW AND FUEL TO AIR CALCULATIONS

mdot25 = mdot/(1+alpha); %after bypass leaves
mdot13 = mdot - mdot25; % bypass mass flow

mdot31 = mdot25*(1-beta - ep1 -ep2); %after bleed and coolant leaves
mdotbeta = mdot25*beta; %bleed air
mdotep1 = mdot25*ep1; %coolant air 1
mdotep2 = mdot25*ep2; %coolant air 2
mdotep = mdotep1+mdotep2; %combined coolant air

mdot_f = mdot31 * f; %mass flow per fuel/air ratio
mdot4 = mdot31 + mdot_f; %mass flow rate post-burner
f4 = mdot_f / mdot31; %fuel/air ratio post-burner

mdot41 = mdot4 + mdotep1; %mass flow rate after addtion of cooling air 1
f41 = f4*mdot4 / mdot41; %fuel/air ratio after addtion of cooling air 1

mdot45 = mdot41 + mdotep2; %mass flow rate after addtion of cooling air 2
f45 = f41*mdot41 / mdot45; %fuel/air ratio after addtion of cooling air 2

mdot6A = mdot45 + mdot13; %mass flow rate after addtion of bypass air
f6A = f45*mdot45 / mdot6A; %fuel/air ratio after addtion of bypass air

f7 = f6A + fAB;
mdot7 = mdot6A * (1+f7)/(1+f6A);

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
state(15,4) = {f6A};
state(15,5) = {mdot6A};
state(16:18,4) = {f7};
state(16:18,5) = {mdot7};

state(19,4) = {0};
state(19,5) = {mdotbeta};
state(20,4) = {0};
state(20,5) = {mdotep};
state(21,4) = {0};
state(21,5) = {mdotep1};
state(22,4) = {0};
state(22,5) = {mdotep2};


end
function [state, component,v0] = off_ambient(state,component,inputs,A0)
alt = inputs{2,2};
M0 = inputs{3,2};
[T0, ~, Pzero, rho0] = atmosisa(alt); %obtain standard atmospheric conditions
state(2,3) = {T0};
[state] = unFAIR3(state,2);
[~,~,T0,~,~,cp0,gamma0,~,~] = state{2,:};
R0 = cp0 - cp0/gamma0;
a0 = sqrt(R0*gamma0*T0); %[m/s]
v0 = M0*a0; %[m/s]


To0 = T0*(1+((M0^2)*((gamma0-1)/2))); %find total temperature using isentropic
state(3,3) = {To0};
[state] = unFAIR3(state,3);


P0 = state{2,2};
Po0 = state{3,2};
pi_r = Po0/P0;
component{2,2} = pi_r;

h0 = state{2,8};
ho0 = state{3,8};
tau_r = ho0/h0;
component{2,4} = tau_r;

% Po = Pzero*((1+((gamma0-1)/2)*M0^2)^(gamma0/(gamma0-1)));
% %mdot0 = rho0*A0*v0;
% mdot0 = (A0*Po)*sqrt(gamma0/(To0*R0))*M0*(1+((gamma0-1)/2)*M0^2)^(-1*(gamma0+1)/(2*(gamma0-1))); %compressible eqn
% state{2,5} = mdot0;
end
function [state,component] = off_inlet(state,component,inputs)
M0 = inputs{3,2};
pi_dmax = component{3,2}(1);
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
function [state,component] = off_fan(state,component,componentR,design)

ho4 = state{9,8};
h0 = state{2,8};
tau_tL = component{12,4};
tau_tH = component{9,4};
tau_cL = component{5,4};
tau_cLR = componentR{5,4};
tau_cH = component{6,4};
eta_mL = .995;
alpha = design{2,2};
beta = design{3,2};
ep1 = design{4,2};
ep2 = design{5,2};
f = state{9,4};
ho0 = state{3,8};
tau_r = ho0/h0;
PtoL = design{6,2};
eta_mPL = 1;
mdot0 = state{2,5};
tau_fR = componentR{4,4};
eta_f = componentR{4,5};

tau_lambda = ho4/h0;
tau_f = 1 ...
    + ( (1-tau_tL)*eta_mL* (((1-beta-ep1-ep2)*(1+f)*tau_lambda*tau_tH/tau_r  +  (ep1*tau_tH + ep2)*tau_cL*tau_cH))   -   (1+alpha)*PtoL/(tau_r*eta_mPL*mdot0*h0)) ...
    / ((tau_cLR - 1)/(tau_fR - 1) + alpha);


statei = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
statei(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
ho2 = state{4,8};
ho13i = ho2*(1+eta_f*(tau_f-1));
statei{5,8} = ho13i;
statei{5,4} = 0;
[statei] = unFAIR3(statei,5);
Pro13i = statei{5,2};
Pro2 = state{4,2};

pi_f = Pro13i/Pro2;

component{4,2} = pi_f;
component{4,4} = tau_f;

state{5,2} = [];
state{5,3} = [];
ho13 = ho2*tau_f;
state{5,8} = ho13;
[state] = unFAIR3(state,5);
end
function [state,component] = off_comp(state,component,design, componentR)
statei = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
statei(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};

% LP Compressor
tau_f = component{4,4};
tau_fR = componentR{4,4};
tau_cLR = componentR{5,4};
eta_cL = component{5,5};

tau_cL = 1 + (tau_f -1)*(tau_cLR - 1)/(tau_fR - 1);

ho2 = state{4,8};
ho25i = ho2*(1+eta_cL*(tau_cL-1));
statei{6,8} = ho25i;
statei{6,4} = 0;
[statei] = unFAIR3(statei,6);
Pro25i = statei{6,2};
Pro2 = state{4,2};
pi_cL = Pro25i/Pro2;

component{5,2} = pi_cL;
component{5,4} = tau_cL;

state{6,2} = [];
state{6,3} = [];
ho25 = ho2*tau_cL;
state{6,8} = ho25;
[state] = unFAIR3(state,6);



% HP Compressor
tau_r = component{2,4};
tau_tH = component{9,4};
tau_cL = component{5,4};
tau_cH = component{6,4};

eta_cH = component{6,5};
etamH = component{9,5};
etamPH = component{10,5};

alpha = design{2,2};
beta = design{3,2};
ep1 = design{4,2};
ep2 = design{5,2};
PtoH = design{7,2};

f = state{9,4};
h0 = state{2,8};
ho4 = state{9,8};
mdot0 = state{2,5};

tau_lambda = ho4/h0;

tau_cH = 1 ...
    + (1-tau_tH)*etamH* (((1-beta-ep1-ep2)*(1+f)*tau_lambda/(tau_r*tau_cL) + ep1*tau_cH))...
    - (1+alpha)*PtoH/(tau_r*tau_cL*etamPH*mdot0*h0);

ho25 = state{6,8};
ho3i = ho25*(1+eta_cH*(tau_cH-1));
statei{7,8} = ho3i;
statei{7,4} = 0;
[statei] = unFAIR3(statei,7);
Pro3i = statei{7,2};
Pro25 = state{6,2};

pi_cH = Pro3i/Pro25;

component{6,2} = pi_cH;
component{6,4} = tau_cH;

state{7,2} = [];
state{7,3} = [];
ho3 = ho25*tau_cL;
state{7,8} = ho3;
[state] = unFAIR3(state,7);
end
function [state,component] = off_burner(state,component,design,fAB)
state(8,2:3) = state(7,2:3);
state(8,6:12) = state(7,6:12);
mdot31 = state{8,3};
h_PR = design{8,2};

f_i = state{9,4};
error = 1;
while error > .01
    state{9,2} = [];
    state{9,8} = [];
    [state] = unFAIR3(state,9);

    ho31 = state{7,8};
    ho4 = state{9,8};
    taub = ho4/ho31;
    component{7,4} = taub;

    mdotf = mdot31*(ho4 - ho31) / h_PR;
    f = mdotf/mdot31;
    state{9,4} = f;
    [state,design] = off_derived_parameters(state,design,fAB);
    error = (f - f_i)/f_i;
    f_i = f;
end

end
function [state,component] = off_turb(state,component,M6,A45_6)
tau_m1 = component{8,4};
tau_tH = component{9,4};
tau_m2 = component{11,4};
ho4 = state{9,8};

gamma = 1.3;
eta_tL = component{12,6};
X = A45_6 * 1/M6 * (2/(gamma+1)    *  (1+(gamma-1)/2*M6^2)  )^  ((gamma+1)/(2*(gamma-1)));

error = 1;
up = 10;
low = 0;
while error > .0000001
    pi_tL = (up+low)/2;
    tau_tL = 1 - eta_tL*(1-pi_tL^((gamma-1)/gamma));
    X_i = pi_tL/sqrt(tau_tL);
    error = norm(X_i - X)/X;
    if X_i > X
       up =  pi_tL;
    else
       low =  pi_tL;
    end
end



%         tau_tL = 0.8598;
%         pi_tL = 0.4831;
%         
%         tau_m1 = 0.9673;
%         tau_tH = 0.8381;
%         tau_m2 = 0.9731;
%         
%         component{8,4} = tau_m1;
%         component{9,4} = tau_tH;
%         component{11,4} = tau_m2;
        

state{10,2} = [];
state{10,3} = [];
ho41 = ho4*tau_m1;
state{10,8} = ho41;
[state] = unFAIR3(state,10);

state{11,2} = [];
state{11,3} = [];
ho44 = ho41*tau_tH;
state{11,8} = ho44;
[state] = unFAIR3(state,11);

state{12,2} = [];
state{12,3} = [];
ho45 = ho44*tau_m2;
state{12,8} = ho45;
[state] = unFAIR3(state,12);

state{13,2} = [];
state{13,3} = [];
ho5 = ho45*tau_tL;
state{13,8} = ho5;
[state] = unFAIR3(state,13);

%         component{9,2} = (state{11,2} / state{10,2})^(1/.89);
%         component{9,2} = 0.4231;


component{12,2} = pi_tL;
component{12,4} = tau_tL;
end
function state = off_afterburner(state,fAB)
    if fAB == 0
        state(16,2:end) = state(15,2:end);
    else
        %Input Tt7 and f_AB in initial conditions
        [state] = unFAIR3(state,16);
        ho6A = state{15,8};
        ho7 = state{16,8};
        tau_AB = ho7/ho6A;
        component{15,4} = tau_AB;
        %This is temporary and solely for reference purposes
    end
end
function [state,component,performance] = off_nozzle(state,component,v0,design)
%Take in inputs    
alpha = design{2,2};
beta = design{3,2};
P0_P9 = 1;
h_PR = design{8,2};
PtoL = design{6,2};
PtoH = design{7,2};
[~,pi_r,pi_d,~,pi_cL,pi_cH,pi_b,~,pi_tH,~,~,pi_tL,~,pi_M,pi_AB,pi_n,~] = component{:,2};

%Calculate pressure drop across nozzle
Pro7 = state{16,2};
pin = component{16,2};
Pro9 = Pro7*pin;
state(17,2) = {Pro9};
[state] = unFAIR3(state,17);

%Calculate static conditions of exhaust
Po9_P0 = pi_r*pi_d(2)*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_M*pi_AB*pi_n;
Po9_P9 = Po9_P0/P0_P9;
Pr9 = Pro9 / Po9_P9;
state(18,2) = {Pr9};
[state] = unFAIR3(state,18);
component{16,4} = state{17,8} / state{16,8};

%Find exhause mach number
[~,~,~,~,~,~,~,ho9] = state{17,:};
[~,~,T9,~,~,cp9,gamma9,h9] = state{18,:};
R9 = cp9 - cp9/gamma9;
a9 = sqrt(R9*gamma9*T9); %m/s
v9 = sqrt(2*(ho9-h9));
M9 = v9 / a9;

%Find misc values
[~,Pr0,~,~,~,cp0,gamma0,~] = state{2,:};
R0 = cp0 - cp0/gamma0;
mdot0 = state{2,5};
f_0 = state{18,4};

%Calculate Performance
F_mdot = (1+f_0-(beta/(1+alpha)))*v9     -   v0  +   (1+f_0-(beta/(1+alpha)))*R9*T9*(1-Pr0/Pr9)/(R0*v9*gamma0);
S = f_0 / F_mdot;

eta_TH = ((v0^2/2*((1+f_0-(beta/(1+alpha)))*(v9/v0)^2 - 1) + (PtoL + PtoH)/mdot0))/...
    (f_0*h_PR);
eta_P = (2*F_mdot/v0)/...
    ((1+f_0-beta/(1+alpha))*(v9/v0)^2 - 1);
eta_o = eta_TH*eta_P;
performance(1,:) = {'Thrust','Specific Fuel Consumption','Propulsive Efficiency','Thermal Efficiency','Overall Efficiency','Exhaust Mach'};
performance(2,:) = {F_mdot,S,eta_TH,eta_P,eta_o,M9};
end
%% Misc Functions
function [state,component] = a(state,component,alt,To4)

%% State 0
[T0, ~, ~, ~] = atmosisa(alt);
state(2,3) = {T0};
[state] = unFAIR3(state,2);

%% State 0o
tau_r = component{2,4};
T_o0 = T0*tau_r;
state(3,3) = {T_o0};
[state] = unFAIR3(state,3);
P_o0 = state{3,2};

%% State 2o
pi_d = component{3,2};
P_o2 = P_o0*pi_d;
state(4,2) = {P_o2};
[state] = unFAIR3(state,4);
T_o2 = state{4,3};

%% State 1.3o
tau_f = component{4,4};
T_o13 = T_o2*tau_f;
state(5,3) = {T_o13};
[state] = unFAIR3(state,5);

%% State 2.5o
tau_cL = component{5,4};
T_o25 = T_o2*tau_cL;
state{6,3} = T_o25;
[state] = unFAIR3(state,6);

%% State 3o
tau_cH = component{6,4};
T_o3 = T_o25*tau_cH;
state{7,3} = T_o3;
[state] = unFAIR3(state,7);

%% State 3.1o
T_o31 = T_o3;
state{8,3} = T_o31;
[state] = unFAIR3(state,8);

%% State 4o
[state] = unFAIR3(state,9);

%% State 4.1o
tau_m1 = component{8,4};
T_o41 = To4*tau_m1;
state{10,3} = T_o41;
[state] = unFAIR3(state,10);

%% State 4.4o
tau_tH = component{9,4};
T_o44 = T_o41*tau_tH;
state{11,3} = T_o44;
[state] = unFAIR3(state,11);

%% State 4.5o
tau_m2 = component{11,4};
T_o45 = T_o44*tau_m2;
state{12,3} = T_o45;
[state] = unFAIR3(state,12);

%% State 5o
tau_tL = component{12,4};
T_o5 = T_o45*tau_tL;
state{13,3} = T_o5;
[state] = unFAIR3(state,13);

%% State 6o
state(14,2:end) = state(13,2:end);

%% State 6Ao
pi_mixer = component{14,2};
Pro6 = state{14,2};
Pro6A = Pro6*pi_mixer;
state{15,2} = Pro6A;
[state] = unFAIR3(state,15);

%% State 7o
[state] = unFAIR3(state,16);
end
function [A45_6] = A_turb(component, M6)
    gamma = 1.3;
    [pi_tL,~,tau_tL] = component{12,2:4};
    A45_6 = pi_tL/sqrt(tau_tL) * M6 * 1/(2/(gamma+1)*(1+(gamma-1)/2*M6^2))^((gamma+1)/(2*(gamma-1)));
end
