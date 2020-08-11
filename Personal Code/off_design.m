function [state,component,design, inputs, performance] = off_design(stateR,componentR,designR, inputsR, performanceR,M0,alt)
%% Read Me
% This function takes in the results from the on-design analysis and
% returns the engine performance for any given flight condition

%% Calculations
state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
state(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o19';'19';'~';'o9';'9';'beta';'eptot';'ep1';'ep2'};
component = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
component(2:17,1) = {'Ram Recovery';'Inlet(Ideal,Actual)';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
design = {'Parameter','Value'};
design(2:8,1) = {'alpha';'beta';'epsilon1';'epsilon2';'PtoL';'PtoH';'h_PR'};
inputs = {'Parameter','Value'};
inputs(2:8,1) = {'Altitude (m)','Mach Number','F/mdot','Mass Flow Rate (kg/s)','SFC (kg/s/N)','Max Burner Temp (K)','Inlet Area (m^2)'};


[state, component,v0] = off_ambient(state,component,inputs,A0);
    [state,design] = off_derived_parameters(state,design,inputs,fAB);
    [state,component] = off_inlet(state,component,inputs);
    [To4] = thetabreak(state,inputs);
    To4
    state(9,3) = {To4};
    check = 0;
    alpha_track = 0;
    while check == 0
        [state,design] = off_derived_parameters(state,design,inputs,fAB);
        [state,component] = off_fan(state,component,componentR,design); %flag
        [state,component] = off_comp(state,component,design, componentR); %flag
        [state,component] = off_burner(state,component,design,inputs,fAB);
        [state,component] = off_turb(state,component,M6,A45_6); %flag
        [state,design,check,M6,alpha_track] = mixer_(state,component,design,A16_6,M6,M6A_ref,alpha_track);
        alpha_track = alpha_track+1;
    end
    
    state = off_afterburner(state,fAB);
    [state,component,performance] = off_nozzle(state,component,v0,design);
    
    
%% Functions




end