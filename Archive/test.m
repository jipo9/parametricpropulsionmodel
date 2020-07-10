clear
clc

state = {'Station','Pressure', ' Temperature', 'Fuel to air ratio','Cp', 'Gamma','Mass Flow','Pressure Ratio','Enthalpy Ratio','Polytropic Efficiency'};
state(2,:) = {'1',2,3,4,5,0,0,0,0,0};
state(3,3) = {500};
state(3,4) = {0};
state(3,1) = {'2'};

[state] = unFAIR3(state)