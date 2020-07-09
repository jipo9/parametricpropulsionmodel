function [R,gamma] = JetA_AirMixProps(f)
% Inputs:
%   f = fuel/air ratio
% Outputs:
%   R = Individual gas constant (J/kg*K)
%   gamma = Ratio of specific heats
%   Cp = Specific heat at constant pressure
%   Cv = Specific heat at constant volume

Ru = 8.3144598; %universal gas constant J/molK

MWjeta = 158.6; %g/mol
MWair = 28.9647; %g/mol

molfuel = f/MWjeta;
molair = (1-f)/MWair;
totalmol = molfuel+molair;

Mmix = ((molfuel/totalmol)*MWjeta)+((molair/totalmol)*MWair);

R = Ru/Mmix;
R = R*1000; %J/kgK

Cpfuel = 2.01; %kJ/kgK %specific heats based on kerosene
Cpair = 1.01; %kJ/kgK
Cvfuel = 1.71; %kJ/kgK
Cvair = .72; %kJ/kgK

Cp = ((molfuel/totalmol)*Cpfuel) + ((molair/totalmol)*Cpair);
Cv = ((molfuel/totalmol)*Cvfuel) + ((molair/totalmol)*Cvair);
gamma = Cp/Cv;


end

