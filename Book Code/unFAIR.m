function [state] = unFAIR(num,f,prop)
% The following is an adaptation of the FAIR method developed in AEDsys
% system and the subsequent adaptation of the unFAIR method developed by
% Captain Andrus in Turbofan Analysis thesis.

% This method takes in fuel to air ratio and the option of 1) temperature,
% 2)enthalpy, or 3) reduced pressure and calculates thermodynamic states
% based on NASA Glenn Thermochemical data

%% Initialize function
%Check if fuel air ratio has an impact on results
if f< 1E-9
    f = 0;
end
g_c = 25037.00; %conversion from BTU lbm to ft^2 s^2

T_min = 0;
T_max = 10000;
error = 100;

%If temperature is being input, proceed. Otherwise, use bisection method to
%get an estimate of temperature within .1% accuracy of the input value before proceeding
if num == 1
    T = prop;
elseif num == 2
    h = prop;
    while norm(error) > .00001
        T= (T_min + T_max)/2;
        [state_i] = unFAIR(1,f,T);
        error = (state_i.h - h)/h;
        if error<0
            T_min = T;
        else
            T_max = T;
        end
    end
elseif num == 3
    Pr = prop;
    while norm(error) > .00001
        T= (T_min + T_max)/2;
        [state_i] = unFAIR(1,f,T);
        error = (state_i.Pr - Pr)/Pr;
        if error<0
            T_min = T;
        else
            T_max = T;
        end
    end
else
    warning('Invalid number input')
end

%% Input data
%Input data for air
air.A0 = 2.5020051E-01;
air.A1 =-5.1536879E-05;
air.A2 = 6.5519486E-08;
air.A3 =-6.7178376E-12;
air.A4 =-1.5128259E-14;
air.A5 = 7.6215767E-18;
air.A6 =-1.4526770E-21;
air.A7 = 1.0115540E-25;
air.h_ref =-1.7558886; %BTU/lbm
air.phi_r1 = 0.0454323; %BTU/lbm

%Input data for products of combustion
vitiated.A0 = 7.3816638E-02;
vitiated.A1 = 1.2258630E-03;
vitiated.A2 =-1.3771901E-06;
vitiated.A3 = 9.9686793E-10;
vitiated.A4 =-4.2051104E-13;
vitiated.A5 = 1.0212913E-16;
vitiated.A6 =-1.3335668E-20;
vitiated.A7 = 7.2678710E-25;
vitiated.h_ref = 30.58153; %BTU/lbm
vitiated.phi_r2 = 0.6483398; %BTU/lbm



%% Calculate Thermodymic State

%Calculate parameters for air for air
air.cp = air.A0 + air.A1*T + air.A2*T^2 + air.A3*T^3 + air.A4*T^4 + air.A5*T^5 + air.A6*T^6 + air.A7*T^7;
air.h = air.h_ref + air.A0*T + air.A1/2*T^2 + air.A2/3*T^3 + air.A3/4*T^4 + air.A4/5*T^5 + air.A5/6*T^6 + air.A6/7*T^7 + air.A7/8*T^8;
air.phi = air.phi_r1 + air.A0*log(T) + air.A1*T + air.A2/2*T^2 + air.A3/3*T^3 + air.A4/4*T^4 + air.A5/5*T^5 + air.A6/6*T^6 + air.A7/7*T^7;

%Calculate parameters for products of combustion
vitiated.cp = vitiated.A0 + vitiated.A1*T + vitiated.A2*T^2 + vitiated.A3*T^3 + vitiated.A4*T^4 + vitiated.A5*T^5 + vitiated.A6*T^6 + vitiated.A7*T^7;
vitiated.h = vitiated.h_ref + vitiated.A0*T + vitiated.A1/2*T^2 + vitiated.A2/3*T^3 + vitiated.A3/4*T^4 + vitiated.A4/5*T^5 + vitiated.A5/6*T^6 + vitiated.A6/7*T^7 + vitiated.A7/8*T^8;
vitiated.phi = vitiated.phi_r2 + vitiated.A0*log(T) + vitiated.A1*T + vitiated.A2/2*T^2 + vitiated.A3/3*T^3 + vitiated.A4/4*T^4 + vitiated.A5/5*T^5 + vitiated.A6/6*T^6 + vitiated.A7/7*T^7;

%Calculate thermo state
state.R =  1.9857117 / (28.97-f*0.946186); %BTU/lbm R
state.cp = ( air.cp + f * vitiated.cp ) / (1+f); %BTU/lbm R
state.h =  ( air.h + f * vitiated.h ) / (1+f); %BTU/lbm
state.phi =( air.phi + f * vitiated.phi ) / (1+f); %BTU/lbm

phi_ref = 1.578420959; %BTU/ lbm R (phi at 492 R)
state.Pr = exp((state.phi - phi_ref) / state.R); %unitless
state.gamma = state.cp / (state.cp-state.R); %unitless
state.a = sqrt(state.gamma*state.R*T*g_c); %ft/s
state.T = T; % R
end

