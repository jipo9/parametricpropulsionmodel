function [cp,gamma] = unFAIR2(f)
% The following is an adaptation of the FAIR method developed in AEDsys
% system and the subsequent adaptation of the unFAIR method developed by
% Captain Andrus in Turbofan Analysis thesis.

% This method takes in fuel to air ratio and assumes 400 R

% Units are in metric


%% Initialize function
%Check if fuel air ratio has an impact on results
if f< 1E-9
    f = 0;
end
T = 400; %Rankine
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
R =  1.9857117 / (28.97-f*0.946186) * 4186.8; % J/kg K from BTU/lbm R
cp = ( air.cp + f * vitiated.cp ) / (1+f) * 4186.8; % J/kg K from BTU/lbm R
gamma = cp / (cp-R); %unitless
end

