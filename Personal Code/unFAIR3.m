function [state] = unFAIR3(state,station)
% The following is an adaptation of the FAIR method developed in AEDsys
% system and the subsequent adaptation of the unFAIR method developed by
% Captain Andrus in Turbofan Analysis thesis.

% This method takes in fuel to air ratio and the option of 1) pressure,
% 2)temperature, or 3) enthalpy and calculates thermodynamic states
% based on NASA Glenn Thermochemical data


%% Initialize function
%Check if fuel air ratio has an impact on results
f = state{station, 4};
if f< 1E-9
    f = 0;
end

T_min = 0;
T_max = 100000;
error = 100;
ref_Pressure = 7.5205 * 6894.76; %Pa from Psi

%If temperature is being input, proceed. Otherwise, use bisection method to
%get an estimate of temperature within .1% accuracy of the input value before proceeding
if size(state{station,2},1) == 1 
% if state{station,2} ~= 0 
    P = state{station,2};
    while norm(error) > .00001
        T= (T_min + T_max)/2;
        state_i = state;
        state_i(station,3) = {T* .5556};
        state_i(station,2) = state_i(station,8);
        [state_i] = unFAIR3(state_i,station);
        error = (state_i{station,2} - P)/P;
        if error<0
            T_min = T;
        else
            T_max = T;
        end
    end    
elseif size(state{station,8},1) == 1 
    h = state{station,8};
    while norm(error) > .00001
        T= (T_min + T_max)/2;
        state_i = state;
        state_i(station,8) = state_i(station,2);
        state_i(station,3) = {T* .5556};
        [state_i] = unFAIR3(state_i,station);
        h_i = state_i{station,8};
        
        %The following two lines are to detect if a temperature value too
        %low is input, which would cause this code to break
%         if h_i <0
%                 h_i = 0;
%         end
        error = (h_i - h)/h;
%         if isnan(error)
%             error  = -1;
%         end
%         
        if error<0
            T_min = T;
        else
            T_max = T;
        end
        
    end    
% elseif size(state{end,3},1) == 1
elseif state{station,3} ~= 0
    T = state{station,3} / .5556; % R;
else
    error()
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
R =  1.9857117 / (28.97-f*0.946186); %BTU/lbm R
cp = ( air.cp + f * vitiated.cp ) / (1+f); %BTU/lbm R
h =  ( air.h + f * vitiated.h ) / (1+f); %BTU/lbm
phi =( air.phi + f * vitiated.phi ) / (1+f); %BTU/lbm

phi_ref = 1.578420959; %BTU/ lbm R (phi at 492 R)
Pr = exp((phi - phi_ref) / R); %unitless
gamma = cp / (cp-R); %unitless
% a = sqrt(gamma*R*T*g_c); %ft/s

% Comversion to metric
T = T* .5556; % K
% P = Pr*ref_Pressure; %Pa
P = Pr; %temporary
cp = cp* 4186.8; %J/kg K
h = h* 4186.8*.5556; %J/kg
phi =( air.phi + f * vitiated.phi ) / (1+f)* 4186.8; %J/kg K


state(station,2:3) = {P,T};
state(station,6:9)  = {cp,gamma,h,phi};
end

