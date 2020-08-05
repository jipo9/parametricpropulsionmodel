clear
clc
close all

%% Temporary inputs
% will eventually take in efficiency and operation inputs
% Rest will dissapear w/ functions
inputs = 1;
state = 1;
design = 1;
component = 1;
mdot =1;



                clear
                clc
                close all

                state = {'Station','Relative Pressure', ' Temperature (K)', 'Fuel to air ratio','Mass Flow (kg/s)','Cp (J/kg-K)', 'Gamma', 'Enthalpy (J/kg)', 'Entropy (J/kg-K)','Gas Constant (m^2/s^2*K)','Relative Density(kg/m^3)','Relative Volume(s*m^3??)'};
                state(2:22,1) = {'0';'o0';'o2';'o13';'o2.5';'o3';'o3.1';'o4';'o4.1';'o4.4';'o4.5';'o5';'o6';'o6A';'o7';'o9';'9';'beta';'eptot';'ep1';'ep2'};
                component = {'Component','Pressure Ratio','Polytropic Efficieny','Enthalpy Ratio', 'Mechanical/Combustion Efficiency'};
                component(2:17,1) = {'Ram Recovery';'Inlet Actual';'Fan';'LP Compressor';'HP Compressor';'Main Burner';'Coolant Mixer 1';'HP Turbine';'HP Takeoff';'Coolant Mixer 2';'LP Turbine';'LP Takeoff';'Mixer';'Afterburner';'Nozzle';'Overall'};
                design = {'Parameter','Value'};
                design(2:8,1) = {'alpha';'beta';'epsilon1';'epsilon2';'PtoL';'PtoH';'h_PR'};
                inputs = {'Parameter','Value'};
                inputs(2:8,1) = {'Altitude (m)','Mach Number','F/mdot','Mass Flow Rate (kg/s)','SFC (kg/s/N)','Max Burner Temp (K)','Po9/P9'};

                % M0 = 1.8000;
                % alt = 40000/3.281; %altitude [m from feet]
                % To4 = 3200*.555556;
                % To7 = 3600*.555556;
                % pi_n = .97;
                % 
                % component(2:16,2) = {5.7458; .9067; 3.0054; 3.0054; 4.7208; .95;    [];    .4231; []; [];    .5023; []; .9735; .95; pi_n};
                % component(2:16,4) = {1.6480; [];    1.4259; 1.4259; 1.6377;  []; .9673;    .8381; []; .9731; .8667; []; .8268;  []; []};
                % 
                % mdot = 188.72*0.45359237;
                % alpha = .530;
                % beta = .01;
                % ep1 = .05;
                % ep2 = .05;
                % f = .02875;
                % fAB = .03371;

                M0 = 1.451;
                alt = 36000/3.281; %altitude [m from feet]
                To4 = 3200*.555556;
                To7 = 3600*.555556;
                pi_n = .97;

                component(2:16,2) = {3.4211; .9354; 3.9000; 3.9000; 5.1282; .95; [];    .4231; []; [];    .4831; []; .9779; .95; pi_n};
                component(2:16,4) = {1.4211; [];    1.5479; 1.5479; 1.6803;  []; .9673; .8381; []; .9731; .8598; []; .8404;  []; []};

                mdot = 188.72*0.45359237;
                alpha = .449;
                beta = .01;
                ep1 = .05;
                ep2 = .05;
                f = .03070;
                fAB = .03353;

                state{2,5} = mdot;
                state{9,4} = f;
                state{9,3} = To4;
                state{16,4} = fAB;
                state{16,3} = To7;

                inputs{2,2} = alt;
                inputs{3,2} = M0;

                design{2,2} = alpha;
                design{3,2} = beta;
                design{4,2} = ep1;
                design{5,2} = ep2;


                [state,design] = derived_parameters_performance(state,inputs,design,component);
                [state,component] = a(state,component,alt,To4);
                
                
                
                A16_6 = .2715;
                M6 = .4;
                M6A_ref = .4188;
                alt = inputs{2,2};
                M0 = inputs{3,2};
                [T0, ~, ~, ~] = atmosisa(alt); %obtain standard atmospheric conditions
                state(2,3) = {T0};
                [state] = unFAIR3(state,2);
                [~,~,T0,~,~,cp0,gamma0,~,~] = state{2,:};
                R0 = cp0 - cp0/gamma0;
                a0 = sqrt(R0*gamma0*T0); %[m/s]
                v0 = M0*a0; %[m/s]
               
%% Calcs
[performance,inputs,state,design,component] = offdesign(inputs,state,design,component,A16_6,M6,M6A_ref,v0);
performance


%% Functions
function [performance,inputs,state,design,component] = offdesign(inputs,state,design,component,A16_6,M6,M6A_ref,v0);
%     derived_parameters_parametric(X);
%     ambient(X);
%     inlet(X);
    check = 0;
    while check == 0
%         derived_parameters_parametric(X);
%         fan(X);
%         LPcomp(X);
%         HPcomp(X);
%         burner(X);
%         HPturbine(X); %w/ mixer
%         LPturbine(X); %w/ mixer
        [state,design,check,M6] = mixer(state,component,design,A16_6,M6,M6A_ref);
    end
    state = afterburner(state);
    [state,component,performance] = nozzle(state,component,v0,design);
end

function state = afterburner(state)
    %Input Tt7 and f_AB in initial conditions
    [state] = unFAIR3(state,16);
    
    %This is temporary and solely for reference purposes
end

function [state,component,performance] = nozzle(state,component,v0,design)
%Take in inputs    
alpha = design{2,2};
beta = design{3,2};
P0_P9 = 1;
% h_PR = design{8,2};
% PtoL = design{6,2};
% PtoH = design{7,2};
[~,pi_r,pi_d,~,pi_cL,pi_cH,pi_b,~,pi_tH,~,~,pi_tL,~,pi_M,pi_AB,pi_n,~] = component{:,2};

%Calculate pressure drop across nozzle
Pro7 = state{16,2};
pin = component{16,2};
Pro9 = Pro7*pin;
state(17,2) = {Pro9};
[state] = unFAIR3(state,17);

%Calculate static conditions of exhaust
Po9_P0 = pi_r*pi_d*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_M*pi_AB*pi_n;
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
%Actual is 110.83*9.806655 = 1.0869e+03
F  = F_mdot*mdot0;
S = f_0 / F_mdot;
%Actual is 1.6941*2.8325e-05 =  4.7985e-05

% eta_TH = ((v0^2/2*((1+f_0-(beta/(1+alpha)))*(v9/v0)^2 - 1) + (PtoL + PtoH)/mdot0))/...
%     (f_0*h_PR);
eta_P = (2*F_mdot/v0)/...
    ((1+f_0-beta/(1+alpha))*(v9/v0)^2 - 1);
% eta_o = eta_TH*eta_P;
performance(1,:) = {'Thrust','Specific Fuel Consumption','Propulsive Efficiency','Thermal Efficiency','Overall Efficiency','Exhaust Mach'};
performance(2,:) = {F_mdot,S,0,eta_P,0,M9};
end

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

function [state,design] = derived_parameters_performance(state,inputs,design,component)
%Derived parameters for performance model, w/ changing values

mdot = state{2,5};
alpha = design{2,2};
beta = design{3,2};
ep1 = design{4,2};
ep2 = design{5,2};
f = state{9,4};
fAB = state{16,4};
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

