clc
clear
close all
%Note: All values in imperial
%GET ACTUAL VALUES FOR 0'S
%On design only
%Note, we are using JP-8 fuel
%% To do
% Troubleshoot all
% Comment out nicer
% Move over to simulink
%% Inputs
%Flight Parameters
M_0 = 1.6;
T_0 = 394.1; %R
P_0 = 3.467; %psia
 
%Aircraft System Parameters
betta = 0.01; %bleed ratio
C_TOL = 0; %Low pressure shaft takeoff power coefficient
C_TOH = 0.015; %High pressure shaft takeoff power coefficient
 
%Design Limitations
    %Fuel heating valve
h_PR = 18400; %Btu/lbm, JP-8
    %Component figures of merit
        %Total pressure ratios
pi.d_max = .96; %inlet, wall friction effects
pi.b = .95; %burner
pi.M_max = .97; %mixer
pi.AB = .95; %afterburner
pi.n = .97; %nozzle
        %Polytropic efficiencies
e.f = .89; %fan
e.cL = .89; %LP compressor
e.cH = .90; %HP compressor
e.tH = .89; %HP turbine
e.tL = .90; %LP turbine
        %Component efficiency
eta.b = .999; %burner
eta.AB = .99; %afterburner
eta.mL = .995; %LP spool
eta.mH = .995; %HP spool
eta.mPL = 1; %Power takeoff - LP spool
eta.mPH = .99; %Power takeoff - HP spool
 
%Design Choices
pi.f = 3.8; %fan
pi.cL = 3.8; %LP compressor
pi.cH = 4.2105; %HP compressor, ACTUALLY AN OUTPUT???, needs to be fixed!!!
alpha = .4; %bypass ratio
T_t4 = 3200; %R, main burner exit max temp
T_t7 = 3600; %R, afterburner
M_6 = .4; %mach at mixer
P_0_over_P_9 = 1;
 
if T_t4 > 2400
    epsilon_1 = (T_t4-2400)/16000; %bypass ratio for mixer 1
    epsilon_2 = epsilon_1; %bypass ratio for mixer 2
else
    epsilon_1 = 0;
    epsilon_2 = 0;
end
g_c = 25037.00; %conversion from BTU/ lbm to ft^2/ s^2
        
%% Inlet

state_0.T = T_0;
[state_0] = unFAIR(1,0,state_0.T);
state_0.M = M_0;
state_0.P = P_0;
state_0.v = state_0.M*state_0.a;
state_t0.h = state_0.h + state_0.v^2/(2*g_c);
tau.r = state_t0.h / state_0.h;
[state_t0] = unFAIR(2,0,state_t0.h);
pi.r = state_t0.Pr/state_0.Pr;
if state_0.M<= 1
    eta.R_spec = 1; %ram recovery of MIL-E-5008B
elseif state_0.M<= 5
    eta.R_spec = 1-.075*(state_0.M-1)^1.35;   
else
    eta.R_spec = 800/(state_0.M^4 + 935);
end
pi.d = pi.d_max*eta.R_spec;
state_t2.h = state_t0.h;
state_t2.Pr = state_t0.Pr;
 
%% Fan

state_t13.Pr = state_t2.Pr*pi.f^(1/e.f);
[state_t13] = unFAIR(3,0,state_t13.Pr);
tau.f = state_t13.h / state_t2.h;
state_t13i.Pr = state_t2.Pr*pi.f;
[state_t13i] = unFAIR(3,0,state_t13i.Pr);
eta.f =  (state_t13i.h-state_t2.h) / (state_t13.h-state_t2.h);

%% Low Pressure Compressor

state_t2pt5.Pr = state_t2.Pr*pi.cL^(1/e.cL);
[state_t2pt5] = unFAIR(3,0,state_t2pt5.Pr);
tau.cL = state_t2pt5.h/state_t2.h;
state_t2pt5i.Pr = state_t2.Pr*pi.cL;
[state_t2pt5i] = unFAIR(3,0,state_t2pt5i.Pr);
eta.cL =  (state_t2pt5i.h-state_t2.h) / (state_t2pt5.h-state_t2.h);

%% High Pressure Compressor

state_t3.Pr = state_t2pt5.Pr*pi.cH^(1/e.cH);
[state_t3] = unFAIR(3,0,state_t3.Pr);
tau.cH = state_t3.h/state_t2pt5.h;
state_t3i.Pr = state_t2pt5.Pr*pi.cH;
[state_t3i] = unFAIR(3,0,state_t3i.Pr);
eta.cH =  (state_t3i.h-state_t2pt5.h) / (state_t3.h-state_t2pt5.h);
%% Main Burner

state_t4.T = T_t4;
f_4i = 0;
delta_f = 1;
while norm(delta_f)>.0001
    [state_t4] = unFAIR(1,f_4i,state_t4.T);
    f = (state_t4.h - state_t3.h)/(eta.b*h_PR - state_t4.h);
    delta_f = f_4i - f;
    f_4i = f;
end
[state_t4] = unFAIR(1,f,state_t4.T);
tau.lambda = state_t4.h / state_0.h;

%% Coolant Mixer 1

tau.m1 = ((1-betta-epsilon_1 - epsilon_2)*(1+f) + epsilon_1*tau.r*tau.cL*tau.cH/tau.lambda)/...
    ((1-betta-epsilon_1 - epsilon_2)*(1+f)+ epsilon_1);

%% High Pressure Turbine

tau.tH = 1-(tau.r*tau.cL*(tau.cH - 1)+(1+alpha)*C_TOH/eta.mPH)/...
    (eta.mH*tau.lambda*((1-betta-epsilon_1 - epsilon_2)*(1+f) + epsilon_1*tau.r*tau.cL*tau.cH/tau.lambda));
state_t4pt1.h = state_t4.h*tau.m1;
f_4pt1 = f/(1+f+epsilon_1/(1-betta-epsilon_1-epsilon_2));
 [state_t4pt1] = unFAIR(2,f_4pt1,state_t4pt1.h);
state_t4pt4.h = state_t4pt1.h*tau.tH;
 [state_t4pt4] = unFAIR(2,f_4pt1,state_t4pt4.h);
pi.tH = (state_t4pt4.Pr/state_t4pt1.Pr)^(1/e.tH);
pi.tH = .4693
state_t4pt4i.Pr = pi.tH*state_t4pt1.Pr;
 [state_t4pt4i] = unFAIR(3,f_4pt1,state_t4pt4i.Pr);
eta.tH = (state_t4pt1.h - state_t4pt4.h)/(state_t4pt1.h-state_t4pt4i.h);

%% Coolant Mixer 2

tau.m2 = ((1-betta-epsilon_1 - epsilon_2)*(1+f) + epsilon_1 + epsilon_2*tau.r*tau.cL*tau.cH/tau.lambda/tau.m1/tau.tH)/...
    ((1-betta-epsilon_1 - epsilon_2)*(1+f)+epsilon_1 + epsilon_2);
%% Low Pressure Turbine

state_t4pt5.h = state_t4pt4.h*tau.m2;
f_4pt5 = f/(1+f+(epsilon_1+ epsilon_2)/(1-betta-epsilon_1-epsilon_2));
 [state_t4pt5] = unFAIR(2,f_4pt5,state_t4pt5.h);
tau.tL = 1-(tau.r*((tau.cL - 1)+alpha*(tau.f - 1))+(1+alpha)*C_TOL/eta.mPL)/...
    (eta.mL*tau.lambda*tau.tH*((1-betta-epsilon_1 - epsilon_2)*(1+f) + (epsilon_1 + epsilon_2/tau.tH)*tau.r*tau.cL*tau.cH/tau.lambda));
state_t5.h = state_t4pt5.h*tau.tL;
 [state_t5] = unFAIR(2,f_4pt5,state_t5.h);
pi.tL = (state_t5.Pr/state_t4pt5.Pr)^(1/e.tL);
pi.tL = .4939
state_t5i.Pr = pi.tL*state_t4pt5.Pr;
 [state_t5i] = unFAIR(3,f_4pt5,state_t5i.Pr);
eta.tL = (state_t4pt5.h-state_t5.h)/(state_t4pt5.h-state_t5i.h);

%% Mixer

state_t6.h = state_t5.h;
state_t6.T = state_t5.T;
f_6 = f_4pt5;
state_t16.h = state_t13.h;
state_t16.T = state_t13.T;
state_t16.Pr = state_t13.Pr;
f_16 = 0;

%check


alpha_prime = alpha/((1-betta-epsilon_1-epsilon_2)*(1+f)+epsilon_1+epsilon_2);
f_6A = f_6/(1+alpha_prime);
state_t6A.h = (state_t6.h + alpha_prime*state_t16.h)/(1+alpha_prime);
tau.M = state_t6A.h / state_t6.h;

%100 check

P_t16_over_P_t6 = pi.f/(pi.cL*pi.cH*pi.b*pi.tH*pi.tL)
P_t16_over_P_t6 = 1.0786

[state_6] = RGCOMPR(1,state_t6.T,f_6,M_6);
state_6.M = M_6;
state_6.T = state_t6.T/state_6.T_t_over_T;
        %Not in book
         [state_6_temp] = unFAIR(1,f_6,state_6.T);
         state_6 = catstruct(state_6,state_6_temp);

         [state_t6A_temp] = unFAIR(2,f_6A,state_t6A.h);
         state_t6A = catstruct(state_t6A,state_t6A_temp);

P_t16_over_P_16 = state_6.P_t_over_P*P_t16_over_P_t6;
state_16.Pr = state_t16.Pr/P_t16_over_P_16;
 [state_16] = unFAIR(3,f_16,state_16.Pr);
state_16.v = sqrt(2*g_c*(state_t16.h - state_16.h));
state_16.M = state_16.v/state_16.a;
state_16.M = .5159
[state_16_temp] = RGCOMPR(1,state_t16.T,f_16,state_16.M);
state_16 = catstruct(state_16,state_16_temp);

A_16_over_A_6 = alpha_prime*sqrt(state_t16.T/state_t6.T)/P_t16_over_P_t6*state_6.MFP/state_16.MFP %idk where some of these came from
%alpha prime or one of the mfps
%error comes 100% from compounding pressure ratios
A_16_over_A_6 = .1844
A_6_over_A_6A = 1/(1+A_16_over_A_6);


constant = sqrt(state_6.R*state_6.T/state_6.gamma)*...
    ((1+state_6.gamma*M_6^2)+A_16_over_A_6*(1+state_16.gamma*state_16.M^2))...
    /(M_6*(1+alpha_prime));
M_6Ai = 0;
delta_M = 1;
while norm(delta_M) > .0001
    [state_6A] = RGCOMPR(1,state_t6A.T,f_6A,M_6Ai);
    state_6A.T = state_t6A.T*state_6A.T_t_over_T;
    [state_6A_temp] = unFAIR(1,f_6A,state_6A.T);
    state_6A = catstruct(state_6A,state_6A_temp);
    state_6A.M = sqrt(state_6A.R*state_6A.T/state_6A.gamma)*(1+state_6A.gamma*M_6Ai^2)/constant;
    delta_M = state_6A.M-M_6Ai;
    M_6Ai = state_6A.M;
end
M_6Ai
M_6Ai = .4331
            [state_6A] = RGCOMPR(1,state_t6A.T,f_6A,M_6Ai);
    state_6A.T = state_t6A.T*state_6A.T_t_over_T;
    [state_6A_temp] = unFAIR(1,f_6A,state_6A.T);
    state_6A = catstruct(state_6A,state_6A_temp);
    state_6A.M = sqrt(state_6A.R*state_6A.T/state_6A.gamma)*(1+state_6A.gamma*M_6Ai^2)/constant;
    delta_M = state_6A.M-M_6Ai;
    M_6Ai = state_6A.M;




pi.M_ideal = (1+alpha_prime)*sqrt(tau.M)*A_6_over_A_6A*state_6.MFP/state_6A.MFP;
pi.M = pi.M_max*pi.M_ideal
pi.M = .9771
%% Afterburner

f_7i = 0;
delta_f = 1;
state_t7.T = T_t7;
while norm(delta_f) > .0001
    [state_t7] = unFAIR(1,f_7i,state_t7.T);
    tau.lambda_AB = state_t7.h/state_0.h;
    f_AB = (1 + f*(1-betta-epsilon_1 - epsilon_2)/(1+alpha-betta))*...
        (tau.lambda_AB - tau.lambda*tau.m1*tau.tH*tau.m2*tau.tL*tau.M)/(h_PR*eta.AB/state_0.h - tau.lambda_AB);
    f_7 = f_6A + f_AB;
    delta_f = f_7i - f_7;
    f_7i = f_7;
end


%% Nozzle

f_0 = f_7;
state_t9.T = state_t7.T;
state_t9.h = state_t7.h;
state_t9.Pr = state_t7.Pr;
P_t9_over_P_9 = (P_0_over_P_9) * pi.r * pi.d * pi.cL * pi.cH * pi.b * pi.tH * pi.tL * pi.M * pi.AB * pi.n;
P_t9_over_P_9 = 12.418
state_9.Pr = state_t9.Pr/P_t9_over_P_9;
[state_9] = unFAIR(3,f_0,state_9.Pr);
state_9.v = sqrt(2*g_c*(state_t9.h-state_9.h));
state_9.M = state_9.v / state_9.a;

%% Overall Performance
g_c = 32.2; %  (lb·ft)/(lbf·s2)
F_over_mdot = state_0.a / g_c * ...
    ((1+f_0 - betta/(1+alpha))*state_9.v/state_0.a - state_0.M + ...
    (1+f_0 - betta/(1+alpha))*state_9.R/state_0.R*state_9.T/state_0.T*state_0.a/state_9.v*(1-state_0.Pr/state_9.Pr)/state_0.gamma)
%check units!!!
F_over_mdot = 110.634
S = f_0 / F_over_mdot * 3600 %lbm / hr lbm
S = 1.6878
eta.P = (2*g_c*state_0.M/state_0.a*F_over_mdot)/...
    ((1+f_0-betta/(1+alpha))*(state_9.v/state_0.a)^2 - state_0.M^2)
eta.P = .4901
g_c = 25037.00; %conversion from BTU/ lbm to ft^2/ s^2
eta.TH = (1/(2*g_c)*((1+f_0-(betta/(1+alpha)))*state_9.v^2 - state_0.v^2) + (C_TOL + C_TOH)*state_0.h)/...
    (f_0*h_PR)
eta.TH = .474
eta.o = eta.TH*eta.P;

%what is thrust??

 
 
 
 
 
 
 
 
 
