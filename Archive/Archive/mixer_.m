function [state,design,check,M6,alpha_track] = mixer_(state,component,design,A16_6,M6,M6A_ref,alpha_track)
%check = 0 if bad and 1 if good
a_dampner = .9;
n = floor(alpha_track/20);
a_dampner = a_dampner * .8^n;

if a_dampner < .1
    a_dampner = .1;
end

[alpha,beta,ep1,ep2] = design{2:5,2};
[pi_f, pi_cL, pi_cH, pi_b, ~, pi_tH, ~, ~, pi_tL] = component{4:12,2};
[~,~,To16,~,~,~,gamma16,~,~,R16,~] = state{5,:};
[~,~,~,f,~,~,~,~,~,~,~] = state{9,:};
[~,~,To6,~,~,~,gamma6,~,~,R6,~] = state{14,:};

Po16_Po6 = pi_f / (pi_cL * pi_cH * pi_b * pi_tH * pi_tL);
%         Po16_Po6 = 1/1.0042;

MFP6 = MFP2(M6, gamma6, R6);

% Calculate M16 and MFP16
[M16] = Kutta_mach(gamma16,M6,gamma6,Po16_Po6);

if ~isreal(M16)
    M16 = 0.1;
end

MFP16 = MFP2(M16, gamma16, R16);

% Find bypass ratios and quantify error
[alpha_prime,alpha_i] = bypass_ratio(MFP16,To16,MFP6,Po16_Po6,To6,A16_6,beta,ep1,ep2,f);

alpha_err = norm((alpha - alpha_i)/alpha);
if alpha_err > .01 %changed from .001 for speed
    design{2,2} = alpha - a_dampner*(alpha - alpha_i);
    check = 0;
    
else
    [M6A_i,state] = mixer_state(state,alpha_prime,A16_6,M16,M6);
    M6A_err = norm((M6A_ref - M6A_i)/M6A_ref);
    if M6A_err > .001
        M6 = M6 + .1*(M6A_ref - M6A_i); %maybe needs a damper in there? (IE multiply difference by .9)
        check = 0;
        alpha_track = 0;
    else
        check = 1;
    end
end
        


function [M16] = Kutta_mach(gamma16,M6,gamma6,Po16_Po6)
%calculates the mach if kutta condition is satisfied between the two states
P_Pt6 = pressure(M6,gamma6);
P_Pt16 = P_Pt6*Po16_Po6;
M16 = sqrt((2/(gamma16 - 1)) * (P_Pt16 ^ ((gamma16 - 1)/gamma16)  - 1));
end

function [alpha_prime,alpha] = bypass_ratio(MFP16,Tt16,MFP6,Pt6_Pt16,Tt6,A16_6,beta,ep1,ep2,f)
%Calculates bypass ratio at mixer
alpha_prime = 1/Pt6_Pt16 * MFP16/MFP6 * sqrt(Tt6/Tt16) * A16_6;
% alpha_prime = MFP16/MFP6 * A16_6 * 1.65;
alpha = alpha_prime * ((1-beta-ep1-ep2)*(1+f) + ep1 +ep2);
end

function [MFP] = MFP1(mdot, Tt, Pt, A)
%Calculates mass flow, stagnation temp and pressure, and Area dependent MFP
MFP = mdot*sqrt(Tt) / (Pt*A); %kg/s*K^.5 / Pa*m^2 = s*K^.5 / m
end

function [MFP] = MFP2(M, gamma, R)
%Calculates Mach, gamma, and R dependent MFP
[P_Pt] = pressure(M,gamma);
[T_Tt] = temperature(M,gamma);
MFP = M*sqrt(gamma/R)/sqrt(T_Tt)*sqrt(P_Pt); %sqrt(s^2*K / m^2) = s*K^.5 / m
end

function [P_Pt] = pressure(M,gamma)
%Calculates static over stagnation pressure for isentropic compressible
%gas
P_Pt = (1 - (gamma - 1)/2*M^2)^(-gamma/(gamma-1));
end

function [T_Tt] = temperature(M,gamma)
%Calculates static over stagnation tempurature for isentropic compressible
%gas
T_Tt = (1 - (gamma - 1)/2*M^2)^-1;
end
        
function [M6A,state] = mixer_state(state,alpha_prime,A16_6,M16,M6)
%Calculates state and ideal pressure ratio of the mixer
[~,~,To16,~,~,cp16,gamma16,ho16,~,R16,~] = state{5,:};
[~,~,To6,~,mdot6,cp6,gamma6,ho6,~,R6,~] = state{14,:};

A6_6A = 1/(1 + A16_6);

ho6A = (ho6 + alpha_prime*ho16)/(1+alpha_prime);
state(15,2) = {[]};
state(15,3) = {[]};
state(15,8) = {ho6A};
[state] = unFAIR3(state,15);
[~,~,To6A,~,mdot6A,cp6A,gamma6A,ho6A,~,R6A,~] = state{15,:};

[T_Tt16] = temperature(M16,gamma16);
[T_Tt6] = temperature(M6,gamma6);
T16 = To16 * T_Tt16;
T6 = To6 * T_Tt6;

const = 1/(1+alpha_prime) * ...
    (sqrt(R6*T6/gamma6)*(1+gamma6*M6^2)/(M6)   +   alpha_prime*sqrt(R16*T16/gamma16)*(1+gamma16*M16^2)/(M16));

M6Ai = M6;
error = 1;

while error > .0001
    [T_Tt6A] = temperature(M6Ai,gamma6A);
    T6A = T_Tt6A * To6A;
    M6A = sqrt(R6A*T6A/gamma6A)*(1+gamma6A*M6Ai^2)/(const);
    error = norm((M6Ai - M6A)/M6Ai);
    M6Ai = M6A;
end

% [MFP6] = MFP2(M6, gamma6, R6);
% [MFP6A] = MFP2(M6A, gamma6A, R6A);
% 
% pi = (1+alpha_prime) * sqrt(To6A/To6) * A6_6A * MFP6 / MFP6A;
end

end        