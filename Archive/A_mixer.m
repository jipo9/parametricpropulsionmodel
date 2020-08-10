function [A16_6,P_ref] = A_mixer(state,component,design,inputs,M6_est)
[alpha,beta,ep1,ep2] = design{2:5,2};
[pi_f, pi_cL, pi_cH, pi_b, ~, pi_tH, ~, ~, pi_tL] = component{4:12,2};
[~,~,To16,~,~,~,gamma16,~,~,R16,~] = state{5,:};
[~,~,~,f,~,~,~,~,~,~,~] = state{9,:};
[~,~,To6,~,~,~,gamma6,~,~,R6,~] = state{14,:};

Po16_Po6 = pi_f / (pi_cL * pi_cH * pi_b * pi_tH * pi_tL);
% Po16_Po6 = pi_b;
% Po16_Po6 = 1.0042;
% Po16_Po6 = 1.0786;

upper_bound = 100;
lower_bound = 0;
error = 1;
while norm(error) > .0001
    A16_6 = (upper_bound + lower_bound)/2;
    MFP6 = MFP2(M6_est, gamma6, R6);
    % Calculate M16 and MFP16
    [M16] = Kutta_mach(gamma16,M6_est,gamma6,Po16_Po6);
    MFP16 = MFP2(M16, gamma16, R16);
    % Find bypass ratios and quantify error
    [alpha_prime,alpha_i] = bypass_ratio(MFP16,To16,MFP6,Po16_Po6,To6,A16_6,beta,ep1,ep2,f);
    error = (alpha_i-alpha)/alpha;
    if error > 0
        upper_bound = A16_6;
    else
        lower_bound = A16_6;
    end
end


alt = inputs{2,2};
[~, ~, P0, ~] = atmosisa(alt); %obtain standard atmospheric conditions
[~,pi_r,pi_d,~,pi_cL,pi_cH,pi_b,~,pi_tH,~,~,pi_tL,~,pi_M,~,~,~] = component{:,2};
Po6A_P0 = pi_r*pi_d(2)*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_M;
P_ref = Po6A_P0 * P0;
    
% [M6A_i,state] = mixer_state(state,alpha_prime,A16_6,M16,M6_est);

%% Functions

function [M16] = Kutta_mach(gamma16,M6,gamma6,pi6_16)
%calculates the mach if kutta condition is satisfied between the two states
P_Pt6 = pressure(M6,gamma6);
P_Pt16 = P_Pt6*pi6_16;
M16 = sqrt((2/(gamma16 - 1)) * (P_Pt16 ^ ((gamma16 - 1)/gamma16)  - 1));
end

function [alpha_prime,alpha] = bypass_ratio(MFP16,Tt16,MFP6,Pt16_Pt6,Tt6,A16_6,beta,ep1,ep2,f)
%Calculates bypass ratio at mixer
alpha_prime = Pt16_Pt6 * MFP16/MFP6 * sqrt(Tt6/Tt16) * A16_6;
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

