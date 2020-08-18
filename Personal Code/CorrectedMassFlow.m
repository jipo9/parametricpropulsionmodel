function [mdotc] = CorrectedMassFlow(state,altR,alt,component)
% This function calculates the corrected mass flow rates at any of the
% compressors/fan

[T0_std, ~, P0_std, ~] = atmosisa(altR); %obtain standard atmospheric conditions at on design
[~, ~, P0, ~] = atmosisa(alt); %obtain standard atmospheric conditions at off design
[~,pi_r,pi_d,pi_f,pi_cL,pi_cH,pi_b,~,pi_tH,~,~,pi_tL,~,~,~,pi_n,~] = component{:,2};

% Fan
pi = pi_r*pi_d*pi_f*P0/P0_std;
T = state{5,3};
tau = T/T0_std;
mdot = state{5,5};
mdotc_f = mdot * sqrt(tau) / (pi);

% LP Compressor
pi = pi_r*pi_d*pi_cL*P0/P0_std;
T = state{6,3};
tau = T/T0_std;
mdot = state{6,5};
mdotc_cL = mdot * sqrt(tau) / (pi);

% HP Compressor
pi = pi_r*pi_d*pi_cL*pi_cH*P0/P0_std;
T = state{6,3};
tau = T/T0_std;
mdot = state{7,5};
mdotc_cH = mdot * sqrt(tau) / (pi);

mdotc = [mdotc_f,mdotc_cL,mdotc_cH];
end

