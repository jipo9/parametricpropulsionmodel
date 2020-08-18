function [mdot_cor] = CorrectedMassFlow(state,altR,alt,num)
% This function calculates the corrected mass flow rates at any of the
% compressors/fan

[T0_std, ~, P0_std, ~] = atmosisa(altR); %obtain standard atmospheric conditions at SL
[~, ~, P0, ~] = atmosisa(alt); %obtain standard atmospheric conditions at SL


[~,pi_r,pi_d,pi_f,pi_cL,pi_cH,pi_b,~,pi_tH,~,~,pi_tL,~,~,~,pi_n,~] = component{:,2};
if num == 1 %fan
    pi = pi_r*pi_d*pi_f*P0/P0_std;
    T = state{5,3};
    tau = T/T0_std;
    mdot = state{5,5};
elseif num == 2 %LP compressor
    pi = pi_r*pi_d*pi_cL*P0/P0_std;
    T = state{6,3};
    tau = T/T0_std;
    mdot = state{6,5};
elseif num == 3 %HP Compressor
    pi = pi_r*pi_d*pi_cL*pi_cH*P0/P0_std;
    T = state{6,3};
    tau = T/T0_std;
    mdot = state{7,5};
else
    error('Return a valid number. 1 for fan, 2 for LP compressor, 3 for HP compressor')
end
        
mdot_cor = mdot * sqrt(tau) / (pi);
end

