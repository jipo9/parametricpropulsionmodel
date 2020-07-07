function [T_t_over_T,P_t_over_P,MFP] = MASSFP(T_t,f,M)
% The following is a transcription of subroutine MASSFP used in AEDsys 
g_c = 25037.00; %conversion from BTU lbm to ft^2 s^2

[stag_state] = unFAIR(1,f,T_t);
v = M*stag_state.a / (1+((stag_state.gamma - 1)/2)*M^2); %ft/s
v_error = 1;

while norm(v_error) > .00001
    state.h = stag_state.h - v^2 / (2*g_c);
    [state] = unFAIR(2,f,state.h);
    v_n = M*state.a;
    if v~= 0 
        v_error = (v - v_n)/v;
    else
        v_error = (v - v_n);
    end
    v = v_n;
end
T_t_over_T = stag_state.T / state.T;
P_t_over_P = stag_state.Pr / state.Pr;
MFP = M/P_t_over_P * sqrt(state.gamma*g_c/state.R*T_t_over_T);
end

