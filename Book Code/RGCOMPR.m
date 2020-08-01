function [state] = RGCOMPR(item,Tt,f,prop)
%The following is a transcription of subroutine RGCOMPR used in AEDsys
%Inputs:
%   item = specified given variable
%       1 = M
%       2 = Tt/T ratio
%       3 = Pt/P ratio
%       4 = MFP with M<=1
%       5 = MFP with M>1
%   Tt = Total Temperature
%   f = Fuel/Air Ratio
%   prop = Given property as specified above
%Outputs:
%   state = thermodynamic state vector including...
%       M, Tt/T , Pt/P, MFP

if item == 1
    M = prop;
    [T_t_over_T,P_t_over_P,MFP] = MASSFP(Tt,f,M);
    state.T_t_over_T = T_t_over_T;
    state.P_t_over_P = P_t_over_P;
    state.MFP = MFP;
    state.M = M;
end

if item == 2 || item == 3
    [statet] = unFAIR(1,f,Tt);
    if item == 2
        T_t_over_T = prop;
        T = statet.T/T_t_over_T;
        [statex] = unFAIR(1,f,T);
    elseif item == 3
        P_t_over_P = item;
        Pr = statet.Pr / P_t_over_P;
        [statex] = unFAIR(3,f,Pr);
    end
    
    g_c = 25037.00; %conversion from BTU lbm to ft^2 s^2
    Vsq = (2*g_c*(statet.h-statex.h));
    if Vsq < 0
        state.M = 0;
        %T = statet.T; %why is this here
    else
        a = statex.a;
        state.M = (sqrt(Vsq))/a;
    end
    [T_t_over_T,P_t_over_P,MFP] = MASSFP(Tt,f,state.M);
    state.T_t_over_T = T_t_over_T;
    state.P_t_over_P = P_t_over_P;
    state.MFP = MFP;
    
end

if item == 4 || 5
    MFP = prop;
    if item == 4
        state.M = 2;
    elseif item == 5
        state.M = .5;
    end
    dM = .1;
    [~,~,MFP0] = MASSFP(Tt,f,state.M);
    MFPerr = 1;
    while MFPerr>.00001
        state.M = state.M+dM;
        [T_t_over_T,P_t_over_P,MFPn] = MASSFP(Tt,f,state.M);
        MFPerr = abs(MFPn-MFP0);
        if MFPerr > .00001
            dM = ((MFP-MFPn)/(MFPn-MFP0))*dM;
            MFP0 = MFPn;
        else
            state.M = state.M;
            state.T_t_over_T = T_t_over_T;
            state.P_t_over_P = P_t_over_P;
            state.MFP = MFP;
        end
    end
    
end


