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
    state.TtTratio = T_t_over_T;
    state.PtPratio = P_t_over_P;
    state.MFP = MFP;
    state.M = M;
end

if item == 2 || 3
    [statet] = unFAIR(f,1,Tt);
    if item == 2
        TtTratio = prop;
        T = statet.T/TtTratio;
        [statex] = unFAIR(f,1,T);
    elseif item == 3
        PtPratio = item;
        Pr = statet.Pr / PtPratio;
        [statex] = unFAIR(f,3,Pr);
    end
    
    gc = 31.174; %lbm*ft/lbf*s^2
    Vsq = (2*gc*(statet.h-statex.h));
    if Vsq < 0
        state.M = 0;
        %T = statet.T; %why is this here
    else
        state.M = (sqrt(Vsq))/a;
    end
    [T_t_over_T,P_t_over_P,MFP] = MASSFP(Tt,f,state.M);
    state.TtTratio = T_t_over_T;
    state.PtPratio = P_t_over_P;
    state.MFP = MFP;
    
end

if item == 4 || 5
    MFP = prop;
    if item == 4
        state.M = 2;
    elseif item == 5
        state.M = .5;
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
                state.TtTratio = T_t_over_T;
                state.PtPratio = P_t_over_P;
                state.MFP = MFP;
            end
        end
        
    end
    
    %%CHANGEETST
    
end

