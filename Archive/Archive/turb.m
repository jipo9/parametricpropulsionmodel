function [a] = turb_R(state,M45,M6)
To45 = state{12,3};
f45 = state{12,4};
[To45_T45,Po45_P45,MFP45] = MASSFP(To45,f45,M45);
M6_M45 = M6/M45;

while T_error > .01
    [To6_T6,Po6_P6,MFP6] = MASSFP(To6,f45,M6);
   
    state(13,2) = {[]};
    state(13,8) = {[]};
    [state] = unFAIR3(state,13);
    
    gamma45 = state{12,7};
    R45 = state{12,10};
    mdot45 = state{12,5};
    To45 = state{12,3};
    Pro45 = state{12,2};
    
    gamma6 = state{13,7};
    R6 = state{13,10};
    mdot6 = state{13,5};
    To6 = state{13,3};
    Pro6 = state{13,2};
    
    fun = @(A6_A45)...
        M45*sqrt(gamma45/R45)*(1+(gamma45-1)/2 * M45^2)^((gamma45 + 1)/(2*(1-gamma45)))...
        /...
        M45*M6_M45*sqrt(gamma6/R6)*(1+(gamma6-1)/2 * (M45*M6_M45)^2)^((gamma6 + 1)/(2*(1-gamma6)))...
        -...
        (mdot45*sqrt(To45)/Pro45)/(mdot6*sqrt(To6)/Pro6)    *    A6_A45;
    A6_A45 = fzero(fun,M6_M45);

    
    pi_tL = MFP45/sqrt(To45)/(MFP6/sqrt(To6))    / A6_A45;
    
    Pro6 = pi_tL*Pro45;
%     %fair(f,Proei)
%     hoe = hoi - eta_t*(hoi - hoei);
%     tau_t = hoe/hoi;
%     %fair(f,hoe) - >Toen
%     T_error = abs(Toe - Toen)
%     Toe = Toen;
end
a = 1;
end

