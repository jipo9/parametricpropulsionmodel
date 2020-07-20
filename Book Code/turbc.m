function [state] = turbc(state,M4,M45,A4,A45,To45R,eta_t)
%take in M4, M45, A5,A45,To45R, eta_t
%fair(f,To4)
ho3 = state{7,8};



To4 = state{9,3};
f = state{9,4};
mdot4 = state{9,5};

[To4_T4,Po45_P4,MFP4] = MASSFP(To4,f,M4);


%fair(f=0,To3)

% mdotf = f*(1-betta-ep1-ep2);
% mdot4 = (1+f)*(1-betta-ep1-ep2);
% mdot41 = (1+f)*(1-betta-ep1-ep2) + ep1;
% mdot45 = (1+f)*(1-betta-ep1-ep2) + ep1 + ep2;
% 
% f41 = mdotf/(mdot41 - mdotf);
% f45 = mdotf/(mdot45 - mdotf);


hoi = state{9,8}; %???
mdot41 = state{10,5};
mdotep1 = state{21,5};
ho41 = (mdot4*hoi + mdotep1*ho3) / (mdot41);
state(10,8) = {ho41};
[state] = unFAIR3(state,10);
Pro41 = state{10,2};

f45 = state{12,4};
To45 = To45R;
mdotep2 = state{21,5};
mdot45 = state{12,5};

error = 1;

while error >.01
    [To45_T45,Po45_P45,MFP45] = MASSFP(To45,f45,M45);
    pi_tH = (MFP4*A4/(mdot4 * sqrt(To4))) /...
        (MFP45*A45/(mdot45 * sqrt(To45)));
    Pro44i = pi_tH*Pro41;
    
    
    state(11,2) = {Pro44i};
    state(11,3) = {[]};
    state(11,8) = {[]};
    [state] = unFAIR3(state,12);
    ho44i = state{11,8};
    
    state(10,3) = {To45};
    state(10,2) = {[]};
    state(10,8) = {[]};
    [state] = unFAIR3(state,10);
    
    ho44 = ho41 - eta_t * (ho41 - ho44i);
    tau_tH = ho44/ho41;
    ho45 = (mdot41*ho44 + mdotep2*ho3)/(mdot45);
    state(10,2) = {ho45};
    state(10,3) = {[]};
    state(10,8) = {[]};
    [state] = unFAIR3(state,10);
    
    To45n = state{10,3};
    error = norm(To45 - To45n);
    To45 = To45n;
end
end

