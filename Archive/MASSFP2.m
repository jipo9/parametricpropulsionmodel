function [state] = MASSFP2(state, stag_station, static_station, M)
% The following is an adaptation of subroutine MASSFP used in AEDsys 

[~,~,T0,~,~,cp0,gamma0,h0] = state{stag_station,:};
R0 = cp0 - cp0/gamma0;
a0 = sqrt(R0*gamma0*T0); %m/s
v = M*a0 / (1+((gamma0 - 1)/2)*M^2); %m/s
v_error = 1;

while norm(v_error) > .00001
    hi = h0 - v^2 /2
    if hi>0
        h = hi;
    else
        h = h/2;
    end
    
    state(static_station,2:3) = {[],[]};
    state(static_station,6:7) = {[],[]};
    state(static_station,8) = {h};
    [state] = unFAIR3(state,static_station);
    [~,~,T,~,~,cp,gamma,~] = state{static_station,:};
    R = cp - cp/gamma;
    a = sqrt(R*gamma*T); %m/s
    v_n = M*a;
    if v~= 0 
        v_error = (v - v_n)/v;
    else
        v_error = (v - v_n);
    end
    v = v_n
end
end

