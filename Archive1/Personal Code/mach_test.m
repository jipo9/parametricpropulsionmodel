 clear
 clc
 close all

 
 
 %% 
 M = linspace(.5,1);
 gamma = 1.4;
 for ii = 1 : size(M,2)
     Mi = M(ii);
     [P_Pt] = pressure(Mi,gamma);
     P(ii) = P_Pt;
 end
 
 %%
P9_P0 = 1;
Po9_Po0 = 11;

M = 1.8;
alt = 40000/3.281;

%% State 0
[T0, P0, ~, ~] = atmosisa(alt);
state(2,3) = {T0};
[state] = unFAIR3(state,2);
gamma0 = state{2,4};
Pr0 = state{2,2};


%% State 0o
To0 = T0*(1+((M0^2)*((gamma0-1)/2))); %find total temperature using isentropic
state(3,3) = {To0};
[state] = unFAIR3(state,3);
Pro0 = state{3,2};

Po0_P0 = Pro0/Pr0;
%% State 9
f = .02;
Pro9 = Pro0 * Po9_Po0;
state(17,2) = {Pro9};
[state] = unFAIR3(state,17);
Pro0 = state{17,2};


%% State 9o
[~,~,T9,~,~,cp9,gamma9,h9] = state{18,:};
R9 = cp9 - cp9/gamma9;
a9 = sqrt(R9*gamma9*T9); %m/s
v9 = sqrt(2*(ho9-h9));
M9 = v9 / a9;


function [P_Pt] = pressure(M,gamma)
%Calculates static over stagnation pressure for isentropic compressible
%gas
P_Pt = (1 - (gamma - 1)/2*M^2)^(-gamma/(gamma-1));
end

