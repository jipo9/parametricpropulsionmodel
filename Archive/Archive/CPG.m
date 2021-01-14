function [output] = CPG(Ti,Pi,gamma, prop, num)
% num = 1 is prop = Pf and output = Tf
% num = 2 is Prop = Tf and output = Pf
if num == 1
    Pf = prop;
    funT = @(Tf) (Tf/Ti)^((gamma)/(gamma-1)) - Pf/Pi;
    Tf= fzero(funT,Ti);
    output = Tf;
    
elseif num == 2
    Tf = prop;
    funP = @(Pf) (Tf/Ti)^((gamma)/(gamma-1)) - Pf/Pi;
    Pf= fzero(funP,Pi);
    output = Pf;
else
    error()
end
end

