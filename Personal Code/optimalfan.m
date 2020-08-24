function [pi_f] = optimalfan(alpha)
%This function calculates the optimal fan pressure ratio
pif_VS_alpha = ...
    [1.1, 1.15, 1.2, 1.25, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0,  2.25, 2.5; 
     30,  20,   17,  14.5, 12,  9.5, 7.5, 6,   5,   4.5, 3.75, 3.05,    3];

 
for ii = 1:size(pif_VS_alpha,2)-1
   alpha1 = pif_VS_alpha(2,ii);
   alpha2 = pif_VS_alpha(2,ii+1);
   if alpha1 >= alpha
       if alpha2 <= alpha
           pif1 = pif_VS_alpha(1,ii);
           pif2 = pif_VS_alpha(1,ii+1);
           m = (pif2-pif1) / (alpha2-alpha1);
           
           pi_f =  pif1 + m*(alpha-alpha1)
       end
   end
end

end

