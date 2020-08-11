function [pi_dmax, e_c, e_f, pi_b, eta_b, e_t, pi_n, To4max] = LOT(year)
%Based on Table 4.4 in Mattingly

if year > 1945 && year <= 1965
    pi_dmax = .90;
    e_c = .80;
    e_f = .78;
    pi_b = .90;
    eta_b = .88;
    e_t = .80;
    pi_n = .95;
    To4max = 1110; %[K]     
elseif year > 1965 && year <= 1985    
    pi_dmax = .95;
    e_c = .84;
    e_f = .82;
    pi_b = .92;
    eta_b = .94;
    e_t = .83;
    pi_n = .97;
    To4max = 1390; %[K]
elseif year > 1985 && year <= 2005        
    pi_dmax = .98;
    e_c = .88;
    e_f = .86;
    pi_b = .94;
    eta_b = .99;
    e_t = .87;
    pi_n = .98;
    To4max = 1780; %[K] 
elseif year > 2005 && year <= 2020    
    pi_dmax = .995;
    e_c = .90;
    e_f = .89;
    pi_b = .96;
    eta_b = .995;
    e_t = .89;
    pi_n = .995;
    To4max = 2000; %[K]   
else 
    error('Enter a valid year between 1945 and 2020')
end
end