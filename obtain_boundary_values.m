function [p_B,p_G] = obtain_boundary_values(r, A_G, A_B ,delta,theta)
%OBTAIN_BOUNDARY_VALUES Summary of this function goes here
%   Detailed explanation goes here

i_B = A_B + sqrt( A_B^2 - 2*(A_B-delta-r)/theta );
i_G = A_G + sqrt( A_G^2 - 2*(A_G-delta-r)/theta );

% Check whether i_n >= 0:
if i_B < 0
    disp(i_B)
end
if i_G < 0
    disp(i_B)
end

p_B = (A_B - i_B) / (1 - i_B - theta/2*i_B^2 - delta); 
p_G = (A_G - i_G) / (1 - i_G - theta/2*i_G^2 - delta); 

