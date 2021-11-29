function [p_B,p_G,i_B,i_G] = obtain_boundary_values(r, A_G, A_B ,delta,theta)
%OBTAIN_BOUNDARY_VALUES Summary of this function goes here
%   Detailed explanation goes here

% Check:
if  A_B^2 + 2*(delta + r - A_B)/theta < 0
    error('Investment not defined!')
    return
end

if  A_G^2 + 2*(delta + r - A_G)/theta < 0
    error('Investment defined!')
    return
end


% Should be the positive roots...
i_B1 = A_B + sqrt( A_B^2 + 2*(delta + r - A_B)/theta );
i_G1 = A_G + sqrt( A_G^2 + 2*(delta + r - A_G)/theta );

% Should be the negative roots...
i_B2 = A_B - sqrt( A_B^2 + 2*(delta+r-A_B)/theta );
i_G2 = A_G - sqrt( A_G^2 + 2*(delta+r-A_G)/theta );

% if i_B1>=0
%     i_B=i_B1;
% else
%     i_B=i_B2;
%     disp('using 2')
% end
% 
% if i_G1>=0
%     i_G=i_G1;
% else
%     i_G=i_G2;
%     disp('using 2')
% end

i_B=i_B1;
i_G=i_G1;

p_B = (A_B - i_B) / (1 - i_B - theta/2*i_B^2 - delta); 
p_G = (A_G - i_G) / (1 - i_G - theta/2*i_G^2 - delta); 

