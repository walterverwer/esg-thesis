function dy = ode_v_fb(z,y,r,sigma_G, sigma_B ,A_G, A_B ,delta,theta)
%ODE_V_FB Function to solve for the ODE of the first best value of V
%   write the second order ODE as a system of first-order equations. Define
%   y(1) := V(z) and y(2) := V'(z)

i_G = 1/theta*(1- 1/(y(1)-z*y(2))); % optimal i_G
i_B = 1/theta*(1-1/(y(1)+(1-z)*y(2))); % optimal i_B

dy = zeros(2,1);  
dy(1) = y(2);
dy(2) = 2/(sigma_B^2+sigma_G^2*z^2*(1-z)^2) * ( r*y(1) - (1-z)*(A_B-i_B) ...
    - z*(A_G-i_G) - (1-z)*(y(1) - z*y(2))*(i_B - theta/2*i_B^2 - delta ) ...
    - z*(y(1) + (1-z)*y(2))*(i_G - theta/2*i_G^2 - delta) )  ;

end

