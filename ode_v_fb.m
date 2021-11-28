function dx = ode_v_fb(z,x,r,sigma_G, sigma_B ,A_G, A_B ,delta,theta)
%ODE_V_FB Function to solve for the ODE of the first best value of V
%   write the second order ODE as a system of first-order equations. Define
%   x(1) := V(z) and x(2) := V'(z)

i_G = 1/theta*(x(1)-z*x(2)-1); % optimal i_G
i_B = 1/theta*(x(1)+(1-z)*x(2)-1); % optimal i_B

dx = zeros(2,1);  
dx(1) = x(2);
dx(2) = 2/(sigma_B^2+sigma_G^2*z^2*(1-z)^2) * ( r*x(1) - (1-z)*(A_B-i_B) ...
    - z*(A_G-i_G) - (1-z)*(x(1) - z*x(2))*(i_B + theta/2*i_B^2 - delta ) ...
    - z*(x(1) + (1-z)*x(2))*(i_G + theta/2*i_G^2 - delta) )  ;

end

