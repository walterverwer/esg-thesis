function dy = ode(z,y,r,mu_B,mu_G,theta,gamma,sigma,expr, exprTheta, exprThetaMinus)
%ODE_V_FB Function to solve for the ODE of the first best value of V
%   write the second order ODE as a system of first-order equations. Define
%   y(1) := j(z) and y(2) := j'(z)

% First, define the optimal effort level (a^*):
a = ( y(2)* exprThetaMinus(z) ) / (1 + gamma*r*sigma^2* (exprTheta(z)));

dy = zeros(2,1);  
dy(1) = y(2);

dy(2) = ( expr(z) ) *  (r*y(1) + (a^2*exprTheta(z))/2 ...
        + (gamma*r)/2 *a^2 *exprTheta(z)^2*sigma^2  - mu_G*z - mu_B*(1-z) ...
        - y(2)*a*z*(1-z) );
end

