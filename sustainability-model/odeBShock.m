function dy = odeBShock(z,y,r,mu_B,mu_G,xi,lambda,expr,a_bar,AShock_fun)
%   ode function to solve for j_s(z)
%   write the second order ODE as a system of first-order equations. Define
%   y(1) := j(z) and y(2) := j'(z)

a = y(2);
if abs(a)>a_bar
    a=sign(a)*a_bar;
end

% ode after a shock
dy = zeros(2,1);  
dy(1) = y(2);

dy(2) = ( expr(z) ) * ((r+lambda)*y(1) + ( a^2*z*(1-z) )/2 ...
            - mu_G*z - mu_B*(1-z) - xi*z - lambda*AShock_fun(z) - y(2)*a*z*(1-z) );

end

