function dy = odeBShockAgency(z,y,r,mu_B,mu_G,omega,lambda,gamma,sigma,expr,expr_a,a_bar,AShock_fun,type)
%   ode function to solve for j_s(z)
%   write the second order ODE as a system of first-order equations. Define
%   y(1) := j(z) and y(2) := j'(z)



if isequal(type,'scaled')
    a = expr_a(z)*y(2);
    if abs(a)>a_bar
        a=sign(a)*a_bar;
    end
    
    % ode after a shock
    dy = zeros(2,1);  
    dy(1) = y(2);

    dy(2) = ( expr(z) ) * ((r+lambda)*y(1) + ( a^2*z*(1-z) )/2 ...
            - mu_G*z - mu_B*(1-z) - omega*z - lambda*AShock_fun(z) - y(2)*a*z*(1-z) ...
            + (gamma*r)/2 * (a^2*z^2*(1-z)^2*sigma^2));

elseif isequal(type,'unscaled')
    a = y(2)*z*(1-z);
    %if abs(a)>a_bar
        %a=sign(a)*a_bar;
    %end
    
    % ode after a shock
    dy = zeros(2,1);  
    dy(1) = y(2);

    dy(2) = ( expr(z) ) * ((r+lambda)*y(1) + ( a^2 )/2 ...
            - mu_G*z - mu_B*(1-z) - omega*z - lambda*AShock_fun(z) - y(2)*a*z*(1-z) );
end

end

