function dy = odeAgency(z,y,r,mu_B,mu_G,gamma,sigma,expr,exprTheta,exprThetaMinus,a_bar,option)
%ODE_V_FB Function to solve for the ODE of the first best value of V
%   write the second order ODE as a system of first-order equations. Define
%   y(1) := j(z) and y(2) := j'(z)

% option is a string. Choose: 'theta', 'scaled', 'unscaled'
% if option=='theta': g(z,a) = a^2*z^theta*(1-z)^theta/2
% if option=='scaled': g(z,a) = a^2*z*(1-z)/2
% if option=='unscaled': g(z,a) = a^2/2

if isequal(option,'scaled')
    % First, define the optimal effort level (a^*), depending on choice
    a = ( y(2)*z*(1-z) ) / ( 1 + gamma*r*sigma^2*z*(1-z) );
    if abs(a)>a_bar
        a=sign(a)*a_bar;
    end
    
    % ode stuff:
    dy = zeros(2,1);  
    dy(1) = y(2);
    
    dy(2) = ( expr(z) ) *  (r*y(1) + ( a^2*z*(1-z) )/2 ...
            + ( gamma*r )/2*a^2*z^2*(1-z)^2  - mu_G*z - mu_B*(1-z) ...
            - y(2)*a*z*(1-z) );
elseif isequal(option,'unscaled')
    % First, define the optimal effort level (a^*), depending on choice
    a = ( y(2)*z*(1-z) ) / ( 1 + gamma*r*sigma^2 );
    if abs(a)>a_bar
        a=sign(a)*a_bar;
    end

    % ode stuff:
    dy = zeros(2,1);  
    dy(1) = y(2);
    
    dy(2) = ( expr(z) ) *  (r*y(1) + ( a^2 )/2 ...
            + ( gamma*r )/2*a^2 - mu_G*z - mu_B*(1-z) ...
            - y(2)*a*z*(1-z) );
else
    % First, define the optimal effort level (a^*), depending on choice
    a = ( y(2)* exprThetaMinus(z) ) / (1 + gamma*r*sigma^2* (exprTheta(z)));
    if abs(a)>a_bar
        a=sign(a)*a_bar;
    end

    % ode stuff:
    dy = zeros(2,1);  
    dy(1) = y(2);
    
    dy(2) = ( expr(z) ) *  (r*y(1) + (a^2*exprTheta(z))/2 ...
            + (gamma*r)/2 *a^2 *exprTheta(z)^2*sigma^2  - mu_G*z - mu_B*(1-z) ...
            - y(2)*a*z*(1-z) );
end
end

