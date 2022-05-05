function [sol_A, sol_B] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau,sigma,a_bar,lambda1,omega1,type)

    % change lambda and omega
    lambda = lambda1;
    omega = omega1;
    
    % first solve after the shock:
    % boundary values after the shock:
    p_BA = (mu_B-tau)/r;
    p_GA = (mu_G+omega)/r;
    
    syms x
    expr = matlabFunction(...
        taylor( 2/( sigma^2* x^2*(1-x)^2 ), x, 'ExpansionPoint',.5,'Order',5));
    
    ode_fun = @(z,y) odeAShock(z,y,r,mu_G,mu_B,tau,omega,expr,a_bar,type);
    bc_fun = @(ya, yb) bc(ya, yb, p_BA, p_GA);
    
    % obtain init
    xmesh = linspace(0,1,5000);
    y0_guess = [p_BA;0]; % start is known, V' is unknown. If code gives error, vary the derivative.
    guessFun = @(z) guessAShock(z,y0_guess,r,mu_G,mu_B,tau,omega,expr,a_bar,type);
    solinit = bvpinit(xmesh, guessFun);
    
    bvpoptions = bvpset(Stats="on",Nmax=150000,AbsTol=1e-4,RelTol=1e-4);
    sol_AShock = bvp5c(ode_fun,bc_fun,solinit,bvpoptions);
    disp('From first best after shock')
    
    % Fit value function after shock
    pol = 17;
    fit_AShock = polyfit(sol_AShock.x, sol_AShock.y(1,:),pol);
    
    syms zi
    AShock_fun = fit_AShock*zi.^(17:-1:0)';
    AShock_fun = matlabFunction(AShock_fun);
    
    % solve the model before the shock, using the polyfit after the shock
    p_BB = (mu_B + lambda*p_BA) / (lambda + r);
    p_GB = (mu_G + omega + lambda*p_GA)/(r+lambda);
    
    ode_fun = @(z,y) odeBShock(z,y,r,mu_B,mu_G,omega,lambda,expr,a_bar,AShock_fun,type);
    bc_fun = @(ya, yb) bc(ya, yb, p_BB, p_GB);
    
    % obtain init
    xmesh = linspace(0,1,5000);
    y0_guess = [p_BB;0]; % start is known, V' is unknown. If code gives error, vary the derivative.
    guessFun = @(z) guessBShock(z,y0_guess,r,mu_B,mu_G,omega,lambda,expr,a_bar,AShock_fun,type);
    solinit = bvpinit(xmesh, guessFun);
    
    bvpoptions = bvpset(Stats="on",Nmax=150000,AbsTol=1e-4,RelTol=1e-4);
    sol_BShock = bvp5c(ode_fun,bc_fun,solinit,bvpoptions);
    disp('From first best before shock')

    sol_A = sol_AShock;
    sol_B = sol_BShock;
end

