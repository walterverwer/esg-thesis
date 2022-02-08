clear all, close all

%% Solve the boundary value problem
% Load in the parameters:
parameters;

% Boundary values, denoted by p_B and p_G:
p_B = mu_B/r;
p_G = mu_G/r;

% Taylor expansion for singularity terms:
syms x
expr = matlabFunction(...
    taylor( 2/( sigma^2* x^2*(1-x)^2 ), x, 'ExpansionPoint',.5,'Order',5));

exprTheta = matlabFunction(...
    taylor( (x*(1-x))^theta, x, 'ExpansionPoint',.5,'Order',7));

exprThetaMinus = matlabFunction(...
    taylor( (x*(1-x))^(1-theta), x, 'ExpansionPoint',.5,'Order',7));


ode_fun = @(z,y) ode(z,y,r,mu_B,mu_G,theta,gamma,sigma,expr, exprTheta, exprThetaMinus);
bc_fun = @(ya, yb) bc(ya, yb, p_B, p_G);

% obtain init
xmesh = linspace(0,1,1000);
y0_guess = [p_B;0.1]; % start is known, V' is unknown. Guess whole line
guessFun = @(z) guess(z,y0_guess,r ,mu_G, mu_B,theta,gamma,sigma,expr, exprTheta, exprThetaMinus);
solinit = bvpinit(xmesh, guessFun);

bvpoptions = bvpset(Stats="on",Nmax=10000,AbsTol=1e-4,RelTol=1e-4);
sol = bvp5c(ode_fun,bc_fun,solinit,bvpoptions);

% Plot results
if mu_G==mu_B
    plot(sol.x,sol.y(1,:))
    grid on
    title('$\mu^G = \mu^B$','interpreter','latex')
    saveas(gca,'sol_equal_fb.eps','epsc')
else
    plot(sol.x,sol.y(1,:))
    grid on
    title('$\mu^G \neq \mu^B$','interpreter','latex')
    saveas(gca,'sol_unequal_fb.eps','epsc')
end









