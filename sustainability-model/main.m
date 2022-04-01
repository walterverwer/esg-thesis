clear all, close all, clc
%% Code options
% Explanation of the code:
% This code solves an ODE regarding the degree of sustainability (z in
% [0,1], z=1 means 100% green and z=0 means 100% brown). It does so for a 
% model without agency frictions (i.e. first best)
% and a model with agency frictions. Agency frictions require a cost of
% effort function of the agent. The variable gFunOption represents the type
% of functions implemented. The most general is
% g(z,a)=a^2*(z*(1-z))^theta/2. In the options this can be enabled via
% 'theta'. If one chooses 'scaled', then theta=1, and if one chooses
% 'unscaled', then theta=0. These two cases represent numerically
% problematic cases, hence they have to be manually selected. The ODE is
% solved via MATLAB's boundary value problem solver (bvp). How this exactly
% works can be found in their documentation. My guess should be a
% relatively robust shooting algorithm. It is based on ode89(.) of MATLAB.

% Load in the parameters:
parameters;

% Boundary values, denoted by p_B and p_G:
p_B = mu_B/r;
p_G = mu_G/r;

% choose from: 'theta': g(theta,z,a), 'scaled': g(theta=1,z,a), 'unscaled': g(theta=0,z,a)
gFunOption = 'scaled'; % I recommend either scaled or unscaled.
%% Solve the boundary value problem with agency frictions

% Taylor expansion for singularity terms:
syms x
expr = matlabFunction(...
    taylor( 2/( sigma^2* x^2*(1-x)^2 ), x, 'ExpansionPoint',.5,'Order',5));

exprTheta = matlabFunction(...
    taylor( (x*(1-x))^theta, x, 'ExpansionPoint',.5,'Order',7));

exprThetaMinus = matlabFunction(...
    taylor( (x*(1-x))^(1-theta), x, 'ExpansionPoint',.5,'Order',7));


ode_fun = @(z,y) odeAgency(z,y,r,mu_B,mu_G,gamma,sigma,expr,exprTheta,exprThetaMinus,a_bar,gFunOption);
bc_fun = @(ya, yb) bc(ya, yb, p_B, p_G);

% obtain init
xmesh = linspace(0,1,1000);
y0_guess = [p_B;0.1]; % start is known, V' is unknown. If code gives error, vary the derivative.
guessFun = @(z) guessAgency(z,y0_guess,r,mu_G,mu_B,gamma,sigma,expr,exprTheta,exprThetaMinus,a_bar,gFunOption);
solinit = bvpinit(xmesh, guessFun);

bvpoptions = bvpset(Stats="on",Nmax=100000,AbsTol=1e-4,RelTol=1e-4);
sol_agency = bvp5c(ode_fun,bc_fun,solinit,bvpoptions);

% Plot results (same basic plots, not necessarily fancy)
plotAgency(sol_agency.x,sol_agency.y(1,:),mu_G,mu_B,gFunOption)



%% Solve the boundary value problem first best case

% Taylor expansion for singularity terms:
syms x
expr = matlabFunction(...
    taylor( 2/( sigma^2* x^2*(1-x)^2 ), x, 'ExpansionPoint',.5,'Order',5));

ode_fun = @(z,y) odeFb(z,y,r,mu_B,mu_G,expr,exprTheta,exprThetaMinus,a_bar,gFunOption);
bc_fun = @(ya, yb) bc(ya, yb, p_B, p_G);

% obtain init
xmesh = linspace(0,1,1000);
y0_guess = [p_B;0]; % start is known, V' is unknown. If code gives error, vary the derivative.
guessFun = @(z) guessFb(z,y0_guess,r,mu_G,mu_B,expr,exprTheta,exprThetaMinus,a_bar,gFunOption);
solinit = bvpinit(xmesh, guessFun);

bvpoptions = bvpset(Stats="on",Nmax=100000,AbsTol=1e-4,RelTol=1e-4);
sol_fb = bvp5c(ode_fun,bc_fun,solinit,bvpoptions);

% Plot results (same basic plots, not necessarily fancy)
plotFb(sol_fb.x,sol_fb.y(1,:),mu_G,mu_B,gFunOption)


%% Stranded Asset Risk and Green Preferences First Best Case
parameters;

% first solve after the shock:
% boundary values after the shock:
p_BA = 0;
p_GA = (mu_G+omega)/r;

syms x
expr = matlabFunction(...
    taylor( 2/( sigma^2* x^2*(1-x)^2 ), x, 'ExpansionPoint',.5,'Order',5));

ode_fun = @(z,y) odeAShock(z,y,r,mu_G,omega,expr,a_bar);
bc_fun = @(ya, yb) bc(ya, yb, p_BA, p_GA);

% obtain init
xmesh = linspace(0,1,2000);
y0_guess = [p_BA;0]; % start is known, V' is unknown. If code gives error, vary the derivative.
guessFun = @(z) guessAShock(z,y0_guess,r,mu_G,omega,expr,a_bar);
solinit = bvpinit(xmesh, guessFun);

bvpoptions = bvpset(Stats="on",Nmax=100000,AbsTol=1e-4,RelTol=1e-4);
sol_AShock = bvp5c(ode_fun,bc_fun,solinit,bvpoptions);

% Plot results (same basic plots, not necessarily fancy)
plot(sol_AShock.x,sol_AShock.y(1,:))
grid on
title('After Shock, $g(z,a)=\frac{a^2z(1-z)}{2}$', 'interpreter','latex')
saveas(gca,[pwd '/figures/various_gFun_plots/after_shock/sol_after_shock.eps'],'epsc')
%saveas(gca,[pwd '/figures/various_gFun_plots/after_shock/sol_after_shock.jpeg'],'jpeg')

% Fit value function after shock
pol = 17;
fit_AShock = polyfit(sol_AShock.x, sol_AShock.y(1,:),pol);

syms zi
AShock_fun = fit_AShock*zi.^(17:-1:0)';
AShock_fun = matlabFunction(AShock_fun);

% solve the model before the shock, using the polyfit after the shock
p_BB = mu_B / (lambda + r);
p_GB = (mu_G + omega + lambda*p_GA)/(r+lambda);

ode_fun = @(z,y) odeBShock(z,y,r,mu_B,mu_G,omega,lambda,expr,a_bar,AShock_fun);
bc_fun = @(ya, yb) bc(ya, yb, p_BB, p_GB);

% obtain init
xmesh = linspace(0,1,2000);
y0_guess = [p_BB;0]; % start is known, V' is unknown. If code gives error, vary the derivative.
guessFun = @(z) guessBShock(z,y0_guess,r,mu_B,mu_G,omega,lambda,expr,a_bar,AShock_fun);
solinit = bvpinit(xmesh, guessFun);

bvpoptions = bvpset(Stats="on",Nmax=100000,AbsTol=1e-4,RelTol=1e-4);
sol_BShock = bvp5c(ode_fun,bc_fun,solinit,bvpoptions);

% Plot results (same basic plots, not necessarily fancy)
plot(sol_BShock.x,sol_BShock.y(1,:))
grid on
title('Before Shock, $g(z,a)=\frac{a^2z(1-z)}{2}$', 'interpreter','latex')
saveas(gca,[pwd '/figures/various_gFun_plots/after_shock/sol_before_shock.eps'],'epsc')
%saveas(gca,[pwd '/figures/various_gFun_plots/after_shock/sol_after_shock.jpeg'],'jpeg')



%% Comparison plots

if isequal(gFunOption,'scaled')
    plot(sol_fb.x,sol_fb.y(1,:))
    hold on
    plot(sol_agency.x,sol_agency.y(1,:))
    grid on
    legend('First Best', 'With Agency Friction')
    title('Comparison first best with agency friction, $\mu^G \neq \mu^B$ for $g(z,a)=\frac{a^2z(1-z)}{2}$', 'interpreter','latex')
    saveas(gca,[pwd '/figures/various_gFun_plots/comparison/sol_unequal_comp_scaled.eps'],'epsc')

elseif isequal(gFunOption,'unscaled')
    plot(sol_fb.x,sol_fb.y(1,:))
    hold on
    plot(sol_agency.x,sol_agency.y(1,:))
    grid on
    legend('First Best', 'With Agency Friction')
    title('Comparison first best with agency friction, $\mu^G \neq \mu^B$ for $g(z,a)=\frac{a^2}{2}$', 'interpreter','latex')
    saveas(gca,[pwd '/figures/various_gFun_plots/comparison/sol_unequal_comp_unscaled.eps'],'epsc')

else
    plot(sol_fb.x,sol_fb.y(1,:))
    hold on
    plot(sol_agency.x,sol_agency.y(1,:))
    grid on
    legend('First Best', 'With Agency Friction')
    title('Comparison first best with agency friction, $\mu^G \neq \mu^B$ for $g(z,a)=\frac{a^2z^\theta(1-z)^\theta}{2}$', 'interpreter','latex')
    saveas(gca,[pwd '/figures/various_gFun_plots/comparison/sol_unequal_comp_theta.eps'],'epsc')
end



%% FB SAR GP - a^* -> d lambda & d omega, before shock

% for fixed omega
lambda_low = 0.01;
[~, sol_lambda_low] = solve_ODE_SAR_FB(r,mu_B,mu_G,sigma,a_bar,lambda_low,omega);

lambda_med = 0.02;
[~, sol_lambda_med] = solve_ODE_SAR_FB(r,mu_B,mu_G,sigma,a_bar,lambda_med,omega);

lambda_high = 0.05;
[~, sol_lambda_high] = solve_ODE_SAR_FB(r,mu_B,mu_G,sigma,a_bar,lambda_high,omega);

figure
plot(sol_lambda_low.x, sol_lambda_low.y(2,:))
hold on
grid on
plot(sol_lambda_med.x, sol_lambda_med.y(2,:))
plot(sol_lambda_high.x, sol_lambda_high.y(2,:))
legend('Low lambda', 'Medium lambda', 'High lambda')


% for fixed lambda
omega_zero = 0.0;
[~, sol_omega_zero] = solve_ODE_SAR_FB(r,mu_B,mu_G,sigma,a_bar,lambda,omega_zero);

omega_low = 0.01;
[~, sol_omega_low] = solve_ODE_SAR_FB(r,mu_B,mu_G,sigma,a_bar,lambda,omega_low);

omega_med = 0.02;
[~, sol_omega_med] = solve_ODE_SAR_FB(r,mu_B,mu_G,sigma,a_bar,lambda,omega_med);

omega_high = 0.05;
[~, sol_omega_high] = solve_ODE_SAR_FB(r,mu_B,mu_G,sigma,a_bar,lambda,omega_high);

figure
plot(sol_omega_zero.x, sol_omega_zero.y(2,:))
hold on
grid on
plot(sol_omega_low.x, sol_omega_low.y(2,:))
plot(sol_omega_med.x, sol_omega_med.y(2,:))
plot(sol_omega_high.x, sol_omega_high.y(2,:))
legend('Zero Omega', 'Low Omega', 'Medium Omega', 'High Omega')

