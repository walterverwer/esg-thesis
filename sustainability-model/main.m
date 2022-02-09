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
gFunOption = 'theta'; % I recommend either scaled or unscaled.
%% Solve the boundary value problem with agency frictions

% Taylor expansion for singularity terms:
syms x
expr = matlabFunction(...
    taylor( 2/( sigma^2* x^2*(1-x)^2 ), x, 'ExpansionPoint',.5,'Order',5));

exprTheta = matlabFunction(...
    taylor( (x*(1-x))^theta, x, 'ExpansionPoint',.5,'Order',7));

exprThetaMinus = matlabFunction(...
    taylor( (x*(1-x))^(1-theta), x, 'ExpansionPoint',.5,'Order',7));


ode_fun = @(z,y) odeAgency(z,y,r,mu_B,mu_G,gamma,sigma,expr,exprTheta,exprThetaMinus,gFunOption);
bc_fun = @(ya, yb) bc(ya, yb, p_B, p_G);

% obtain init
xmesh = linspace(0,1,1000);
y0_guess = [p_B;0.1]; % start is known, V' is unknown. If code gives error, vary the derivative.
guessFun = @(z) guessAgency(z,y0_guess,r,mu_G,mu_B,gamma,sigma,expr,exprTheta,exprThetaMinus,gFunOption);
solinit = bvpinit(xmesh, guessFun);

bvpoptions = bvpset(Stats="on",Nmax=100000,AbsTol=1e-4,RelTol=1e-4);
sol = bvp5c(ode_fun,bc_fun,solinit,bvpoptions);

% Plot results (same basic plots, not necessarily fancy)
plotAgency(sol.x,sol.y(1,:),mu_G,mu_B,gFunOption)



%% Solve the boundary value problem first best case

% Taylor expansion for singularity terms:
syms x
expr = matlabFunction(...
    taylor( 2/( sigma^2* x^2*(1-x)^2 ), x, 'ExpansionPoint',.5,'Order',5));

ode_fun = @(z,y) odeFb(z,y,r,mu_B,mu_G,expr,exprTheta,exprThetaMinus,gFunOption);
bc_fun = @(ya, yb) bc(ya, yb, p_B, p_G);

% obtain init
xmesh = linspace(0,1,1000);
y0_guess = [p_B;0]; % start is known, V' is unknown. If code gives error, vary the derivative.
guessFun = @(z) guessFb(z,y0_guess,r,mu_G,mu_B,expr,exprTheta,exprThetaMinus,gFunOption);
solinit = bvpinit(xmesh, guessFun);

bvpoptions = bvpset(Stats="on",Nmax=100000,AbsTol=1e-4,RelTol=1e-4);
sol = bvp5c(ode_fun,bc_fun,solinit,bvpoptions);

% Plot results (same basic plots, not necessarily fancy)
plotFb(sol.x,sol.y(1,:),mu_G,mu_B,gFunOption)
