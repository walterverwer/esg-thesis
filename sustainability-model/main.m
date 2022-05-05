clear all, close all, clc

%% Code options
% Explanation of the code:
% This code solves an ODE regarding the degree of sustainability (z in
% [0,1], z=1 means 100% green and z=0 means 100% brown).The ODE is
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
type = 'scaled'; %% ONLY WORKS FOR SCALED YET!!!

% first solve after the shock:
% boundary values after the shock:
p_BA = 0;
p_GA = (mu_G+omega)/r;

syms x
expr = matlabFunction(...
    taylor( 2/( sigma^2* x^2*(1-x)^2 ), x, 'ExpansionPoint',.5,'Order',5));

ode_fun = @(z,y) odeAShock(z,y,r,mu_G,omega,expr,a_bar,type);
bc_fun = @(ya, yb) bc(ya, yb, p_BA, p_GA);

% obtain init
xmesh = linspace(0,1,2000);
y0_guess = [p_BA;0]; % start is known, V' is unknown. If code gives error, vary the derivative.
guessFun = @(z) guessAShock(z,y0_guess,r,mu_G,omega,expr,a_bar,type);
solinit = bvpinit(xmesh, guessFun);

bvpoptions = bvpset(Stats="on",Nmax=100000,AbsTol=1e-4,RelTol=1e-4);
sol_AShock = bvp5c(ode_fun,bc_fun,solinit,bvpoptions);

% Plot results (same basic plots, not necessarily fancy)
figure;
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

ode_fun = @(z,y) odeBShock(z,y,r,mu_B,mu_G,omega,lambda,expr,a_bar,AShock_fun,type);
bc_fun = @(ya, yb) bc(ya, yb, p_BB, p_GB);

% obtain init
xmesh = linspace(0,1,2000);
y0_guess = [p_BB;0]; % start is known, V' is unknown. If code gives error, vary the derivative.
guessFun = @(z) guessBShock(z,y0_guess,r,mu_B,mu_G,omega,lambda,expr,a_bar,AShock_fun,type);
solinit = bvpinit(xmesh, guessFun);

bvpoptions = bvpset(Stats="on",Nmax=100000,AbsTol=1e-4,RelTol=1e-4);
sol_BShock = bvp5c(ode_fun,bc_fun,solinit,bvpoptions);

% Plot results (same basic plots, not necessarily fancy)
figure;
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
parameters;
type = 'scaled'; %% ONLY WORKS FOR SCALED YET!!!

% for fixed omega
lambda_low = 0.01;
[~, sol_lambda_low] = solve_ODE_SAR_FB(r,mu_B,mu_G,sigma,a_bar,lambda_low,omega,type);

lambda_med = 0.02;
[~, sol_lambda_med] = solve_ODE_SAR_FB(r,mu_B,mu_G,sigma,a_bar,lambda_med,omega,type);

lambda_high = 0.05;
[~, sol_lambda_high] = solve_ODE_SAR_FB(r,mu_B,mu_G,sigma,a_bar,lambda_high,omega,type);

y_lambda_low = obtain_a(sol_lambda_low.y(2,:), a_bar);
y_lambda_med = obtain_a(sol_lambda_med.y(2,:), a_bar);
y_lambda_high = obtain_a(sol_lambda_high.y(2,:), a_bar);

figure
plot(sol_lambda_low.x, y_lambda_low)
hold on
grid on
plot(sol_lambda_med.x, y_lambda_med)
plot(sol_lambda_high.x, y_lambda_high)
legend('Low lambda', 'Medium lambda', 'High lambda')
title('Comparison for Different Values of $\lambda$ on $a^*_{ns}$', 'interpreter','latex')
saveas(gca,[pwd '/figures/fb_cross_derivatives/comparison_lambda.eps'],'epsc')

% for fixed lambda
omega_zero = 0.0;
[~, sol_omega_zero] = solve_ODE_SAR_FB(r,mu_B,mu_G,sigma,a_bar,lambda,omega_zero,type);

omega_low = 0.01;
[~, sol_omega_low] = solve_ODE_SAR_FB(r,mu_B,mu_G,sigma,a_bar,lambda,omega_low,type);

omega_med = 0.02;
[~, sol_omega_med] = solve_ODE_SAR_FB(r,mu_B,mu_G,sigma,a_bar,lambda,omega_med,type);

omega_high = 0.05;
[~, sol_omega_high] = solve_ODE_SAR_FB(r,mu_B,mu_G,sigma,a_bar,lambda,omega_high,type);


y_omega_zero = obtain_a(sol_omega_zero.y(2,:), a_bar);
y_omega_low = obtain_a(sol_omega_low.y(2,:), a_bar);
y_omega_med = obtain_a(sol_omega_med.y(2,:), a_bar);
y_omega_high = obtain_a(sol_omega_high.y(2,:), a_bar);

figure
plot(sol_omega_zero.x, y_omega_zero)
hold on
grid on
plot(sol_omega_low.x, y_omega_low)
plot(sol_omega_med.x, y_omega_med)
plot(sol_omega_high.x, y_omega_high)
legend('Zero omega', 'Low omega', 'Medium omega', 'High omega')
title('Comparison for Different Values of $\omega$ on $a^*_{ns}$', 'interpreter','latex')
saveas(gca,[pwd '/figures/fb_cross_derivatives/comparison_omega.eps'],'epsc')


%% Plots from analytical derivatives (wrong)
parameters; 
% From the model we know a^*_i = j'_i(z) or a^*_i = j'_i(z)z(1-z)

type = 'scaled'; % scaled: g(z,a), unscaled: g(a)
[sol_AShock, sol_BShock] = solve_ODE_SAR_FB(r,mu_B,mu_G,sigma,a_bar,lambda,omega,type);


% Fit value function after shock
pol = 25;
fit_AShock = polyfit(sol_AShock.x, sol_AShock.y(1,:),pol);

syms zi
AShock_fun = fit_AShock*zi.^(pol:-1:0)';
AShock_fun = matlabFunction(AShock_fun);

% Fit value function before shock
fit_BShock = polyfit(sol_BShock.x, sol_BShock.y(1,:),pol);

syms zi
BShock_fun = fit_BShock*zi.^(pol:-1:0)';
BShock_fun = matlabFunction(BShock_fun);

% Fit derivative of the value function before shock
fit_dBShock = polyfit(sol_BShock.x, sol_BShock.y(2,:),pol);

syms zi
dBShock_fun = fit_dBShock*zi.^(pol:-1:0)';
dBShock_fun = matlabFunction(dBShock_fun);

if isequal(type,'unscaled')
    % Plot da_b / dlambda for g(a)
    figure('name', 'Effect Lambda');
    fplot(@(z)  (BShock_fun(z) - AShock_fun(z)) / ( dBShock_fun(z)*z*(1-z) ), [0.05 0.95] )
    title('Effect of Stranded Asset Risk on effort, before the shock: $\frac{\partial a^*_{ns}}{\partial \lambda}$', 'interpreter','latex')
    saveas(gca,[pwd '/figures/fb_cross_derivatives/dadlambda_before_shock.eps'],'epsc')
    
    
    % Plot da_b / domega for g(a)
    figure('name', 'Effect Omega');
    fplot(@(z) ( (1+r) / (dBShock_fun(z)*(z-1)*r) ), [0.05 0.95] )
    title('Effect of a Preference for Sustainability on effort, before the shock: $\frac{\partial a^*_{ns}}{\partial \omega}$', 'interpreter','latex')
    saveas(gca,[pwd '/figures/fb_cross_derivatives/dadomega_before_shock.eps'],'epsc')
    
elseif isequal(type,'scaled')  
    % Plot da_b / dlambda for g(z,a)
    figure('name', 'Effect Lambda');
    fplot(@(z)  (BShock_fun(z) - AShock_fun(z)) / ( dBShock_fun(z)*z*(1-z) ), [0.05 0.95] )
    title('Effect of Stranded Asset Risk on effort, before the shock: $\frac{\partial a^*_{ns}}{\partial \lambda}$', 'interpreter','latex')
    saveas(gca,[pwd '/figures/fb_cross_derivatives/dadlambda_before_shock.eps'],'epsc')
    
    
    % Plot da_b / domega for g(z,a)
    figure('name', 'Effect Omega');
    fplot(@(z) ((1+r) / (dBShock_fun(z)*(z-1)*r ) ), [0.05 0.95] )
    title('Effect of a Preference for Sustainability on effort, before the shock: $\frac{\partial a^*_{ns}}{\partial \omega}$', 'interpreter','latex')
    saveas(gca,[pwd '/figures/fb_cross_derivatives/dadomega_before_shock.eps'],'epsc')
end

%% check
figure('name', 'A shock')
fplot(@(z) AShock_fun(z), [0 1])

figure('name', 'B shock')
fplot(@(z) BShock_fun(z), [0 1])

figure('name', 'dB shock')
fplot(@(z) dBShock_fun(z), [0 1])


%% Solve analytical ODE fb
parameters;

[sol_AShock, sol_BShock] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,type);

figure;
plot(sol_BShock.x, sol_BShock.y(1,:))
title('before')

figure;
plot(sol_AShock.x, sol_AShock.y(1,:))
title('after')



%% Approximate (cross) derivative with finite difference
% Cross derivative tau omega

parameters;
% Apply (central) finite difference method to obtain an approximation of
% the cross derivative for omega and lambda on effort (a).

% Step 1: define step sizes
h_omega = omega*0.5; % change in omega
h_tau = tau*0.1; % change in tau

% Step 2: compute the four value functions for different parameters
type = 'scaled'; % scaled: g(z,a)

% Filling in the formula for finite difference (see notes)
[~, sol_BShock_1] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau+h_tau,sigma,a_bar,lambda,omega+h_omega,type);
[~, sol_BShock_2] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau+h_tau,sigma,a_bar,lambda,omega-h_omega,type);
[~, sol_BShock_3] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau-h_tau,sigma,a_bar,lambda,omega+h_omega,type);
[~, sol_BShock_4] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau-h_tau,sigma,a_bar,lambda,omega-h_omega,type);

% Step 3: fit the derivative to obtain a comparable 'function'
pol = 25;
fit_dBShock1 = polyfit(sol_BShock_1.x, sol_BShock_1.y(2,:),pol);
fit_dBShock2 = polyfit(sol_BShock_2.x, sol_BShock_2.y(2,:),pol);
fit_dBShock3 = polyfit(sol_BShock_3.x, sol_BShock_3.y(2,:),pol);
fit_dBShock4 = polyfit(sol_BShock_4.x, sol_BShock_4.y(2,:),pol);

% Step 4: obtain functional form:
syms zi
BShock_fun1 = fit_dBShock1*zi.^(pol:-1:0)';
BShock_fun1 = matlabFunction(BShock_fun1);

BShock_fun2 = fit_dBShock2*zi.^(pol:-1:0)';
BShock_fun2 = matlabFunction(BShock_fun2);

BShock_fun3 = fit_dBShock3*zi.^(pol:-1:0)';
BShock_fun3 = matlabFunction(BShock_fun3);

BShock_fun4 = fit_dBShock4*zi.^(pol:-1:0)';
BShock_fun4 = matlabFunction(BShock_fun4);

% Step 5: compute central finite difference
cross_tau_omega = @(z) (BShock_fun1(z) - BShock_fun2(z) - BShock_fun3(z) + BShock_fun4(z)) / 4*h_omega*h_tau;

figure('name', 'Numerical cross derivative')
fplot(cross_tau_omega, [0 1])

% Try with makima to get cross derivative:
xq = 0:0.0001:1;
m_ShockB_1 = makima(sol_BShock_1.x, sol_BShock_1.y(2,:), xq);
m_ShockB_2 = makima(sol_BShock_2.x, sol_BShock_2.y(2,:), xq);
m_ShockB_3 = makima(sol_BShock_3.x, sol_BShock_3.y(2,:), xq);
m_ShockB_4 = makima(sol_BShock_4.x, sol_BShock_4.y(2,:), xq);

cross_tau_omega_fb = (m_ShockB_1 - m_ShockB_2 - m_ShockB_3 + m_ShockB_4) / 4*h_tau*h_omega;

figure('name', '(MAKIMA) Numerical cross derivative tau omega fb')
plot(xq,cross_tau_omega_fb)
title('Numerical cross derivative tau omega (first best)')
saveas(gca,[pwd '/figures/agency_derivatives/makima_cross_tau_omega_fb.eps'],'epsc')

save('cross_derivative_tau_omega_fb.mat','cross_tau_omega_fb')


%% Cross derivative lambda omega

parameters;
% Apply (central) finite difference method to obtain an approximation of
% the cross derivative for omega and lambda on effort (a).

% Step 1: define step sizes
h_omega = omega*0.5; % change in omega
h_lambda = lambda*0.5; % change in lambda

% Step 2: compute the four value functions for different parameters
type = 'scaled'; % scaled: g(z,a)

% Filling in the formula for finite difference (see notes)
[~, sol_BShock_1] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau,sigma,a_bar,lambda+h_lambda,omega+h_omega,type);
[~, sol_BShock_2] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau,sigma,a_bar,lambda+h_lambda,omega-h_omega,type);
[~, sol_BShock_3] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau,sigma,a_bar,lambda-h_lambda,omega+h_omega,type);
[~, sol_BShock_4] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau,sigma,a_bar,lambda-h_lambda,omega-h_omega,type);

% Step 3: fit the derivative to obtain a comparable 'function'
pol = 25;
fit_dBShock1 = polyfit(sol_BShock_1.x, sol_BShock_1.y(2,:),pol);
fit_dBShock2 = polyfit(sol_BShock_2.x, sol_BShock_2.y(2,:),pol);
fit_dBShock3 = polyfit(sol_BShock_3.x, sol_BShock_3.y(2,:),pol);
fit_dBShock4 = polyfit(sol_BShock_4.x, sol_BShock_4.y(2,:),pol);

% Step 4: obtain functional form:
syms zi
BShock_fun1 = fit_dBShock1*zi.^(pol:-1:0)';
BShock_fun1 = matlabFunction(BShock_fun1);

BShock_fun2 = fit_dBShock2*zi.^(pol:-1:0)';
BShock_fun2 = matlabFunction(BShock_fun2);

BShock_fun3 = fit_dBShock3*zi.^(pol:-1:0)';
BShock_fun3 = matlabFunction(BShock_fun3);

BShock_fun4 = fit_dBShock4*zi.^(pol:-1:0)';
BShock_fun4 = matlabFunction(BShock_fun4);

% Step 5: compute central finite difference
cross_lambda_omega = @(z) (BShock_fun1(z) - BShock_fun2(z) - BShock_fun3(z) + BShock_fun4(z)) / 4*h_omega*h_lambda;

figure('name', 'Numerical cross derivative')
fplot(cross_lambda_omega, [0.1 0.9])

% Try with makima to get cross derivative:
xq = 0:0.0001:1;
m_ShockB_1 = makima(sol_BShock_1.x, sol_BShock_1.y(2,:), xq);
m_ShockB_2 = makima(sol_BShock_2.x, sol_BShock_2.y(2,:), xq);
m_ShockB_3 = makima(sol_BShock_3.x, sol_BShock_3.y(2,:), xq);
m_ShockB_4 = makima(sol_BShock_4.x, sol_BShock_4.y(2,:), xq);

cross_lambda_omega_fb = (m_ShockB_1 - m_ShockB_2 - m_ShockB_3 + m_ShockB_4) / 4*h_lambda*h_omega;

figure('name', '(MAKIMA) Numerical cross derivative lambda omega fb')
plot(xq,cross_lambda_omega_fb)
title('Numerical cross derivative (first best)')
saveas(gca,[pwd '/figures/agency_derivatives/makima_cross_fb.eps'],'epsc')

save('cross_derivative_fb.mat','cross_lambda_omega_fb')

% Derivative lambda
parameters;
% Step 1: define step sizes
h_lambda = lambda*0.5; % change in lambda

% Step 2: compute the four value functions for different parameters
type = 'scaled'; % scaled: g(z,a)

% Filling in the formula for finite difference (see notes)
[~, sol_BShock_1] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau,sigma,a_bar,lambda+2*h_lambda,omega,type);
[~, sol_BShock_2] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau,sigma,a_bar,lambda+h_lambda,omega,type);
[~, sol_BShock_3] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau,sigma,a_bar,lambda-h_lambda,omega,type);
[~, sol_BShock_4] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau,sigma,a_bar,lambda-2*h_lambda,omega,type);

% Step 3: fit the derivative to obtain a comparable 'function'
pol = 25;
fit_dBShock1 = polyfit(sol_BShock_1.x, sol_BShock_1.y(2,:),pol);
fit_dBShock2 = polyfit(sol_BShock_2.x, sol_BShock_2.y(2,:),pol);
fit_dBShock3 = polyfit(sol_BShock_3.x, sol_BShock_3.y(2,:),pol);
fit_dBShock4 = polyfit(sol_BShock_4.x, sol_BShock_4.y(2,:),pol);

% Step 4: obtain functional form:
syms zi
BShock_fun1 = fit_dBShock1*zi.^(pol:-1:0)';
BShock_fun1 = matlabFunction(BShock_fun1);

BShock_fun2 = fit_dBShock2*zi.^(pol:-1:0)';
BShock_fun2 = matlabFunction(BShock_fun2);

BShock_fun3 = fit_dBShock3*zi.^(pol:-1:0)';
BShock_fun3 = matlabFunction(BShock_fun3);

BShock_fun4 = fit_dBShock4*zi.^(pol:-1:0)';
BShock_fun4 = matlabFunction(BShock_fun4);

% Step 5: compute central finite difference
derivative_lambda = @(z) (-BShock_fun1(z) + 8*BShock_fun2(z) - 8*BShock_fun3(z) + BShock_fun4(z)) / 12*h_lambda;

figure('name', 'Numerical derivative lambda')
fplot(derivative_lambda, [0.1 0.9])

% Try with makima to get cross derivative:

xq = 0:0.0001:1;
m_ShockB_1 = makima(sol_BShock_1.x, sol_BShock_1.y(2,:), xq);
m_ShockB_2 = makima(sol_BShock_2.x, sol_BShock_2.y(2,:), xq);
m_ShockB_3 = makima(sol_BShock_3.x, sol_BShock_3.y(2,:), xq);
m_ShockB_4 = makima(sol_BShock_4.x, sol_BShock_4.y(2,:), xq);

derivative_lambda_fb = (-m_ShockB_1 + 8*m_ShockB_2 - 8*m_ShockB_3 + m_ShockB_4) / 12*h_lambda;

figure('name', '(MAKIMA) Numerical derivative lambda agency')
plot(xq,derivative_lambda_fb)
title('Numerical derivative lambda (first best)')
saveas(gca,[pwd '/figures/agency_derivatives/makima_lambda_fb.eps'],'epsc')

save('derivative_lambda_fb.mat','derivative_lambda_fb')

% Derivative omega
parameters;
% Step 1: define step sizes
h_omega = omega*0.5; % change in lambda

% Step 2: compute the four value functions for different parameters
type = 'scaled'; % scaled: g(z,a)

% Filling in the formula for finite difference (see notes)
[~, sol_BShock_1] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega+2*h_omega,type);
[~, sol_BShock_2] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega+h_omega,type);
[~, sol_BShock_3] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega-h_omega,type);
[~, sol_BShock_4] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega-2*h_omega,type);

% Step 3: fit the derivative to obtain a comparable 'function'
pol = 25;
fit_dBShock1 = polyfit(sol_BShock_1.x, sol_BShock_1.y(2,:),pol);
fit_dBShock2 = polyfit(sol_BShock_2.x, sol_BShock_2.y(2,:),pol);
fit_dBShock3 = polyfit(sol_BShock_3.x, sol_BShock_3.y(2,:),pol);
fit_dBShock4 = polyfit(sol_BShock_4.x, sol_BShock_4.y(2,:),pol);

% Step 4: obtain functional form:
syms zi
BShock_fun1 = fit_dBShock1*zi.^(pol:-1:0)';
BShock_fun1 = matlabFunction(BShock_fun1);

BShock_fun2 = fit_dBShock2*zi.^(pol:-1:0)';
BShock_fun2 = matlabFunction(BShock_fun2);

BShock_fun3 = fit_dBShock3*zi.^(pol:-1:0)';
BShock_fun3 = matlabFunction(BShock_fun3);

BShock_fun4 = fit_dBShock4*zi.^(pol:-1:0)';
BShock_fun4 = matlabFunction(BShock_fun4);

% Step 5: compute central finite difference
derivative_omega = @(z) (-BShock_fun1(z) + 8*BShock_fun2(z) - 8*BShock_fun3(z) + BShock_fun4(z)) / 12*h_omega;

figure('name', 'Numerical derivative omega')
fplot(derivative_omega, [0.1 0.9])

% Try with makima to get cross derivative:
xq = 0:0.0001:1;
m_ShockB_1 = makima(sol_BShock_1.x, sol_BShock_1.y(2,:), xq);
m_ShockB_2 = makima(sol_BShock_2.x, sol_BShock_2.y(2,:), xq);
m_ShockB_3 = makima(sol_BShock_3.x, sol_BShock_3.y(2,:), xq);
m_aShockB_4 = makima(sol_BShock_4.x, sol_BShock_4.y(2,:), xq);

derivative_omega_fb = (-m_ShockB_1 + 8*m_ShockB_2 - 8*m_ShockB_3 + m_ShockB_4) / 12*h_omega;

figure('name', '(MAKIMA) Numerical derivative omega first best')
plot(xq,derivative_omega_fb)
title('Numerical derivative omega (first best)')
saveas(gca,[pwd '/figures/agency_derivatives/makima_omega_fb.eps'],'epsc')

save('derivative_omega_fb.mat','derivative_omega_fb')

%% Solve Agency Friction Model With Stranded Asset Risk and Green Preferences
parameters;
type='scaled';
[sol_A_Agency, sol_B_Agency] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);

figure;
plot(sol_B_Agency.x,sol_B_Agency.y(1,:))
grid on
title('Before Shock Agency')

figure;
plot(sol_A_Agency.x,sol_A_Agency.y(1,:))
grid on
title('After Shock Agency')


%% Compute Numerical Derivatives Lambda, Omega, Cross for Agency Frictions

% Compute Numerical Derivative Lambda Agency Friction Model

parameters;
% Step 1: define step sizes
h_lambda = lambda*0.5; % change in lambda

% Step 2: compute the four value functions for different parameters
type = 'scaled'; % scaled: g(z,a)

% Filling in the formula for finite difference (see notes)
[~, sol_BAgency_1] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau,sigma,a_bar,lambda+2*h_lambda,omega,gamma,type);
[~, sol_BAgency_2] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau,sigma,a_bar,lambda+h_lambda,omega,gamma,type);
[~, sol_BAgency_3] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau,sigma,a_bar,lambda-h_lambda,omega,gamma,type);
[~, sol_BAgency_4] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau,sigma,a_bar,lambda-2*h_lambda,omega,gamma,type);

% Step 3: fit the derivative to obtain a comparable 'function'
pol = 17;
fit_dBAgency1 = polyfit(sol_BAgency_1.x, sol_BAgency_1.y(2,:),pol);
fit_dBAgency2 = polyfit(sol_BAgency_2.x, sol_BAgency_2.y(2,:),pol);
fit_dBAgency3 = polyfit(sol_BAgency_3.x, sol_BAgency_3.y(2,:),pol);
fit_dBAgency4 = polyfit(sol_BAgency_4.x, sol_BAgency_4.y(2,:),pol);

% Step 4: obtain functional form:
syms zi
BAgency_fun1 = fit_dBAgency1*zi.^(pol:-1:0)';
BAgency_fun1 = matlabFunction(BAgency_fun1);

BAgency_fun2 = fit_dBAgency2*zi.^(pol:-1:0)';
BAgency_fun2 = matlabFunction(BAgency_fun2);

BAgency_fun3 = fit_dBAgency3*zi.^(pol:-1:0)';
BAgency_fun3 = matlabFunction(BAgency_fun3);

BAgency_fun4 = fit_dBAgency4*zi.^(pol:-1:0)';
BAgency_fun4 = matlabFunction(BAgency_fun4);

% Step 5: compute central finite difference
derivative_lambda = @(z) (-BAgency_fun1(z) + 8*BAgency_fun2(z) - 8*BAgency_fun3(z) + BAgency_fun4(z)) / 12*h_lambda;

figure('name', 'Numerical derivative lambda agency')
fplot(derivative_lambda, [0 1])

% Try with makima to get cross derivative:
xq = 0:0.0001:1;
m_agencyB_1 = makima(sol_BAgency_1.x, sol_BAgency_1.y(2,:), xq);
m_agencyB_2 = makima(sol_BAgency_2.x, sol_BAgency_2.y(2,:), xq);
m_agencyB_3 = makima(sol_BAgency_3.x, sol_BAgency_3.y(2,:), xq);
m_agencyB_4 = makima(sol_BAgency_4.x, sol_BAgency_4.y(2,:), xq);

derivative_lambda_agency = (-m_agencyB_1 + 8*m_agencyB_2 - 8*m_agencyB_3 + m_agencyB_4) / 12*h_lambda;

figure('name', '(MAKIMA) Numerical derivative lambda agency')
plot(xq,derivative_lambda_agency)
title('Numerical derivative lambda (agency)')
saveas(gca,[pwd '/figures/agency_derivatives/makima_lambda_agency.eps'],'epsc')

save('derivative_lambda_agency.mat','derivative_lambda_agency')

% Compute Numerical Derivative Omega Agency Friction Model

parameters;
% Step 1: define step sizes
h_omega = omega*0.5; % change in omega

% Step 2: compute the four value functions for different parameters
type = 'scaled'; % scaled: g(z,a)

% Filling in the formula for finite difference (see notes)
[~, sol_BAgency_1] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega+2*h_omega,gamma,type);
[~, sol_BAgency_2] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega+h_omega,gamma,type);
[~, sol_BAgency_3] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega-h_omega,gamma,type);
[~, sol_BAgency_4] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega-2*h_omega,gamma,type);

% Step 3: fit the derivative to obtain a comparable 'function'
pol = 17;
fit_dBAgency1 = polyfit(sol_BAgency_1.x, sol_BAgency_1.y(2,:),pol);
fit_dBAgency2 = polyfit(sol_BAgency_2.x, sol_BAgency_2.y(2,:),pol);
fit_dBAgency3 = polyfit(sol_BAgency_3.x, sol_BAgency_3.y(2,:),pol);
fit_dBAgency4 = polyfit(sol_BAgency_4.x, sol_BAgency_4.y(2,:),pol);

% Step 4: obtain functional form:
syms zi
BAgency_fun1 = fit_dBAgency1*zi.^(pol:-1:0)';
BAgency_fun1 = matlabFunction(BAgency_fun1);

BAgency_fun2 = fit_dBAgency2*zi.^(pol:-1:0)';
BAgency_fun2 = matlabFunction(BAgency_fun2);

BAgency_fun3 = fit_dBAgency3*zi.^(pol:-1:0)';
BAgency_fun3 = matlabFunction(BAgency_fun3);

BAgency_fun4 = fit_dBAgency4*zi.^(pol:-1:0)';
BAgency_fun4 = matlabFunction(BAgency_fun4);

% Step 5: compute central finite difference
derivative_omega = @(z) (-BAgency_fun1(z) + 8*BAgency_fun2(z) - 8*BAgency_fun3(z) + BAgency_fun4(z)) / 12*h_omega;

figure('name', 'Numerical derivative omega agency')
fplot(derivative_omega, [0 1])

% Try with makima to get cross derivative:
xq = 0:0.0001:1;
m_agencyB_1 = makima(sol_BAgency_1.x, sol_BAgency_1.y(2,:), xq);
m_agencyB_2 = makima(sol_BAgency_2.x, sol_BAgency_2.y(2,:), xq);
m_agencyB_3 = makima(sol_BAgency_3.x, sol_BAgency_3.y(2,:), xq);
m_agencyB_4 = makima(sol_BAgency_4.x, sol_BAgency_4.y(2,:), xq);

derivative_omega_agency = (-m_agencyB_1 + 8*m_agencyB_2 - 8*m_agencyB_3 + m_agencyB_4) / 12*h_omega;

figure('name', '(MAKIMA) Numerical derivative omega agency')
plot(xq,derivative_omega_agency)
title('Numerical derivative omega (agency)')
saveas(gca,[pwd '/figures/agency_derivatives/makima_omega_agency.eps'],'epsc')

save('derivative_omega_agency.mat','derivative_omega_agency')


%% Compute Numerical Cross Derivative Tau Omega Agency Friction Model

parameters;
% Step 1: define step sizes
h_tau = tau*0.1; % change in lambda
h_omega = omega*0.5; % change in omega

% Step 2: compute the four value functions for different parameters
type = 'scaled'; % scaled: g(z,a)

% Filling in the formula for finite difference (see notes)
[~, sol_BAgency_1] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau+h_tau,sigma,a_bar,lambda,omega+h_omega,gamma,type);
[~, sol_BAgency_2] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau-h_tau,sigma,a_bar,lambda,omega+h_omega,gamma,type);
[~, sol_BAgency_3] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau+h_tau,sigma,a_bar,lambda,omega-h_omega,gamma,type);
[~, sol_BAgency_4] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau-h_tau,sigma,a_bar,lambda,omega-h_omega,gamma,type);

% Step 3: fit the derivative to obtain a comparable 'function'
pol = 17;
fit_dBAgency1 = polyfit(sol_BAgency_1.x, sol_BAgency_1.y(2,:),pol);
fit_dBAgency2 = polyfit(sol_BAgency_2.x, sol_BAgency_2.y(2,:),pol);
fit_dBAgency3 = polyfit(sol_BAgency_3.x, sol_BAgency_3.y(2,:),pol);
fit_dBAgency4 = polyfit(sol_BAgency_4.x, sol_BAgency_4.y(2,:),pol);

% Step 4: obtain functional form:
syms zi
BAgency_fun1 = fit_dBAgency1*zi.^(pol:-1:0)';
BAgency_fun1 = matlabFunction(BAgency_fun1);

BAgency_fun2 = fit_dBAgency2*zi.^(pol:-1:0)';
BAgency_fun2 = matlabFunction(BAgency_fun2);

BAgency_fun3 = fit_dBAgency3*zi.^(pol:-1:0)';
BAgency_fun3 = matlabFunction(BAgency_fun3);

BAgency_fun4 = fit_dBAgency4*zi.^(pol:-1:0)';
BAgency_fun4 = matlabFunction(BAgency_fun4);

% Step 5: compute central finite difference
cross_lambda_omega = @(z) (BAgency_fun1(z) - BAgency_fun2(z) - BAgency_fun3(z) + BAgency_fun4(z)) / 4*h_omega*h_tau;

figure('name', 'Numerical cross derivative')
fplot(cross_lambda_omega, [0 1])

% Try with makima to get derivative:
xq = 0:0.0001:1;
m_agencyB_1 = makima(sol_BAgency_1.x, sol_BAgency_1.y(2,:), xq);
m_agencyB_2 = makima(sol_BAgency_2.x, sol_BAgency_2.y(2,:), xq);
m_agencyB_3 = makima(sol_BAgency_3.x, sol_BAgency_3.y(2,:), xq);
m_agencyB_4 = makima(sol_BAgency_4.x, sol_BAgency_4.y(2,:), xq);

cross_tau_omega_agency = (m_agencyB_1 - m_agencyB_2 - m_agencyB_3 + m_agencyB_4) / 4*h_omega*h_tau;

figure('name', '(MAKIMA) Numerical cross derivative tau omega agency')
plot(xq,cross_tau_omega_agency)
title('Numerical Cross Derivative tau omega (agency)')
saveas(gca,[pwd '/figures/agency_derivatives/makima_cross_tau_omega_agency.eps'],'epsc')

save('cross_derivative_tau_omega_agency.mat','cross_tau_omega_agency')

% Comparison tau omega cross derivatives
load('cross_derivative_tau_omega_agency')
load('cross_derivative_tau_omega_fb')

figure; 
plot(xq,cross_tau_omega_fb)
hold on
plot(xq,cross_tau_omega_agency)
title('Comparison Cross Derivatives Tau Omega')
legend('First Best', 'Moral Hazard')
grid on

% Compute Numerical Cross Derivative Lambda Omega Agency Friction Model

parameters;
% Step 1: define step sizes
h_lambda = lambda*0.5; % change in lambda
h_omega = omega*0.5; % change in omega

% Step 2: compute the four value functions for different parameters
type = 'scaled'; % scaled: g(z,a)

% Filling in the formula for finite difference (see notes)
[~, sol_BAgency_1] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau,sigma,a_bar,lambda+h_lambda,omega+h_omega,gamma,type);
[~, sol_BAgency_2] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau,sigma,a_bar,lambda-h_lambda,omega+h_omega,gamma,type);
[~, sol_BAgency_3] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau,sigma,a_bar,lambda+h_lambda,omega-h_omega,gamma,type);
[~, sol_BAgency_4] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau,sigma,a_bar,lambda-h_lambda,omega-h_omega,gamma,type);

% Step 3: fit the derivative to obtain a comparable 'function'
pol = 17;
fit_dBAgency1 = polyfit(sol_BAgency_1.x, sol_BAgency_1.y(2,:),pol);
fit_dBAgency2 = polyfit(sol_BAgency_2.x, sol_BAgency_2.y(2,:),pol);
fit_dBAgency3 = polyfit(sol_BAgency_3.x, sol_BAgency_3.y(2,:),pol);
fit_dBAgency4 = polyfit(sol_BAgency_4.x, sol_BAgency_4.y(2,:),pol);

% Step 4: obtain functional form:
syms zi
BAgency_fun1 = fit_dBAgency1*zi.^(pol:-1:0)';
BAgency_fun1 = matlabFunction(BAgency_fun1);

BAgency_fun2 = fit_dBAgency2*zi.^(pol:-1:0)';
BAgency_fun2 = matlabFunction(BAgency_fun2);

BAgency_fun3 = fit_dBAgency3*zi.^(pol:-1:0)';
BAgency_fun3 = matlabFunction(BAgency_fun3);

BAgency_fun4 = fit_dBAgency4*zi.^(pol:-1:0)';
BAgency_fun4 = matlabFunction(BAgency_fun4);

% Step 5: compute central finite difference
cross_lambda_omega = @(z) (BAgency_fun1(z) - BAgency_fun2(z) - BAgency_fun3(z) + BAgency_fun4(z)) / 4*h_omega*h_lambda;

figure('name', 'Numerical cross derivative')
fplot(cross_lambda_omega, [0 1])

% Try with makima to get derivative:
xq = 0:0.0001:1;
m_agencyB_1 = makima(sol_BAgency_1.x, sol_BAgency_1.y(2,:), xq);
m_agencyB_2 = makima(sol_BAgency_2.x, sol_BAgency_2.y(2,:), xq);
m_agencyB_3 = makima(sol_BAgency_3.x, sol_BAgency_3.y(2,:), xq);
m_agencyB_4 = makima(sol_BAgency_4.x, sol_BAgency_4.y(2,:), xq);

cross_lambda_omega_agency = (m_agencyB_1 - m_agencyB_2 - m_agencyB_3 + m_agencyB_4) / 4*h_omega*h_lambda;

figure('name', '(MAKIMA) Numerical cross derivative lambda omega agency')
plot(xq,cross_lambda_omega_agency)
title('Numerical Cross Derivative (agency)')
saveas(gca,[pwd '/figures/agency_derivatives/makima_cross_agency.eps'],'epsc')

save('cross_derivative_agency.mat','cross_lambda_omega_agency')

%% Plot comparison cross derivative under fb and agency
load('cross_derivative_agency')
load('cross_derivative_fb')

figure; 
plot(xq,cross_lambda_omega_fb)
hold on
plot(xq,cross_lambda_omega_agency)
title('Comparison Cross Derivatives')
legend('First Best', 'Moral Hazard')
grid on

%% Comparison Derivatives Omega
load('derivative_omega_agency')
load('derivative_omega_fb')

figure; 
plot(xq,derivative_omega_fb)
hold on
plot(xq,derivative_omega_agency)
title('Comparison Derivative Omega')
legend('First Best', 'Moral Hazard')
grid on

%% Comparison Derivatives Lambda
load('derivative_lambda_agency')
load('derivative_lambda_fb')

figure; 
plot(xq,derivative_lambda_fb)
hold on
plot(xq,derivative_lambda_agency)
title('Comparison Derivative Lambda')
legend('First Best', 'Moral Hazard')
grid on

%% Comparison Derivatives Tau
load('derivative_tau_agency')
load('derivative_tau_fb')

figure; 
plot(xq,derivative_tau_fb)
hold on
plot(xq,derivative_tau_agency)
title('Comparison Derivative Tau')
legend('First Best', 'Moral Hazard')
grid on

%% Higher accuracy Cross derivative 
parameters;
% define step sizes
h_lambda = 0.01; % change in lambda
h_omega = 0.01; % change in omega
h=0.01;

% compute the four value functions for different parameters
type = 'scaled'; % scaled: g(z,a)

% Filling in the formula for finite difference (see notes)


%[~, sol_BAgency_1] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,sigma,a_bar,lambda+h_lambda,omega+h_omega,gamma,type);


% positive
f11 = fii(1,1, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);
f12 = fii(1,2, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);
f21 = fii(2,1, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);
f22 = fii(2,2, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);

% negative (m)
fm1m1 = fii(-1,-1, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);
fm1m2 = fii(-1,-2, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);
fm2m1 = fii(-2,-1, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);
fm2m2 = fii(-2,-2, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);

% positive and negative:
f1m1 = fii(1,-1, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);
fm11 = fii(-1,1, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);

fm21 = fii(-2,1, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);
f1m2 = fii(1,-2, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);

f2m1 = fii(2,-1, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);
fm12 = fii(-1,2, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);

f2m2 = fii(2,-2, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);
fm22 = fii(-2,2, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);

%%%%%%%%%%

% positive
f11 = makima(f11.x, f11.y(2,:), xq);
f12 = makima(f12.x, f12.y(2,:), xq);
f21 = makima(f21.x, f21.y(2,:), xq);
f22 = makima(f22.x, f22.y(2,:), xq);

% negative (m)
fm1m1 = makima(fm1m1.x, fm1m1.y(2,:), xq);
fm1m2 = makima(fm1m2.x, fm1m2.y(2,:), xq);
fm2m1 = makima(fm2m1.x, fm2m1.y(2,:), xq);
fm2m2 = makima(fm2m2.x, fm2m2.y(2,:), xq);

% positive and negative:
f1m1 = makima(f1m1.x, f1m1.y(2,:), xq);
fm11 = makima(fm11.x, fm11.y(2,:), xq);

fm21 = makima(fm21.x, fm21.y(2,:), xq);
f1m2 = makima(f1m2.x, f1m2.y(2,:), xq);

f2m1 = makima(f2m1.x, f2m1.y(2,:), xq);
fm12 = makima(fm12.x, fm12.y(2,:), xq);

f2m2 = makima(f2m2.x, f2m2.y(2,:), xq);
fm22 = makima(fm22.x, fm22.y(2,:), xq);

cross_lambda_omega = (8*(f1m2 + f2m1 + fm21 + fm12) - 8*(fm1m2 + fm2m1 + f12 + f21) ...
                       - (f2m2 + fm22 - fm2m2 - f22) + 64*(fm1m1 + f11 - f1m1 - fm11)) / 144*h^2;


figure('name', 'Cross Agency 4th order error')
plot(xq,cross_lambda_omega)
title('Agency Cross')
saveas(gca,[pwd '/figures/agency_derivatives/makima_cross_agency_fourth_order.eps'],'epsc')



%% Comparison fb and agency firm value
parameters;
type='scaled';

[sol_AShock, sol_BShock] = solve_ODE_SAR_FB(r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,type);
[sol_A_Agency, sol_B_Agency] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type);

figure;
plot(sol_BShock.x, sol_BShock.y(1,:))
hold on
plot(sol_B_Agency.x, sol_B_Agency.y(1,:))
legend('First Best', 'Moral Hazard')
title('before')

figure;
plot(sol_AShock.x, sol_AShock.y(1,:))
hold on
plot(sol_A_Agency.x, sol_A_Agency.y(1,:))
legend('First Best', 'Moral Hazard')
title('after')

%% Some basic functions

function y = obtain_a(sol_y_prime,a_bar)
    for i=1:length(sol_y_prime)
        if abs(sol_y_prime(i)) > a_bar
            sol_y_prime(i) = sign(sol_y_prime(i))*a_bar;
        end
    end
    y = sol_y_prime;
end


function sol_ii = fii(k_l, k_o, h, r,mu_B,mu_G,tau,sigma,a_bar,lambda,omega,gamma,type)
    [~, sol_ii] = solve_ODE_SAR_AGENCY(r,mu_B,mu_G,tau,sigma,a_bar,lambda+h*k_l,omega+h*k_o,gamma,type);
end
