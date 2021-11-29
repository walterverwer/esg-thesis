clear all, close all

%% Boundary value problem solver:
parameters;

%S = [0 0; 0 2/(sigma_G^2+sigma_B^2)];
%options = bvpset('SingularTerm',S);

[p_B,p_G,i_B,i_G] = obtain_boundary_values(r, A_G, A_B ,delta,theta);

ode = @(z,y) ode_v_fb(z,y,r,sigma_G, sigma_B ,A_G, A_B ,delta,theta);
bc = @(ya, yb) bc_fb(ya, yb, p_B, p_G);


% obtain init
lmd = 0;
xmesh = linspace(0+lmd,1-lmd,1000);
y0_guess = [p_B;0]; % start is known, V' is unknown. Guess whole line
guessFun = @(z) guess(z,y0_guess,r,sigma_G, sigma_B ,A_G, A_B ,delta,theta);
solinit = bvpinit(xmesh, guessFun);

bvpoptions = bvpset(Stats="on",Nmax=10000,AbsTol=1e-3,RelTol=1e-3);
fb = bvp4c(ode,bc,solinit,bvpoptions);
found_y = fb.y(:,1);

options = odeset(RelTol=1e-4,AbsTol=1e-4,Events=@isOptim);
sol = ode45(ode, [0 1], found_y, options);

%% Plotting


plot(fb.x,fb.y(1,:))
hold on
%plot(sol.x,sol.y(1,:))
grid on


%%



