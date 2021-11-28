clear all, close all

%% Basic parameters
parameters;
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

y0 = ones(2,1);

% Return the function V(z) and plot it
[z,Vz] = ode45(@ode_v_fb ,[0 1], y0, options, r,sigma_G, sigma_B ,A_G, A_B ,delta,theta);
%[z,Vz] = bvp4c(@ode_v_fb , @bc_fb, solinit_fb);
plot(z,Vz(:,1));