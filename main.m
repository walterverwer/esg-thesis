clear all, close all

%% Basic parameters
parameters;
options = odeset('RelTol',1e-8,'AbsTol',1e-8, 'MassSingular','yes');

%% HJB fb

%y0 = ones(2,1);

% Return the function V(z) and plot it, THIS IS INCORRECT! This is an
% initial value problem, I need to do 2 boundary problem!
[p_B, p_G] = obtain_boundary_values(r, A_G, A_B ,delta,theta);
[z,Vz] = ode45(@ode_v_fb ,[0 1], [p_B 1], options, r,sigma_G, sigma_B ,A_G, A_B ,delta,theta);

plot(z,Vz(:,1));




%% Boundary value problem solver:
parameters;

S = [0 0; 0 2/(sigma_G^2+sigma_B^2)];
options = bvpset('SingularTerm',S);


xmesh = linspace(0,1,10);
solinit_fb = bvpinit(xmesh, @guess);
%sol_fb = bvp4c(@ode_v_fb , @bc_fb, solinit_fb, options, r,sigma_G, sigma_B ,A_G, A_B ,delta,theta);



% Use a boundary value problem solver:
%[p_B, p_G] = obtain_boundary_values(A_G, A_B ,delta,theta);
%xmesh = linspace(0,1,100);
%solinit_fb = bvpinit(xmesh, @guess, r,sigma_G, sigma_B ,A_G, A_B ,delta,theta);
%sol_fb = bvp4c(@ode_v_fb , @bc_fb, solinit_fb);
%plot(sol_fb.x,sol_fb.y);






