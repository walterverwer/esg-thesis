function [sol,p_B,p_G,i_B,i_G] = bv_solver(AB, AG)
% Solve the boundary value problem for the FB, model type a la Eberly and
% Wang (2010). Given the singularity, I do a series expansion for
% values of z close to 0 and 1.
parameters;
A_B = AB;
A_G = AG;

[p_B,p_G,i_B,i_G] = obtain_boundary_values(r, A_G, A_B ,delta,theta);

ode = @(z,y) ode_v_fb(z,y,r,sigma_G, sigma_B ,A_G, A_B ,delta,theta);
bc = @(ya, yb) bc_fb(ya, yb, p_B, p_G);


% obtain init
lmd = 0;
xmesh = linspace(0+lmd,1-lmd,1000);
y0_guess = [p_B;0]; % start is known, V' is unknown. Guess whole line
guessFun = @(z) guess(z,y0_guess,r,sigma_G, sigma_B ,A_G, A_B ,delta,theta);
solinit = bvpinit(xmesh, guessFun);

bvpoptions = bvpset(Stats="on",Nmax=10000,AbsTol=1e-2,RelTol=1e-2);
sol = bvp5c(ode,bc,solinit,bvpoptions);
end