function [sol,p_B,p_G,i_B,i_G] = bv_solver(AB, AG)
% Solve the boundary value problem for the FB, model type a la Eberly and
% Wang (2010). Given the singularity, I do a series expansion for
% values of z close to 0 and 1.
parameters;
A_B = AB;
A_G = AG;

[p_B,p_G,i_B,i_G] = obtain_boundary_values(r, A_G, A_B ,delta,theta);

syms x
expr = matlabFunction(...
    taylor( 2/( (sigma_B^2+sigma_G^2)* x^2*(1-x)^2 ), x, 'ExpansionPoint',.5,'Order',4));

ode = @(z,y) ode_v_fb(z,y,r,A_G, A_B ,delta,theta,expr);
bc = @(ya, yb) bc_fb(ya, yb, p_B, p_G);


% obtain init
xmesh = linspace(0,1,1000);
y0_guess = [p_B;0]; % start is known, V' is unknown. Guess whole line
guessFun = @(z) guess(z,y0_guess,r ,A_G, A_B ,delta,theta,expr);
solinit = bvpinit(xmesh, guessFun);

bvpoptions = bvpset(Stats="on",Nmax=10000,AbsTol=1e-2,RelTol=1e-2);
sol = bvp5c(ode,bc,solinit,bvpoptions);
end