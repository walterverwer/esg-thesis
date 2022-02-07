function y = guess(z,y0_guess,r ,mu_G, mu_B,theta,gamma,sigma,expr,exprTheta, exprThetaMinus)
    if z==0
        y=y0_guess;
        return
    end
    ode_fun = @(z,y) ode(z,y,r,mu_B,mu_G,theta,gamma,sigma,expr, exprTheta, exprThetaMinus);

    options = odeset(RelTol=1e-4,AbsTol=1e-4);
    sol = ode89(ode_fun, [0 z], y0_guess,options);
    y = sol.y(:,end);
end
