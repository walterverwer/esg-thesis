function y = guessFb(z,y0_guess,r,mu_G,mu_B,expr,exprTheta,exprThetaMinus,a_bar,option)
    if z==0
        y=y0_guess;
        return
    end
    ode_fun = @(z,y) odeFb(z,y,r,mu_B,mu_G,expr,exprTheta,exprThetaMinus,a_bar,option);

    options = odeset(RelTol=1e-4,AbsTol=1e-4);
    sol = ode89(ode_fun, [0 z], y0_guess,options);
    y = sol.y(:,end);
end
