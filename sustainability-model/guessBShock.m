function y = guessBShock(z,y0_guess,r,mu_B,mu_G,omega,lambda,expr,a_bar,AShock_fun,type)
    if z==0
        y=y0_guess;
        return
    end
    ode_fun = @(z,y) odeBShock(z,y,r,mu_B,mu_G,omega,lambda,expr,a_bar,AShock_fun,type);

    options = odeset(RelTol=1e-4,AbsTol=1e-4);
    sol = ode89(ode_fun, [0 z], y0_guess,options);
    y = sol.y(:,end);
end


