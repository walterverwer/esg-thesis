function y = guessAShock(z,y0_guess,r,mu_G,xi,expr,a_bar)
    if z==0
        y=y0_guess;
        return
    end
    ode_fun = @(z,y) odeAShock(z,y,r,mu_G,xi,expr,a_bar);

    options = odeset(RelTol=1e-4,AbsTol=1e-4);
    sol = ode89(ode_fun, [0 z], y0_guess,options);
    y = sol.y(:,end);
end


