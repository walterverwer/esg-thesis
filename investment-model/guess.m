function y = guess(z,y0_guess,r ,A_G, A_B ,delta,theta,expr)
    if z==0
        y=y0_guess;
        return
    end
    ode = @(z,y) ode_v_fb(z,y,r, A_G, A_B ,delta,theta,expr);

    options = odeset(RelTol=1e-4,AbsTol=1e-4);
    sol = ode89(ode, [0 z], y0_guess,options);
    y = sol.y(:,end);
end

