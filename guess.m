function y = guess(z,y0_guess,r,sigma_G, sigma_B ,A_G, A_B ,delta,theta)
    if z==0
        y=y0_guess;
        return
    end
    ode = @(z,y) ode_v_fb(z,y,r,sigma_G, sigma_B ,A_G, A_B ,delta,theta);

    options = odeset(RelTol=1e-2,AbsTol=1e-2,MassSingular="yes");
    sol = ode45(ode, [0 z], y0_guess,options);
    y = sol.y(:,end);
end

