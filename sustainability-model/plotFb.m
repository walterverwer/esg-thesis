function plotFb(x,y,mu_G,mu_B,option)

if isequal(option,'scaled')
    if mu_B==mu_G
        plot(x,y)
        grid on
        title('First best, $\mu^G = \mu^B$ for $g(z,a)=\frac{a^2z(1-z)}{2}$', 'interpreter','latex')
        saveas(gca,'/figures/various_gFun_plots/first_best/sol_equal_fb_scaled.eps','epsc')
    else
        plot(x,y)
        grid on
        title('First best, $\mu^G \neq \mu^B$ for $g(z,a)=\frac{a^2z(1-z)}{2}$', 'interpreter','latex')
        saveas(gca,'/figures/various_gFun_plots/first_best/sol_unequal_fb_scaled.eps','epsc')
    end
end


if isequal(option,'unscaled')
    if mu_B==mu_G
        plot(x,y)
        grid on
        title('First best, $\mu^G = \mu^B$ for $g(z,a)=\frac{a^2}{2}$', 'interpreter','latex')
        saveas(gca,'/figures/various_gFun_plots/first_best/sol_equal_fb_unscaled.eps','epsc')
    else
        plot(x,y)
        grid on
        title('First best, $\mu^G \neq \mu^B$ for $g(z,a)=\frac{a^2}{2}$', 'interpreter','latex')
        saveas(gca,'/figures/various_gFun_plots/first_best/sol_unequal_fb_unscaled.eps','epsc')
    end
end


if isequal(option,'theta')
    if mu_B==mu_G
        plot(x,y)
        grid on
        title('First best, $\mu^G = \mu^B$ for $g(z,a)=\frac{a^2z^\theta(1-z)^\theta}{2}$', 'interpreter','latex')
        saveas(gca,'/figures/various_gFun_plots/first_best/sol_equal_fb_theta.eps','epsc')
    else
        plot(x,y)
        grid on
        title('First best, $\mu^G \neq \mu^B$ for $g(z,a)=\frac{a^2z^\theta(1-z)^\theta}{2}$', 'interpreter','latex')
        saveas(gca,'/figures/various_gFun_plots/first_best/sol_unequal_fb_theta.eps','epsc')
    end
end

end

