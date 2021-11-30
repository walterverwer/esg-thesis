clear all, close all

%% Boundary value problem solver:
parameters;
[sol,p_B,p_G,i_B,i_G] = bv_solver(A_B,A_G);

% Plot results
plot(sol.x,sol.y(1,:))
legend(sprintf('A_B=A_G=%.1f',A_B),'Location','southeast')
grid on
saveas(gca,'fb.eps','epsc')

%% Plotting: plot for different As
parameters;

% Baseline: equal A
AB1 = A_B;
AG1 = A_G;
[sol1,~,~,~,~] = bv_solver(AB1, AG1);

% A_B>A_G
AB2 = 0.45;
AG2 = 0.35;
[sol2,~,~,~,~] = bv_solver(AB2, AG2);

% A_B<A_G
AB3 = 0.35;
AG3 = 0.45;
[sol3,~,~,~,~] = bv_solver(AB3, AG3);


plot(sol1.x,sol1.y(1,:))
hold on
plot(sol2.x,sol2.y(1,:))
hold on
plot(sol3.x,sol3.y(1,:))
legend(sprintf('A_B=A_G=%.2f',AB1), sprintf('A_B=%.2f & A_G=%.2f',AB2, AG2), sprintf('A_B=%.2f & A_G=%.2f',AB3,AG3),'Location','southeast')
grid on
saveas(gca,'fb_comparison.eps','epsc')
















