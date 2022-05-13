tic

% parameters
parameters.r        = 0.05;         % discounting
parameters.sigma    = 0.25;         % volatility
parameters.mu_G     = 0.10;         % productivity of green capital
parameters.mu_B     = 0.20;         % productivity of brown capital
parameters.gamma    = 5;            % risk aversion coefficient
parameters.theta    = 10;           % constant scaling adjustment cost
parameters.a_bar    = 0.4;          % maximum value of effort to change z
parameters.lambda   = 0.015;        % poisson intensity parameter
parameters.omega    = 0.010;        % green preference parameter of principal           
parameters.tau      = 0.20;         % tax on brown/ subsidy on green
parameters.xi       = 0.02;         % green preference parameter of agent
parameters.pol      = 17;           % polinomial length

% construct waitbar
H = multiwaitbar(2,[0 0],{'Types Progress','Derivatives Progress'});

% model types
types = ["firstBest", "agency"];

% parameters for derivatives
paramsDerivatives = ["lambda", "tau", "omega"];

% solution cells
modelSolutions = {};
modelDerivatives = {};
modelCrossDerivatives = {};

L = length(paramsDerivatives);

% solve model, derivatives and cross derivatives for effort
parfor i=1:length(types)
    model = transitionModel( types(i), parameters );
    [sol_A, sol_B] = model.solve;
    modelSolutions{i} = model;
    
    for j=1:L
        % update waitbar
        multiwaitbar(2,[i j],{'Computing','Computing'},H);

        
            % compute single derivative of effort wrt parameter j
            model = transitionModel( types(i), parameters );
            var1 = paramsDerivatives(j); relativeStep1 = 0.1;
            [derivativeAfterShock, derivativeBeforeShock] = ...
            model.derivative(var1, relativeStep1);
            
            % write single derivative to solution cell
            modelDerivatives{i,j} = [derivativeAfterShock, derivativeBeforeShock];

        if paramsDerivatives(j) ~= "omega"
            % compute cross derivative of effort wrt parameters j (not omega) and omega
            model = transitionModel( types(i), parameters );
            var1 = paramsDerivatives(j); relativeStep1 = 0.1;
            var2 = "omega"; relativeStep2 = 0.1;
            [crossDerivativeAfterShock, crossDerivativeBeforeShock] = ...
                model.crossDerivative(var1, var2, relativeStep1, relativeStep2);
    
            % write cross derivative to solution cell
            modelCrossDerivatives{i,j} = [crossDerivativeAfterShock, crossDerivativeBeforeShock];
        end
    end
end

% close waitbar
delete(H.figure)
clear('H')

% save workspace
save('results.mat')

toc
%% Plot all results

% firm value plots
figure;
plot(modelSolutions{1}.sol_BShock.x, modelSolutions{1}.sol_BShock.y(1,:))
hold on
plot(modelSolutions{2}.sol_BShock.x, modelSolutions{2}.sol_BShock.y(1,:))
title('Firm Value Before Shock', 'interpreter','latex')
grid on
legend('First Best', 'Moral Hazard')
saveas(gca,[pwd '/results/firmValue.eps'],'epsc')

% firm value derivatives
figure;
plot(modelSolutions{1}.sol_BShock.x, modelSolutions{1}.sol_BShock.y(2,:))
hold on
plot(modelSolutions{2}.sol_BShock.x, modelSolutions{2}.sol_BShock.y(2,:))
title('Firm Value Before Shock', 'interpreter','latex')
grid on
legend('First Best', 'Moral Hazard')
saveas(gca,[pwd '/results/firmValueDerivative.eps'],'epsc')