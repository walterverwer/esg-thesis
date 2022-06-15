%% Environemntal Investment Model
% Solve main model. Cross derivatives for omega
tic

% parameters
parameters.r        = 0.05;         % discounting
parameters.sigma    = 1.00;         % volatility
parameters.mu_G     = 0.20;         % productivity of green capital
parameters.mu_B     = 0.30;         % productivity of brown capital
parameters.gamma    = 5;            % risk aversion coefficient
parameters.theta    = 10;           % constant scaling adjustment cost
parameters.a_bar    = 0.4;          % maximum value of effort to change z
parameters.lambda   = 0.020;        % poisson intensity parameter
parameters.omega    = 0.010;        % green preference parameter of principal           
parameters.tau      = 0.30;         % tax on brown/ subsidy on green
parameters.xi       = 0.00;         % green preference parameter of agent (on z)
parameters.pol      = 15;           % polinomial length
parameters.phi      = 0.00;         % manager's preference for doing green effort (on a)

% model types
types = "firstBest";

% parameters for derivatives
paramsDerivatives = ["lambda", "tau", "omega"];
bigParams = ["sigma", "theta", "gamma"];
derivativeType = "effort";
L = length(paramsDerivatives);

% solution cells
modelSolutions = {};
modelDerivatives = {};
modelCrossDerivatives = {};

% param to take cross derivative with
crossDerivParam = "omega";

% solve model, derivatives and cross derivatives for effort
parfor i=1:length(types)
    model = transitionModel( types(i), parameters );
    [sol_A, sol_B] = model.solve;
    modelSolutions{i} = model;
    
    for j=1:L        
        % compute single derivative of effort wrt parameter j
        model = transitionModel( types(i), parameters );
        var1 = paramsDerivatives(j); 
        
        % some params are rather large, so need smaller step size
        if  any(contains(bigParams, paramsDerivatives(j)))
            relativeStep1 = 0.0025;
        else
            relativeStep1 = 0.025;
        end

        [derivativeAfterShock, derivativeBeforeShock] = ...
        model.derivative(var1, relativeStep1, derivativeType);
        
        % write single derivative to solution cell
        modelDerivatives{i,j} = [derivativeAfterShock, derivativeBeforeShock];

        if paramsDerivatives(j) ~= crossDerivParam
            % compute cross derivative of effort wrt parameters j (not omega) and omega
            model = transitionModel( types(i), parameters );
            var1 = paramsDerivatives(j); 
            
            % some params are rather large, so need smaller step size
            if  any(contains(bigParams, paramsDerivatives(j)))
                relativeStep1 = 0.0025;
            else
                relativeStep1 = 0.025;
            end
            
            var2 = crossDerivParam; relativeStep2 = 0.025;

            [crossDerivativeAfterShock, crossDerivativeBeforeShock] = ...
            model.crossDerivative(var1, var2, relativeStep1, relativeStep2);
    
            % write cross derivative to solution cell
            modelCrossDerivatives{i,j} = [crossDerivativeAfterShock, crossDerivativeBeforeShock];
        end
    end
end


% write results to structure
results = struct;
results.solutions = modelSolutions;
results.derivatives = modelDerivatives;
results.crossDerivatives = modelCrossDerivatives;

% save results
save('transitionModel.mat')

toc

%% Plot results: firm value + effort
load('transitionModel.mat')

% create tile plots
t=tiledlayout(2,2);

% Tile 1
nexttile
plot(modelSolutions{1}.sol_BShock.x, modelSolutions{1}.sol_BShock.y(1,:), 'b')
hold on
title('(a) Firm Value Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 2
nexttile
plot(modelSolutions{1}.sol_AShock.x, modelSolutions{1}.sol_AShock.y(1,:), 'b')
hold on
title('(b) Firm Value After Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 3
nexttile
plot(modelSolutions{1}.sol_BShock.x, modelSolutions{1}.effort.aBeforeShock, 'b')
hold on
title('(c) Effort Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 4
nexttile
plot(modelSolutions{1}.sol_AShock.x, modelSolutions{1}.effort.aAfterShock, 'b')
hold on
title('(d) Effort After Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

exportgraphics(t,[pwd '/results/transitionModel/firmValueWithEffort.eps'], 'ContentType', 'vector')

%% derivatives: paramsDerivatives = ["lambda", "tau", "omega"];

x = linspace(0,1,10000);

t=tiledlayout(1,3);

nexttile
plot(x,modelDerivatives{1,1}(:,2), 'b')
hold on
title('(a) $\frac{\partial a}{\partial \lambda}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,2}(:,2), 'b')
hold on
title('(b) $\frac{\partial a}{\partial \tau}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,3}(:,2), 'b')
hold on
title('(c) $\frac{\partial a}{\partial \omega}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

exportgraphics(t,[pwd '/results/transitionModel/derivatives.eps'], 'ContentType', 'vector')

%% cross derivatives all for omega

x = linspace(0,1,10000);

t=tiledlayout(1,2);

nexttile
plot(x,modelCrossDerivatives{1,1}(:,2), 'b')
hold on
title('(a) $\frac{\partial^2 a}{\partial \lambda \partial \omega}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelCrossDerivatives{1,2}(:,2), 'b')
hold on
title('(b) $\frac{\partial^2 a}{\partial \tau \partial \omega}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

exportgraphics(t,[pwd '/results/transitionModel/crossDerivativesOmega.eps'], 'ContentType', 'vector')

