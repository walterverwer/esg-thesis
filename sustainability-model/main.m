%% Environmental Transition Model Solver

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
types = ["firstBest", "agency"];

% parameters for derivatives
paramsDerivatives = ["lambda", "tau", "omega", "sigma"];
derivativeType = "effort";
L = length(paramsDerivatives);

% solution cells
modelSolutions = {};
modelDerivatives = {};
modelCrossDerivatives = {};

% solve model, derivatives and cross derivatives for effort
parfor i=1:length(types)
    model = transitionModel( types(i), parameters );
    [sol_A, sol_B] = model.solve;
    modelSolutions{i} = model;
    
    for j=1:L        
        % compute single derivative of effort wrt parameter j
        model = transitionModel( types(i), parameters );
        var1 = paramsDerivatives(j); 
        
        % gamma is rather large, so need smaller step size
        if paramsDerivatives(j) == "sigma"
            relativeStep1 = 0.0025;
        else
            relativeStep1 = 0.025;
        end

        [derivativeAfterShock, derivativeBeforeShock] = ...
        model.derivative(var1, relativeStep1, derivativeType);
        
        % write single derivative to solution cell
        modelDerivatives{i,j} = [derivativeAfterShock, derivativeBeforeShock];

        if paramsDerivatives(j) ~= "omega"
            % compute cross derivative of effort wrt parameters j (not omega) and omega
            model = transitionModel( types(i), parameters );
            var1 = paramsDerivatives(j); 
            
            % gamma is rather large, so need smaller step size
            if paramsDerivatives(j) == "sigma"
                relativeStep1 = 0.0025;
            elseif paramsDerivatives(j) == "theta"
                relativeStep1 = 0.0025
            else
                relativeStep1 = 0.025;
            end
            
            var2 = "omega"; relativeStep2 = 0.025;

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

% save workspace
save('results_noPhi_noXi.mat')

toc


%% Plot all results without xi and phi
load('results_noPhi_noXi.mat')

% firm value + effort
% create tile plots
t=tiledlayout(2,2);

% Tile 1
nexttile
plot(modelSolutions{1}.sol_BShock.x, modelSolutions{1}.sol_BShock.y(1,:), 'b')
hold on
plot(modelSolutions{2}.sol_BShock.x, modelSolutions{2}.sol_BShock.y(1,:), 'r')
title('(a) Firm Value Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 2
nexttile
plot(modelSolutions{1}.sol_AShock.x, modelSolutions{1}.sol_AShock.y(1,:), 'b')
hold on
plot(modelSolutions{2}.sol_AShock.x, modelSolutions{2}.sol_AShock.y(1,:), 'r')
title('(b) Firm Value After Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 3
nexttile
plot(modelSolutions{1}.sol_BShock.x, modelSolutions{1}.effort.aBeforeShock, 'b')
hold on
plot(modelSolutions{2}.sol_BShock.x, modelSolutions{2}.effort.aBeforeShock, 'r')
title('(c) Effort Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 4
nexttile
plot(modelSolutions{1}.sol_AShock.x, modelSolutions{1}.effort.aAfterShock, 'b')
hold on
plot(modelSolutions{2}.sol_AShock.x, modelSolutions{2}.effort.aAfterShock, 'r')
title('(d) Effort After Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

leg = legend('First Best', 'Moral Hazard');
leg.Layout.Tile = 'south';
%saveas(gca,[pwd '/results/tileFirmValueWithEffort.eps'],'epsc')

exportgraphics(t,[pwd '/results/noPhiXi/tileFirmValueWithEffort.eps'], 'ContentType', 'vector')

%% derivatives: paramsDerivatives = ["lambda", "tau","omega", "gamma"];

x = linspace(0,1,10000);

t=tiledlayout(2,2);

nexttile
plot(x,modelDerivatives{1,1}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,1}(:,2), 'r')
title('(a) $\frac{\partial a}{\partial \lambda}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,2}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,2}(:,2), 'r')
title('(b) $\frac{\partial a}{\partial \tau}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,3}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,3}(:,2), 'r')
title('(c) $\frac{\partial a}{\partial \omega}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,4}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,4}(:,2), 'r')
title('(d) $\frac{\partial a}{\partial \sigma}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

leg = legend('First Best', 'Moral Hazard');
leg.Layout.Tile = 'south';
%saveas(gca,[pwd '/results/tileDerivatives_noXiPhi.eps'],'epsc')

exportgraphics(t,[pwd '/results/noPhiXi/tileDerivatives_noXiPhi.eps'], 'ContentType', 'vector')

%% cross derivatives

x = linspace(0,1,10000);

% cross derivatives: paramsDerivatives = ["lambda", "tau","omega", "gamma"];
t=tiledlayout(2,2);

nexttile
plot(x,modelCrossDerivatives{1,1}(:,2), 'b')
hold on
plot(x,modelCrossDerivatives{2,1}(:,2), 'r')
title('(a) $\frac{\partial^2 a}{\partial \lambda \partial \omega}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelCrossDerivatives{1,2}(:,2), 'b')
hold on
plot(x,modelCrossDerivatives{2,2}(:,2), 'r')
title('(b) $\frac{\partial^2 a}{\partial \tau \partial \omega}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile([1 2])
plot(x,modelCrossDerivatives{1,4}(:,2), 'b')
hold on
plot(x,modelCrossDerivatives{2,4}(:,2), 'r')
title('(c) $\frac{\partial^2 a}{\partial \sigma \partial \omega}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

leg = legend('First Best', 'Moral Hazard');
leg.Layout.Tile = 'south';
%saveas(gca,[pwd '/results/tileCrossDerivatives_noXiPhi.eps'],'epsc')

exportgraphics(t,[pwd '/results/noPhiXi/tileCrossDerivatives_noXiPhi.eps'], 'ContentType', 'vector')



%% Do analysis with xi

tic

% parameters
parameters.r        = 0.05;         % discounting
parameters.sigma    = 1.00;         % volatility
parameters.mu_G     = 0.20;         % productivity of green capital
parameters.mu_B     = 0.30;         % productivity of brown capital
parameters.gamma    = 5;            % risk aversion coefficient
parameters.theta    = 10;           % constant scaling adjustment cost
parameters.a_bar    = 0.4;          % maximum value of effort to change z
parameters.lambda   = 0.025;        % poisson intensity parameter
parameters.omega    = 0.005;        % green preference parameter of principal           
parameters.tau      = 0.25;         % tax on brown/ subsidy on green
parameters.xi       = 0.005;         % green preference parameter of agent (on z)
parameters.pol      = 15;           % polinomial length
parameters.phi      = 0.00;         % manager's preference for doing green effort (on a)

% model types
types = ["firstBest", "agency"];

% parameters for derivatives
paramsDerivatives = ["lambda", "tau", "omega", "gamma", "xi"];
derivativeType = "effort";
L = length(paramsDerivatives);

% solution cells
modelSolutions = {};
modelDerivatives = {};
modelCrossDerivatives = {};

% solve model, derivatives and cross derivatives for effort
parfor i=1:length(types)
    model = transitionModel( types(i), parameters );
    [sol_A, sol_B] = model.solve;
    modelSolutions{i} = model;
    
    for j=1:L        
        % compute single derivative of effort wrt parameter j
        model = transitionModel( types(i), parameters );
        var1 = paramsDerivatives(j); 
        
        % gamma is rather large, so need smaller step size
        if paramsDerivatives(j) == "gamma"
            relativeStep1 = 0.001;
        else
            relativeStep1 = 0.01;
        end

        [derivativeAfterShock, derivativeBeforeShock] = ...
        model.derivative(var1, relativeStep1, derivativeType);
        
        % write single derivative to solution cell
        modelDerivatives{i,j} = [derivativeAfterShock, derivativeBeforeShock];

        if paramsDerivatives(j) ~= "omega"
            % compute cross derivative of effort wrt parameters j (not omega) and omega
            model = transitionModel( types(i), parameters );
            var1 = paramsDerivatives(j); 
            
            % gamma is rather large, so need smaller step size
            if paramsDerivatives(j) == "gamma"
                relativeStep1 = 0.0001;
            else
                relativeStep1 = 0.01;
            end
            
            var2 = "omega"; relativeStep2 = 0.01;

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

% save workspace
save('results_WithXi.mat')

toc

%% Plot all results with xi
load('results_WithXi.mat')

% firm value + effort
% create tile plots
t=tiledlayout(2,2);

% Tile 1
nexttile
plot(modelSolutions{1}.sol_BShock.x, modelSolutions{1}.sol_BShock.y(1,:), 'b')
hold on
plot(modelSolutions{2}.sol_BShock.x, modelSolutions{2}.sol_BShock.y(1,:), 'r')
title('(a) Firm Value Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 2
nexttile
plot(modelSolutions{1}.sol_AShock.x, modelSolutions{1}.sol_AShock.y(1,:), 'b')
hold on
plot(modelSolutions{2}.sol_AShock.x, modelSolutions{2}.sol_AShock.y(1,:), 'r')
title('(b) Firm Value After Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 3
nexttile
plot(modelSolutions{1}.sol_BShock.x, modelSolutions{1}.effort.aBeforeShock, 'b')
hold on
plot(modelSolutions{2}.sol_BShock.x, modelSolutions{2}.effort.aBeforeShock, 'r')
title('(c) Effort Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 4
nexttile
plot(modelSolutions{1}.sol_AShock.x, modelSolutions{1}.effort.aAfterShock, 'b')
hold on
plot(modelSolutions{2}.sol_AShock.x, modelSolutions{2}.effort.aAfterShock, 'r')
title('(d) Effort After Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

leg = legend('First Best', 'Moral Hazard');
leg.Layout.Tile = 'south';
% saveas(gca,[pwd '/results/tileFirmValueWithEffort.eps'],'epsc')

exportgraphics(t,[pwd '/results/withPhiXi/tileFirmValueWithEffortXi.eps'], 'ContentType', 'vector')

x = linspace(0,1,10000);


%% derivatives: paramsDerivatives = ["lambda", "tau", "omega", "gamma", "xi"];
t=tiledlayout(3,2);

nexttile
plot(x,modelDerivatives{1,1}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,1}(:,2), 'r')
title('(a) $\frac{\partial a}{\partial \lambda}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,2}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,2}(:,2), 'r')
title('(b) $\frac{\partial a}{\partial \tau}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,3}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,3}(:,2), 'r')
title('(c) $\frac{\partial a}{\partial \omega}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,4}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,4}(:,2), 'r')
title('(d) $\frac{\partial a}{\partial \gamma}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile([1 2])
plot(x,modelDerivatives{1,5}(:,2), 'b') %= 0
hold on
plot(x,modelDerivatives{2,5}(:,2), 'r')
title('(e) $\frac{\partial a}{\partial \xi}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

leg = legend('First Best', 'Moral Hazard');
leg.Layout.Tile = 'south';
%saveas(gca,[pwd '/results/tileDerivatives_noXiPhi.eps'],'epsc')

exportgraphics(t,[pwd '/results/withPhiXi/tileDerivativesXi.eps'], 'ContentType', 'vector')


%% cross derivatives: paramsDerivatives = ["lambda", "tau", "omega", "gamma", "xi"];
x = linspace(0,1,10000);
t=tiledlayout(2,2);

nexttile
plot(x,modelCrossDerivatives{1,1}(:,2), 'b')
hold on
plot(x,modelCrossDerivatives{2,1}(:,2), 'r')
title('$\frac{\partial^2 a}{\partial \lambda \partial \omega}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelCrossDerivatives{1,2}(:,2), 'b')
hold on
plot(x,modelCrossDerivatives{2,2}(:,2), 'r')
title('$\frac{\partial^2 a}{\partial \tau \partial \omega}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelCrossDerivatives{1,4}(:,2), 'b')
hold on
plot(x,modelCrossDerivatives{2,4}(:,2), 'r')
title('$\frac{\partial^2 a}{\partial \gamma \partial \omega}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelCrossDerivatives{1,5}(:,2), 'b')
hold on
plot(x,modelCrossDerivatives{2,5}(:,2), 'r')
title('$\frac{\partial^2 a}{\partial \xi \partial \omega}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")


leg = legend('First Best', 'Moral Hazard');
leg.Layout.Tile = 'south';
%saveas(gca,[pwd '/results/tileCrossDerivatives_noXiPhi.eps'],'epsc')

exportgraphics(t,[pwd '/results/withPhiXi/tileCrossDerivativesXi.eps'], 'ContentType', 'vector')


%% Quick compution of firm value + effort (for testing)

tic

% parameters
parameters.r        = 0.05;         % discounting
parameters.sigma    = 1.00;         % volatility
parameters.mu_G     = 0.20;         % productivity of green capital
parameters.mu_B     = 0.30;         % productivity of brown capital
parameters.gamma    = 10;            % risk aversion coefficient
parameters.theta    = 10;           % constant scaling adjustment cost
parameters.a_bar    = 0.4;          % maximum value of effort to change z
parameters.lambda   = 0.020;        % poisson intensity parameter
parameters.omega    = 0.010;        % green preference parameter of principal           
parameters.tau      = 0.30;         % tax on brown/ subsidy on green
parameters.xi       = 0.00;         % green preference parameter of agent (on z)
parameters.pol      = 15;           % polinomial length
parameters.phi      = 0.00;         % manager's preference for doing green effort (on a)

% model types
types = ["firstBest","agency"];
modelSolutions = {};

% solve model, derivatives and cross derivatives for effort
parfor i=1:length(types)
    model = transitionModel( types(i), parameters );
    [sol_A, sol_B] = model.solve;
    modelSolutions{i} = model;
end

% write results to structure
results = struct;
results.solutions = modelSolutions;

%% create tile plots
t=tiledlayout(2,2);

% Tile 1
nexttile
plot(modelSolutions{1}.sol_BShock.x, modelSolutions{1}.sol_BShock.y(1,:), 'b')
hold on
plot(modelSolutions{2}.sol_BShock.x, modelSolutions{2}.sol_BShock.y(1,:), 'r')
title('Firm Value Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 2
nexttile
plot(modelSolutions{1}.sol_AShock.x, modelSolutions{1}.sol_AShock.y(1,:), 'b')
hold on
plot(modelSolutions{2}.sol_AShock.x, modelSolutions{2}.sol_AShock.y(1,:), 'r')
title('Firm Value After Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 3
nexttile
plot(modelSolutions{1}.sol_BShock.x, modelSolutions{1}.effort.aBeforeShock, 'b')
hold on
plot(modelSolutions{2}.sol_BShock.x, modelSolutions{2}.effort.aBeforeShock, 'r')
title('Effort Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 4
nexttile
plot(modelSolutions{1}.sol_AShock.x, modelSolutions{1}.effort.aAfterShock, 'b')
hold on
plot(modelSolutions{2}.sol_AShock.x, modelSolutions{2}.effort.aAfterShock, 'r')
title('Effort After Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

leg = legend('First Best', 'Moral Hazard');
leg.Layout.Tile = 'south';

toc


%%
% Tile 4
load('results_WithXi.mat')
nexttile
plot(modelSolutions{1}.sol_AShock.x, modelSolutions{1}.effort.aAfterShock, 'b')
hold on
plot(modelSolutions{2}.sol_AShock.x, modelSolutions{2}.effort.aAfterShock, 'g')
title('(d) Effort After Shock', 'interpreter','latex')
load('results_noPhi_noXi.mat')
plot(modelSolutions{2}.sol_AShock.x, modelSolutions{2}.effort.aAfterShock, 'r')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

%% cross derivative test
% parameters
parameters.r        = 0.05;         % discounting
parameters.sigma    = 0.80;         % volatility
parameters.mu_G     = 0.20;         % productivity of green capital
parameters.mu_B     = 0.30;         % productivity of brown capital
parameters.gamma    = 15;           % risk aversion coefficient
parameters.theta    = 7;            % constant scaling adjustment cost
parameters.a_bar    = 0.4;          % maximum value of effort to change z
parameters.lambda   = 0.025;        % poisson intensity parameter
parameters.omega    = 0.010;        % green preference parameter of principal           
parameters.tau      = 0.25;         % tax on brown/ subsidy on green
parameters.xi       = 0.00;         % green preference parameter of agent (on z)
parameters.pol      = 15;           % polinomial length
parameters.phi      = 0.00;         % manager's preference for doing green effort (on a)



% compute cross derivative of effort wrt parameters j (not omega) and omega
model = transitionModel( "firstBest", parameters );
var1 = "lambda"; relativeStep1 = 0.25;
var2 = "omega"; relativeStep2 = 0.25;

[crossDerivativeAfterShock, crossDerivativeBeforeShock] = ...
model.crossDerivative(var1,var2, relativeStep1,relativeStep2);

figure;
plot(linspace(0,1,10000),crossDerivativeBeforeShock)
hold on
plot(linspace(0,1,10000),crossDerivativeAfterShock)
title('Cross derivative a wrt lambda omega')
legend('before', 'after')


%% Compute second derivative of effort and parameters

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

types = ["firstBest", "agency"];

% parameters for derivatives
paramsDerivatives = ["lambda", "tau", "omega", "sigma"];
L = length(paramsDerivatives);

% solution cell
modelDerivatives = {};

parfor i = 1:length(types)
    for j = 1:L
        model = transitionModel(types(i), parameters);

        var1 = paramsDerivatives(j);

        if paramsDerivatives(j) == "sigma"
            relativeStep1 = 0.0025;
        else
            relativeStep1 = 0.025;
        end

        [secondDerivativeAfterShock, secondDerivativeBeforeShock] = ...
        model.secondDerivative(var1, relativeStep1);
        modelDerivatives{i,j} = [secondDerivativeAfterShock, secondDerivativeBeforeShock]
    end
end



%%
x = linspace(0,1,10000);
t=tiledlayout(2,2);

nexttile
plot(x,modelDerivatives{1,1}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,1}(:,2), 'r')
title('$\frac{\partial^2 a}{\partial \lambda^2}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,2}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,2}(:,2), 'r')
title('$\frac{\partial^2 a}{\partial \tau^2}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,3}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,3}(:,2), 'r')
title('$\frac{\partial^2 a}{\partial \omega^2}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,4}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,4}(:,2), 'r')
title('$\frac{\partial^2 a}{\partial \sigma^2}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")


leg = legend('First Best', 'Moral Hazard');
leg.Layout.Tile = 'south';

exportgraphics(t,[pwd '/results/noPhiXi/tileSecondDerivatives.eps'], 'ContentType', 'vector')
%% Compute firm value derivatives with parameters
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
types = ["firstBest", "agency"];

% parameters for derivatives
paramsDerivatives = ["lambda", "tau", "omega", "sigma"];
derivativeType = "firmValue";
L = length(paramsDerivatives);

% solution cells
modelSolutions = {};
modelDerivatives = {};
modelCrossDerivatives = {};

% solve model, derivatives and cross derivatives for effort
parfor i=1:length(types)
    model = transitionModel( types(i), parameters );
    [sol_A, sol_B] = model.solve;
    modelSolutions{i} = model;
    
    for j=1:L        
        % compute single derivative of effort wrt parameter j
        model = transitionModel( types(i), parameters );
        var1 = paramsDerivatives(j); 
        
        % sigma is rather large, so need smaller step size
        if paramsDerivatives(j) == "sigma"
            relativeStep1 = 0.0025;
        else
            relativeStep1 = 0.025;
        end

        [derivativeAfterShock, derivativeBeforeShock] = ...
        model.derivative(var1, relativeStep1, derivativeType);
        
        % write single derivative to solution cell
        modelDerivatives{i,j} = [derivativeAfterShock, derivativeBeforeShock];
    end
end

toc
%%

figure;
t=tiledlayout(2,2);

nexttile
plot(x,modelDerivatives{1,2}(:,1), 'b')
hold on
plot(x,modelDerivatives{2,2}(:,1), 'r')
title('(b) $\frac{\partial V_1(z)}{\partial \tau}$', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,3}(:,1), 'b')
hold on
plot(x,modelDerivatives{2,3}(:,1), 'r')
title('(c) $\frac{\partial V_1(z)}{\partial \omega}$', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile([1 2])
plot(x,modelDerivatives{1,4}(:,1), 'b')
hold on
plot(x,modelDerivatives{2,4}(:,1), 'r')
title('(d) $\frac{\partial V_1(z)}{\partial \sigma}$', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")


leg = legend('First Best', 'Moral Hazard');
leg.Layout.Tile = 'south';
exportgraphics(t,[pwd '/results/noPhiXi/tileFirmValueDerivativesAfterShock.eps'], 'ContentType', 'vector')


figure;
t=tiledlayout(2,2);

nexttile
plot(x,modelDerivatives{1,1}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,1}(:,2), 'r')
title('(a) $\frac{\partial V_0(z)}{\partial \lambda}$', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,2}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,2}(:,2), 'r')
title('(b) $\frac{\partial V_0(z)}{\partial \tau}$', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,3}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,3}(:,2), 'r')
title('(c) $\frac{\partial V_0(z)}{\partial \omega}$', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,4}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,4}(:,2), 'r')
title('(d) $\frac{\partial V_0(z)}{\partial \sigma}$', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")


leg = legend('First Best', 'Moral Hazard');
leg.Layout.Tile = 'south';
exportgraphics(t,[pwd '/results/noPhiXi/tileFirmValueDerivativesBeforeShock.eps'], 'ContentType', 'vector')

%% Solve First Best for Different gammas
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
gammas = [2, 5, 10];

% solution cells
modelSolutions = {};

% solve model, derivatives and cross derivatives for effort
parfor i=1:length(gammas)
    model = transitionModel( "agency", parameters );
    model.parameters.gamma = gammas(i);
    [sol_A, sol_B] = model.solve;
    modelSolutions{i} = model;
end

% write results to structure
results = struct;
results.solutions = modelSolutions;

% save workspace
save('differentGammasFirstBest.mat')

toc


%% firm value + effort
% create tile plots
load('differentGammasFirstBest.mat')
t=tiledlayout(2,2);

% Tile 1
nexttile
plot(modelSolutions{1}.sol_BShock.x, modelSolutions{1}.sol_BShock.y(1,:))
hold on
plot(modelSolutions{2}.sol_BShock.x, modelSolutions{2}.sol_BShock.y(1,:))
plot(modelSolutions{3}.sol_BShock.x, modelSolutions{3}.sol_BShock.y(1,:))
title('(a) Firm Value Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 2
nexttile
plot(modelSolutions{1}.sol_AShock.x, modelSolutions{1}.sol_AShock.y(1,:))
hold on
plot(modelSolutions{2}.sol_AShock.x, modelSolutions{2}.sol_AShock.y(1,:))
plot(modelSolutions{3}.sol_AShock.x, modelSolutions{3}.sol_AShock.y(1,:))
title('(b) Firm Value After Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 3
nexttile
plot(modelSolutions{1}.sol_BShock.x, modelSolutions{1}.effort.aBeforeShock)
hold on
plot(modelSolutions{2}.sol_BShock.x, modelSolutions{2}.effort.aBeforeShock)
plot(modelSolutions{3}.sol_BShock.x, modelSolutions{3}.effort.aBeforeShock)
title('(c) Effort Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 4
nexttile
plot(modelSolutions{1}.sol_AShock.x, modelSolutions{1}.effort.aAfterShock)
hold on
plot(modelSolutions{2}.sol_AShock.x, modelSolutions{2}.effort.aAfterShock)
plot(modelSolutions{3}.sol_AShock.x, modelSolutions{3}.effort.aAfterShock)
title('(d) Effort After Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

leg = legend('$\gamma=2$','$\gamma=5$', '$\gamma=10$','interpreter', 'latex');
leg.Layout.Tile = 'south';

exportgraphics(t,[pwd '/results/noPhiXi/tilesDifferentGammas.eps'], 'ContentType', 'vector')

%% Solve First Best for Different omegas
tic

% parameters
parameters.r        = 0.05;         % discounting
parameters.sigma    = 1.00;         % volatility
parameters.mu_G     = 0.20;         % productivity of green capital
parameters.mu_B     = 0.30;         % productivity of brown capital
parameters.gamma    = 5;            % risk aversion coefficient
parameters.theta    = 10;           % constant scaling adjustment cost
parameters.a_bar    = 0.4;          % maximum value of effort to change z
parameters.lambda   = 0.03;        % poisson intensity parameter
parameters.omega    = 0.005;        % green preference parameter of principal           
parameters.tau      = 0.25;         % tax on brown/ subsidy on green
parameters.xi       = 0.00;         % green preference parameter of agent (on z)
parameters.pol      = 15;           % polinomial length
parameters.phi      = 0.00;         % manager's preference for doing green effort (on a)

% model types
omegas = [0.0, 0.005, 0.01];

% solution cells
modelSolutions = {};

% solve model, derivatives and cross derivatives for effort
parfor i=1:length(omegas)
    model = transitionModel( "firstBest", parameters );
    model.parameters.omega = omegas(i);
    [sol_A, sol_B] = model.solve;
    modelSolutions{i} = model;
end

% write results to structure
results = struct;
results.solutions = modelSolutions;

% save workspace
save('differentOmegasFirstBest.mat')

toc


%% firm value + effort
% create tile plots
load('differentOmegasFirstBest.mat')
t=tiledlayout(2,2);

% Tile 1
nexttile
plot(modelSolutions{1}.sol_BShock.x, modelSolutions{1}.sol_BShock.y(1,:))
hold on
plot(modelSolutions{2}.sol_BShock.x, modelSolutions{2}.sol_BShock.y(1,:))
plot(modelSolutions{3}.sol_BShock.x, modelSolutions{3}.sol_BShock.y(1,:))
title('(a) Firm Value Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 2
nexttile
plot(modelSolutions{1}.sol_AShock.x, modelSolutions{1}.sol_AShock.y(1,:))
hold on
plot(modelSolutions{2}.sol_AShock.x, modelSolutions{2}.sol_AShock.y(1,:))
plot(modelSolutions{3}.sol_AShock.x, modelSolutions{3}.sol_AShock.y(1,:))
title('(b) Firm Value After Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 3
nexttile
plot(modelSolutions{1}.sol_BShock.x, modelSolutions{1}.effort.aBeforeShock)
hold on
plot(modelSolutions{2}.sol_BShock.x, modelSolutions{2}.effort.aBeforeShock)
plot(modelSolutions{3}.sol_BShock.x, modelSolutions{3}.effort.aBeforeShock)
title('(c) Effort Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

% Tile 4
nexttile
plot(modelSolutions{1}.sol_AShock.x, modelSolutions{1}.effort.aAfterShock)
hold on
plot(modelSolutions{2}.sol_AShock.x, modelSolutions{2}.effort.aAfterShock)
plot(modelSolutions{3}.sol_AShock.x, modelSolutions{3}.effort.aAfterShock)
title('(d) Effort After Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

leg = legend('$\omega=0.0$','$\omega=0.005$', '$\omega=0.01$','interpreter', 'latex');
leg.Layout.Tile = 'south';

exportgraphics(t,[pwd '/results/differentOmegas/tilesDifferentOmegas.eps'], 'ContentType', 'vector')


%% Environmental Transition Model Solver

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
types = ["firstBest", "agency"];

% parameters for derivatives
paramsDerivatives = ["gamma", "tau", "omega", "sigma"];
derivativeType = "effort";
L = length(paramsDerivatives);

% solution cells
modelSolutions = {};
modelDerivatives = {};
modelCrossDerivatives = {};

% solve model, derivatives and cross derivatives for effort
parfor i=1:length(types)
    model = transitionModel( types(i), parameters );
    [sol_A, sol_B] = model.solve;
    modelSolutions{i} = model;
    
    for j=1:L        
        % compute single derivative of effort wrt parameter j
        model = transitionModel( types(i), parameters );
        var1 = paramsDerivatives(j); 
        
        % gamma is rather large, so need smaller step size
        if paramsDerivatives(j) == "sigma"
            relativeStep1 = 0.0025;
        else
            relativeStep1 = 0.025;
        end

        [derivativeAfterShock, derivativeBeforeShock] = ...
        model.derivative(var1, relativeStep1, derivativeType);
        
        % write single derivative to solution cell
        modelDerivatives{i,j} = [derivativeAfterShock, derivativeBeforeShock];

        if paramsDerivatives(j) ~= "lambda"
            % compute cross derivative of effort wrt parameters j (not omega) and omega
            model = transitionModel( types(i), parameters );
            var1 = paramsDerivatives(j); 
            
            % gamma is rather large, so need smaller step size
            if paramsDerivatives(j) == "sigma"
                relativeStep1 = 0.0025;
            elseif paramsDerivatives(j) == "theta"
                relativeStep1 = 0.0025
            else
                relativeStep1 = 0.025;
            end
            
            var2 = "lambda"; relativeStep2 = 0.025;

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

% save workspace
save('results_cross_lambda.mat')

toc

%% derivatives: paramsDerivatives = ["lambda", "tau","omega", "gamma"];
load('results_cross_lambda.mat')
x = linspace(0,1,10000);

t=tiledlayout(2,2);

nexttile
plot(x,modelDerivatives{1,1}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,1}(:,2), 'r')
title('(a) $\frac{\partial a}{\partial \gamma}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,2}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,2}(:,2), 'r')
title('(b) $\frac{\partial a}{\partial \tau}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,3}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,3}(:,2), 'r')
title('(c) $\frac{\partial a}{\partial \omega}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelDerivatives{1,4}(:,2), 'b')
hold on
plot(x,modelDerivatives{2,4}(:,2), 'r')
title('(d) $\frac{\partial a}{\partial \sigma}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

leg = legend('First Best', 'Moral Hazard');
leg.Layout.Tile = 'south';
%saveas(gca,[pwd '/results/tileDerivatives_noXiPhi.eps'],'epsc')

exportgraphics(t,[pwd '/results/noPhiXi/tileDerivativesLambda.eps'], 'ContentType', 'vector')

%% cross derivatives

x = linspace(0,1,10000);

% cross derivatives: paramsDerivatives = ["lambda", "tau","omega", "gamma"];
t=tiledlayout(2,2);

nexttile
plot(x,modelCrossDerivatives{1,1}(:,2), 'b')
hold on
plot(x,modelCrossDerivatives{2,1}(:,2), 'r')
title('(a) $\frac{\partial^2 a}{\partial \gamma \partial \lambda}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile
plot(x,modelCrossDerivatives{1,2}(:,2), 'b')
hold on
plot(x,modelCrossDerivatives{2,2}(:,2), 'r')
title('(b) $\frac{\partial^2 a}{\partial \tau \partial \lambda}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

nexttile([1 2])
plot(x,modelCrossDerivatives{1,4}(:,2), 'b')
hold on
plot(x,modelCrossDerivatives{2,4}(:,2), 'r')
title('(c) $\frac{\partial^2 a}{\partial \sigma \partial \lambda}$ Before Shock', 'interpreter','latex')
xlabel('z')
grid on
ylim("padded")
xlim("tight")

leg = legend('First Best', 'Moral Hazard');
leg.Layout.Tile = 'south';
%saveas(gca,[pwd '/results/tileCrossDerivatives_noXiPhi.eps'],'epsc')

exportgraphics(t,[pwd '/results/noPhiXi/tileCrossDerivativesLambda.eps'], 'ContentType', 'vector')