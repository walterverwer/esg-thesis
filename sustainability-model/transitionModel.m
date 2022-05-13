classdef transitionModel < handle
    %TRANSITIONMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type % 'firstBest', 'agency'
        parameters
%         parameters = struct( ...
%             'r',        0.05, ...       discounting
%             'sigma',    0.25, ...       volatility
%             'mu_G',     0.20, ...       productivity of green capital
%             'mu_B',     0.30, ...       productivity of brown capital
%             'gamma',    5, ...          risk aversion coefficient
%             'theta',    10, ...         constant scaling adjustment cost
%             'a_bar',    0.3, ...        maximum value of effort to change z
%             'lambda',   0.02, ...       Poisson intensity parameter
%             'omega',    0.02, ...       green preference parameter of principal
%             'tau',      0.30, ...       tax on brown/ subsidy on green
%             'xi',       0.02, ...       green preference parameter of agent
%             'pol',      17, ...         polinomial length
%             'phi'       0.05)          % manager's preference for doing green effort (on a)
       
        boundaryValues
        % ^--- boundaryValues.p_BA
        %      boundaryValues.p_GA

        taylorExpression
        sol_AShock
        sol_BShock
        AShock_fun

        effort
        % ^--- effort.effortAShock
        %      effort.effortBShock
    end

    methods
        function obj = transitionModel(type,parameters)
            obj.type = type;
            if ~isempty(parameters)
                obj.parameters = parameters;
            end
        end
        
        function [r, sigma, mu_G, mu_B, gamma, theta, a_bar, lambda, omega, tau, xi, pol, phi] = getParams(obj)
            % getParams extracts the parameters from obj.parameters.
            r =         obj.parameters.r;
            sigma =     obj.parameters.sigma;
            mu_G =      obj.parameters.mu_G;
            mu_B =      obj.parameters.mu_B;
            gamma =     obj.parameters.gamma;
            theta =     obj.parameters.theta;
            a_bar =     obj.parameters.a_bar;
            lambda =    obj.parameters.lambda;
            omega =     obj.parameters.omega;
            tau =       obj.parameters.tau;
            xi =        obj.parameters.xi;
            pol =       obj.parameters.pol;
            phi =       obj.parameters.phi;
        end


        function residuals = bc(~, ya, yb, p_B, p_G)
            % Write boundary conditions in the form g(ya, yb) = 0, where ya denotes
            % V(0) and yb denotes V(1).
            % get states
                % z = 0
                Va  = ya(1);
                dVa = ya(2);
    
                % z = 1
                Vb  = yb(1);
                dVb = yb(2);

            residuals = [   Va - p_B       % @ z = 0, p_BA = Va
                            Vb - p_G ];    % @ z = 1, p_GA = Vb
        end


        function y = guessAShock(obj, z, y0_hat)
            % when z=0, then y = y0_hat
            if z==0
                y=y0_hat;
                return
            end

            % nummerical integration options
            options = odeset( 'RelTol', 1e-4, ...
                              'AbsTol', 1e-4 );
            
            % solve numerical integration
            sol = ode89( @(z,y) obj.odeAShock(z,y), [0 z], y0_hat, options );
            y = sol.y( :, end );
        end


        function y = guessBShock(obj, z, y0_hat)
            % when z=0, then y = y0_hat
            if z==0
                y=y0_hat;
                return
            end

            % nummerical integration options
            options = odeset( 'RelTol', 1e-4, ...
                              'AbsTol', 1e-4 );
            
            % solve numerical integration
            sol = ode89( @(z,y) obj.odeBShock(z,y), [0 z], y0_hat, options );
            y = sol.y( :, end );
        end


        function dy = odeAShock(obj, z, y)
            %   ode function to solve for j_s(z)
            %   write the second order ODE as a system of first-order equations. Define
            %   y(1) := j(z) and y(2) := j'(z)

            % get parameters
            [r, sigma, mu_G, mu_B, gamma, theta, a_bar, lambda, omega, tau, xi, pol, phi] = obj.getParams;
            
            % get taylor expression
            expr = obj.taylorExpression;

            % get current state
            V   = y(1);
            dV  = y(2);
            
            % ode for first best case
            if obj.type=="firstBest"
                % calculate optimal control effort
                a = dV / theta;
    
                % bound a to [-a_bar, a_bar]
                if abs( a ) > a_bar
                    a = sign( a ) * a_bar;
                end
                
                % ode after a shock
            
                ddV = expr( z ) * ( ...
                        r * V + theta * ( a^2 * z * (1-z) ) / 2 ...
                        - mu_G * z - ( mu_B - tau ) * (1-z) - omega * z - dV * a * z * (1-z) ...
                        );
            
            % ode for agency friction case
            elseif obj.type=="agency"

                % calculate optimal control effort
                a = 1 / ( theta + theta^2 * gamma * r * sigma^2 * z * (1-z) ) * ( dV + phi + theta * gamma * r * sigma^2 * z * (1-z) );

                % bound a to [-a_bar, a_bar]
                if abs( a ) > a_bar
                    a = sign( a ) * a_bar;
                end
                
                % ode after a shock

                ddV = expr( z ) * ( ...
                        r * V + theta * ( a^2 * z * (1-z) )/2 - phi * a * z * (1-z) ...
                        - mu_G * z - ( mu_B - tau ) * (1-z) - omega * z - xi * z ...
                        - dV *a * z * (1-z) ...
                        + ( gamma*r ) / 2 * ( (theta * a - phi)^2 * z^2 * (1-z)^2 *sigma^2 ) ...
                        );
            end
            
            % write to system of first order differential equations
            dy = [ dV; ddV ];
        
        end


        function dy = odeBShock(obj, z, y)
            %   ode function to solve for V(z)
            %   write the second order ODE as a system of first-order equations. Define

            % get parameters
            [r, sigma, mu_G, mu_B, gamma, theta, a_bar, lambda, omega, tau, xi, pol, phi] = obj.getParams;
            
            % get taylor expression
            expr = obj.taylorExpression;

            % get current state
            V   = y(1);
            dV  = y(2);
            
            % ode for first best case
            if obj.type=="firstBest"

                % calculate optimal control effort
                a = dV / theta;
    
                % bound a to [-a_bar, a_bar]
                if abs( a ) > a_bar
                    a = sign( a ) * a_bar;
                end
                
                % ode before a shock
    
                ddV = expr(z) * ( ...
                        ( r + lambda ) * V + theta * ( a^2 * z * (1-z) )/2 ...
                        - mu_G * z - mu_B * (1-z) - omega * z  ...
                        - lambda * obj.AShock_fun(z) - dV * a * z * (1-z) ...
                        );
            
            % ode for agency friction case
            elseif obj.type=="agency"

                a = 1 / ( theta + theta^2 * gamma * r * sigma^2 * z * (1-z) ) * ( dV + phi + theta * gamma * r * sigma^2 * z * (1-z) );


                % bound a to [-a_bar, a_bar]
                if abs( a ) > a_bar
                    a = sign( a ) * a_bar;
                end
                
                % ode before a shock

                ddV = expr(z) * ( ...
                        ( r + lambda ) * V + theta * ( a^2*z*(1-z) ) / 2 - phi * a * z * (1-z) ...
                        - mu_G * z - mu_B * (1-z) - omega * z - xi * z...
                        - lambda * obj.AShock_fun(z) - dV * a * z * (1-z) ...
                        + ( gamma * r ) / 2 * ( (theta * a - phi)^2 * z^2 * (1-z)^2 * sigma^2) ...
                        );
            end

            % write to system of first order differential equations
            dy = [ dV; ddV ];
            
        end


        function [ solAfterShock, solBeforeShock ] = solve(obj)
                % get parameters
                [r, sigma, mu_G, mu_B, gamma, theta, a_bar, lambda, omega, tau, xi, pol, phi] = obj.getParams;
                
                % get boundary values after shock
                if obj.type == "firstBest"
                    p_BA = ( mu_B - tau ) / r;
                    p_GA = ( mu_G + omega ) / r;

                elseif obj.type == "agency"
                    p_BA = ( mu_B - tau ) / r;
                    p_GA = ( mu_G + omega + xi ) / r;

                end

                % store boundary values to properties
                obj.boundaryValues.p_BA = p_BA;
                obj.boundaryValues.p_GA = p_GA;
                
                % create taylor expression for z in {0, 1} for
                % singularities
                % define a symbolic variable for z --> x
                syms x
                obj.taylorExpression = matlabFunction(...
                        taylor(... 
                            2 / ( sigma^2 * x^2 * ( 1 - x )^2 ), ...
                            x, ...
                            'ExpansionPoint',   0.5, ...
                            'Order',    9 ...
                        ) ...
                    );

                % create mesh for boundary value solver before and after
                % shock
                xmesh = linspace(0, 1, 5000);

                % initial guess
                y0_hat = [ p_BA; 0 ];
                
                % get solution for initial guess
                solinit = bvpinit(xmesh, @(z) obj.guessAShock(z,y0_hat));
                
                bvpoptions = bvpset( Stats="on", Nmax=150000, AbsTol=1e-4, RelTol=1e-4);
                solAfterShock = bvp5c(...
                                @(z,y) obj.odeAShock(z,y),...
                                @(ya, yb) obj.bc(ya, yb, p_BA, p_GA),...
                                solinit,...
                                bvpoptions);
                obj.sol_AShock = solAfterShock;
                
                % Fit value function after shock
                fit_AShock = polyfit(solAfterShock.x, solAfterShock.y(1,:),pol);
                
                syms zi
                AShock_funFit = fit_AShock*zi.^(pol:-1:0)';
                obj.AShock_fun = matlabFunction(AShock_funFit);


                % solve the model before the shock, using the polyfit after the shock

                % obtain boundary values before the shock
                if obj.type == "firstBest"
                    p_BB = ( mu_B + lambda * p_BA ) / ( lambda + r );
                    p_GB = ( mu_G + omega + lambda * p_GA ) / ( r + lambda );

                elseif obj.type == "agency"
                    p_BB = ( mu_B + lambda * p_BA ) / ( lambda + r );
                    p_GB = ( mu_G + omega + xi + lambda * p_GA ) / ( r + lambda );

                end

                % store boundary values to properties
                obj.boundaryValues.p_BB = p_BB;
                obj.boundaryValues.p_GB = p_GB;
                
                % initial guess
                y0_hat = [p_BB;0]; % start is known, V' is unknown. If code gives error, vary the derivative.
                solinit = bvpinit(xmesh, @(z) obj.guessBShock(z,y0_hat));
                
                bvpoptions = bvpset(Stats="on",Nmax=150000,AbsTol=1e-4,RelTol=1e-4);
                solBeforeShock = bvp5c( @(z,y) obj.odeBShock(z,y), ...
                                    @(ya,yb) obj.bc(ya, yb, p_BB, p_GB), ...
                                    solinit, ...
                                    bvpoptions);
                obj.sol_BShock = solBeforeShock;

                % store effort before and after shock
                [aAfterShock, aBeforeShock] = getEffort(obj);
                obj.effort.aAfterShock = aAfterShock;
                obj.effort.aBeforeShock = aBeforeShock;
        end


        function [aAfterShock, aBeforeShock] = getEffort(obj)
            % get parameters
            [r, sigma, mu_G, mu_B, gamma, theta, a_bar, lambda, omega, tau, xi, pol, phi] = obj.getParams;

            % effort for first best
            if obj.type == "firstBest"
                aAfterShock = obj.sol_AShock.y(2,:) / theta;
                aBeforeShock = obj.sol_BShock.y(2,:) / theta;
            
            % effort for agency
            elseif obj.type == "agency"
                aAfterShock = 1 ./ ( theta + theta^2 * gamma * r * sigma^2 .* obj.sol_AShock.x .* (1-obj.sol_AShock.x) ) .* ( obj.sol_AShock.y(2,:) + phi + theta * gamma * r * sigma^2 .* obj.sol_AShock.x .* (1-obj.sol_AShock.x) );

                aBeforeShock = 1 ./ ( theta + theta^2 * gamma * r * sigma^2 .* obj.sol_BShock.x .* (1-obj.sol_BShock.x) ) .* ( obj.sol_BShock.y(2,:) + phi + theta * gamma * r * sigma^2 .* obj.sol_BShock.x .* (1-obj.sol_BShock.x) );
            end

            
        end

        function [crossDerivativeAfterShock, crossDerivativeBeforeShock] = crossDerivative(obj, var1, var2, relativeStep1, relativeStep2)
            % get parameters
            [r, sigma, mu_G, mu_B, gamma, theta, a_bar, lambda, omega, tau, xi, pol, phi] = obj.getParams;
            
            % define stepsize
            switch var1
                case "r"
                    stepVar1 = r * relativeStep1;
                case "sigma"
                    stepVar1 = sigma * relativeStep1;
                case "mu_G"
                    stepVar1 = mu_G * relativeStep1;
                case "mu_B"
                    stepVar1 = mu_B * relativeStep1;
                case "gamma"
                    stepVar1 = gamma * relativeStep1;
                case "theta"
                    stepVar1 = theta * relativeStep1;
                case "a_bar"
                    stepVar1 = a_bar * relativeStep1;
                case "lambda"
                    stepVar1 = lambda * relativeStep1;
                case "omega"
                    stepVar1 = omega * relativeStep1;
                case "tau"
                    stepVar1 = tau * relativeStep1;
                case "xi"
                    stepVar1 = xi * relativeStep1;
                case "phi"
                    stepVar1 = phi * relativeStep1;
            end

            switch var2
                case "r"
                    stepVar2 = r * relativeStep2;
                case "sigma"
                    stepVar2 = sigma * relativeStep2;
                case "mu_G"
                    stepVar2 = mu_G * relativeStep2;
                case "mu_B"
                    stepVar2 = mu_B * relativeStep2;
                case "gamma"
                    stepVar2 = gamma * relativeStep2;
                case "theta"
                    stepVar2 = theta * relativeStep2;
                case "a_bar"
                    stepVar2 = a_bar * relativeStep2;
                case "lambda"
                    stepVar2 = lambda * relativeStep2;
                case "omega"
                    stepVar2 = omega * relativeStep2;
                case "tau"
                    stepVar2 = tau * relativeStep2;
                case "xi"
                    stepVar2 = xi * relativeStep2;
                case "phi"
                    stepVar2 = phi * relativeStep2;
            end

            situations = [  stepVar1, stepVar2;
                            stepVar1, -stepVar2;
                            -stepVar1, stepVar2;
                            -stepVar1, -stepVar2;];

            xq = linspace(0,1,10000);
            solutionsAfterShock = zeros(length(xq),size(situations,1));
            solutionsBeforeShock = zeros(length(xq),size(situations,1));

            for i=1:size(situations,1)
                switch var1
                    case "r"
                        obj.parameters.r = r + situations(i,1);
                    case "sigma"
                        obj.parameters.sigma = sigma + situations(i,1);
                    case "mu_G"
                        obj.parameters.mu_G = mu_G + situations(i,1);
                    case "mu_B"
                        obj.parameters.mu_B = mu_B + situations(i,1);
                    case "gamma"
                        obj.parameters.gamma = gamma + situations(i,1);
                    case "theta"
                        obj.parameters.theta = theta + situations(i,1);
                    case "a_bar"
                        obj.parameters.a_bar = a_bar + situations(i,1);
                    case "lambda"
                        obj.parameters.lambda = lambda + situations(i,1);
                    case "omega"
                        obj.parameters.omega = omega + situations(i,1);
                    case "tau"
                        obj.parameters.tau = tau + situations(i,1);
                    case "xi"
                        obj.parameters.xi = xi + situations(i,1);
                    case "phi"
                        obj.parameters.phi = phi + situations(i,1);
                end

                switch var2
                    case "r"
                        obj.parameters.r = r + situations(i,2);
                    case "sigma"
                        obj.parameters.sigma = sigma + situations(i,2);
                    case "mu_G"
                        obj.parameters.mu_G = mu_G + situations(i,2);
                    case "mu_B"
                        obj.parameters.mu_B = mu_B + situations(i,2);
                    case "gamma"
                        obj.parameters.gamma = gamma + situations(i,2);
                    case "theta"
                        obj.parameters.theta = theta + situations(i,2);
                    case "a_bar"
                        obj.parameters.a_bar = a_bar + situations(i,2);
                    case "lambda"
                        obj.parameters.lambda = lambda + situations(i,2);
                    case "omega"
                        obj.parameters.omega = omega + situations(i,2);
                    case "tau"
                        obj.parameters.tau = tau + situations(i,2);
                    case "xi"
                        obj.parameters.xi = xi + situations(i,2);
                    case "phi"
                        obj.parameters.phi = phi + situations(i,2);
                end
                [ solAfterShock, solBeforeShock ] = obj.solve;    

                % obtain effort
                [aAfterShock, aBeforeShock] = getEffort(obj);
                
                % apply makima to obtain fit of solution
                solutionsAfterShock(:,i) = makima(solAfterShock.x, aAfterShock, xq);
                solutionsBeforeShock(:,i) = makima(solBeforeShock.x, aBeforeShock, xq);

                obj.parameters.r = r;
                obj.parameters.sigma = sigma;
                obj.parameters.mu_G = mu_G;
                obj.parameters.mu_B = mu_B;
                obj.parameters.gamma = gamma;
                obj.parameters.theta = theta;
                obj.parameters.a_bar = a_bar;
                obj.parameters.lambda = lambda;
                obj.parameters.omega = omega;
                obj.parameters.tau = tau;
                obj.parameters.xi = xi;
                obj.parameters.phi = phi;
            end
            crossDerivativeAfterShock = ...
            (   solutionsAfterShock(:,1) ...
                - solutionsAfterShock(:,2) ...
                - solutionsAfterShock(:,3) ...
                + solutionsAfterShock(:,4)  ) ...
            / ( 4 * stepVar1 * stepVar2 );

            crossDerivativeBeforeShock = ...
            (   solutionsBeforeShock(:,1) ...
                - solutionsBeforeShock(:,2) ...
                - solutionsBeforeShock(:,3) ...
                + solutionsBeforeShock(:,4)  ) ...
            / ( 4 * stepVar1 * stepVar2 );
        end


        function [derivativeAfterShock, derivativeBeforeShock] = derivative(obj, var1, relativeStep1)
            % get parameters
            [r, sigma, mu_G, mu_B, gamma, theta, a_bar, lambda, omega, tau, xi, pol, phi] = obj.getParams;
            
            % define stepsize
            switch var1
                case "r"
                    stepVar1 = r * relativeStep1;
                case "sigma"
                    stepVar1 = sigma * relativeStep1;
                case "mu_G"
                    stepVar1 = mu_G * relativeStep1;
                case "mu_B"
                    stepVar1 = mu_B * relativeStep1;
                case "gamma"
                    stepVar1 = gamma * relativeStep1;
                case "theta"
                    stepVar1 = theta * relativeStep1;
                case "a_bar"
                    stepVar1 = a_bar * relativeStep1;
                case "lambda"
                    stepVar1 = lambda * relativeStep1;
                case "omega"
                    stepVar1 = omega * relativeStep1;
                case "tau"
                    stepVar1 = tau * relativeStep1;
                case "xi"
                    stepVar1 = xi * relativeStep1;
                case "phi"
                    stepVar1 = phi * relativeStep1;
            end

            situations = [  2*stepVar1;
                            stepVar1;
                            -stepVar1;
                            -2*stepVar1];

            xq = linspace(0,1,10000);
            solutionsAfterShock = zeros(length(xq),size(situations,1));
            solutionsBeforeShock = zeros(length(xq),size(situations,1));

            for i=1:size(situations,1)
                switch var1
                    case "r"
                        obj.parameters.r = r + situations(i,1);
                    case "sigma"
                        obj.parameters.sigma = sigma + situations(i,1);
                    case "mu_G"
                        obj.parameters.mu_G = mu_G + situations(i,1);
                    case "mu_B"
                        obj.parameters.mu_B = mu_B + situations(i,1);
                    case "gamma"
                        obj.parameters.gamma = gamma + situations(i,1);
                    case "theta"
                        obj.parameters.theta = theta + situations(i,1);
                    case "a_bar"
                        obj.parameters.a_bar = a_bar + situations(i,1);
                    case "lambda"
                        obj.parameters.lambda = lambda + situations(i,1);
                    case "omega"
                        obj.parameters.omega = omega + situations(i,1);
                    case "tau"
                        obj.parameters.tau = tau + situations(i,1);
                    case "xi"
                        obj.parameters.xi = xi + situations(i,1);
                    case "phi"
                        obj.parameters.phi = phi + situations(i,1);
                end

                [ solAfterShock, solBeforeShock ] = obj.solve;   

                % obtain effort
                [aAfterShock, aBeforeShock] = getEffort(obj);
                
                % apply makima to obtain fit of solution
                solutionsAfterShock(:,i) = makima(solAfterShock.x, aAfterShock, xq);
                solutionsBeforeShock(:,i) = makima(solBeforeShock.x, aBeforeShock, xq);

                obj.parameters.r = r;
                obj.parameters.sigma = sigma;
                obj.parameters.mu_G = mu_G;
                obj.parameters.mu_B = mu_B;
                obj.parameters.gamma = gamma;
                obj.parameters.theta = theta;
                obj.parameters.a_bar = a_bar;
                obj.parameters.lambda = lambda;
                obj.parameters.omega = omega;
                obj.parameters.tau = tau;
                obj.parameters.xi = xi;
                obj.parameters.phi = phi;
            end

            derivativeAfterShock = ...
            (   -solutionsAfterShock(:,1) ...
                + 8 * solutionsAfterShock(:,2) ...
                - 8 * solutionsAfterShock(:,3) ...
                + solutionsAfterShock(:,4)  ) ...
            / ( 12 * stepVar1 );

            derivativeBeforeShock = ...
            (   - solutionsBeforeShock(:,1) ...
                + 8 * solutionsBeforeShock(:,2) ...
                - 8 * solutionsBeforeShock(:,3) ...
                + solutionsBeforeShock(:,4)  ) ...
            / ( 12 * stepVar1  );

        end
    end
end

