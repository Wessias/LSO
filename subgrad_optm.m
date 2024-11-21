function [h_best, pi_opt, iter, okcom, newnl] = subgrad_optm(dimX, dimY, k, com, max_iter, theta_init)
    % SUBGRAD_OPTIMIZATION - Solves the Lagrangean dual for the VLSI problem.
    %
    % INPUTS:
    %   dimX, dimY   - Dimensions of the grid
    %   k            - Number of contact pairs
    %   com          - Contact pairs [k x 2 matrix]
    %   max_iter     - Maximum number of iterations
    %   theta_init   - Initial step size multiplier
    %
    % OUTPUTS:
    %   h_best       - Best dual objective value found
    %   pi_opt       - Optimized dual variables
    %   iter         - Total iterations performed
    %   okcom        - Feasible node pairs
    %   newnl        - Feasible node routes

    % Initialize dual variables
    pi = zeros(dimX * dimY * 2, 1);  % Dual variables
    h_best = inf;  % Best dual objective value found (minimize h)
    theta = theta_init;  % Initial step size multiplier
    h_lbd = 0;
    s_k = 10 % for s^k rule

    
    
    % To track progress over iterations
    h_values = zeros(max_iter, 1);  % Store h(pi) values for each iteration
    
    % Start of ergodic sequence
    x_ergodic = zeros(2 * dimX * dimY, 2* dimX * dimY, k)
    x_ergodic_old = zeros(2 * dimX * dimY, 2* dimX * dimY, k)
    x_old = zeros(2 * dimX * dimY, 2* dimX * dimY, k); 

    % Iterate for subgradient optimization
    for iter = 1:max_iter

        % Initialize x_{ijl} and tracker of feasible routes
        x_new = zeros(2 * dimX * dimY, 2* dimX * dimY, k);
        feasible_routes = cell(1, k);
        feasible_routes(:) = {0};



        % Step 1: Solve Lagrangean subproblem using gsp
        nl = gsp(dimX, dimY, pi, k, com);
        %fprintf('Iteration %d: Solved Lagrangean subproblem.\n', iter);
        
        % Step 2: Calculate h(pi) and extract feasible routes
        h_pi = sum(pi);  % Start with sum of dual variables
        okcom = [];      % Indices of feasible contact pairs
        newnl = [];      % Feasible routes (node indices)
        last = 0;        % Initialize route traversal
        
        for i = 1:k
            % Find route for contact pair
            first = last + 1;
            slask = find(nl(last+1:length(nl)) == com(i,1));
            last = slask(1) + first - 1;
            route = nl(first:last);  % Extract route nodes
            
            % Calculate cost of the route
            cost = sum(pi(nl(first:last)));
            %fprintf('Contact pair %d: Cost = %f, Route = [%s]\n', i, cost, num2str(route'));
            
            % If route is feasible, update h(pi) and store details
            if cost < 1
                okcom = [okcom, i];         % Add feasible pair index
                newnl = [newnl; route];    % Add route to feasible routes
                h_pi = h_pi + (1 - cost);  % Reward feasibility
                feasible_routes{i} = route;
            end
        end
        
        % Update best dual value
        h_best = min(h_best, h_pi);  % Minimize dual objective
        h_values(iter) = h_pi;  % Store current h(pi) for plotting
        %fprintf('Iteration %d: h(pi) = %.4f, h_best = %.4f\n', iter, h_pi, h_best);
        
        % Step 3: Compute subgradients
        gamma = zeros(size(pi));  % Initialize subgradient
        for i = 1:length(pi)
            gamma(i) = 1 - sum(newnl == i);  % Count routes passing through node i
        end
        %fprintf('Subgradient norm at iteration %d: %.4e\n', iter, norm(gamma));
        
        % Step 4: Calculate step length
        epsilon = 1e-6;  % Safeguard for numerical stability
        step_length = max(1e-4, theta * (h_pi - h_lbd) / (norm(gamma)^2 + epsilon));
        %fprintf('Step length at iteration %d: %.4e\n', iter, step_length);
        
        % Step 5: Update dual variables
        pi = max(0, pi - step_length * gamma);  % Ensure pi >= 0 (projection)
        %fprintf('Dual variables (pi) at iteration %d:\n', iter);
        %disp(pi');  % Display dual variables as a row vector
        %fprintf('Updated dual variables (pi) at iteration %d: [%s]\n', iter, num2str(pi'));

        
        % Step 6: Update theta periodically
        if mod(iter, 10) == 0
            theta = theta * 0.95;  % Decay step size multiplier
        end
        
        %HEURISTIC AND ERGODIC SEQUENCE STUFF

        %Populate x based on routes and valid paths.
        for l = 1:k       
            if ~isequal(feasible_routes{l}, 0)
                route = feasible_routes{l};
                %Populate based on route for pair l
                for idx = 1:length(route)-1
                    i = route(idx);       % Start node
                    j = route(idx+1);     % End node
                    x_new(i, j, l) = 1;       % Set x_{ijl} = 1
                end
                x_new(com(l, 1), com(l, 2), l) = 1; % Set logical path from start node to termination node = 1
            end
        end

        %Calculate s^k rule stuff
        if iter - 1 ~= 0
            % Calculate denominators and numerators for s^k-rule
            denom1 = 0;
            for s = 0:iter-1-2
                denom1 = denom1 + (s+1)^s_k;
            end

            num1 = 0;
            for s = 0:iter-1-1
                num1 = num1 + (s+1)^s_k;
            end
    
            % Update x_ergodic
            x_ergodic = (denom1 / num1) * x_ergodic + (iter^(s_k) / num1) * x_old;
        else
            % First iteration: Initialize x_ergodic to x_new (x_0)
            x_ergodic = x_new; %For first iteration
        end

        % Update x_old for the next iteration
        x_old = x_new;
        

    end

    x_binary = x_ergodic > 0.9;
    
    % Plot convergence of h(pi)
    figure;
    plot(1:iter, h_values(1:iter), '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
    title('Convergence of h(pi)');
    xlabel('Iteration');
    ylabel('h(pi)');
    grid on;
    
    % Output results
    pi_opt = pi;  % Optimized dual variables
    fprintf('Optimization completed in %d iterations.\n', iter);
    fprintf('Best dual objective value: %.4f\n', h_best);
end
