function [h_best, pi_opt, iter, okcom, newnl, h_lbd, best_route] = subgrad_optm(dimX, dimY, k, com, max_iter, theta_init)
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

    
    
    % To track progress over iterations
    h_values = zeros(max_iter, 1);  % Store h(pi) values for each iteration
    z_values = zeros(max_iter, 1);

    best_route = [];
    bestnl = [];

    % Iterate for subgradient optimization
    for iter = 1:max_iter

        % Initialize tracker of feasible routes
        
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
        % Initialize a tracker for node usage
        dim = dimX * dimY * 2; % Total number of nodes in the graph (including logical nodes)
        node_usage = zeros(dim, 1); % Tracks how many routes use each node
        
        % Iterate through each route in feasible_routes
        for i = 1:length(feasible_routes)
            curRoute = feasible_routes{i};
        
            % Skip if the current route is empty or invalid
            if isempty(curRoute) || isequal(curRoute, 0)
                continue;
            end
        
            % Check for conflicts with previously used nodes
            conflict = false;
            for j = 1:length(curRoute)
                node = curRoute(j);
                if node_usage(node) > 0
                    conflict = true;
                    break;
                end
            end
        
            % If there is a conflict, prune the route
            if conflict
                %fprintf('Conflict detected. Pruning route for contact pair %d: [%s]\n', i, num2str(curRoute));
                feasible_routes{i} = []; % Mark route as invalid
            else
                % If no conflict, update node usage
                for j = 1:length(curRoute)
                    node = curRoute(j);
                    node_usage(node) = node_usage(node) + 1;
                end
            end
        end

        num_feasible_routes = 0;

        for i = 1:length(feasible_routes)
            if ~isempty(feasible_routes{i}) && ~isequal(feasible_routes{i},0)
                num_feasible_routes = num_feasible_routes + 1;
            end
        end

        z_values(iter) = num_feasible_routes;

        if num_feasible_routes > h_lbd
            h_lbd = num_feasible_routes;
            best_route = feasible_routes;
        end


        




        

    end

    
    
    % Plot convergence of h(pi) and z
    figure;
    hold on; % Allows multiple plots on the same figure
    plot(1:iter, h_values(1:iter), '-o', 'LineWidth', 0.5, 'MarkerSize', 3, 'Color', 'b'); % h(pi) in blue
    plot(1:iter, z_values(1:iter), 'o', 'MarkerSize', 1, 'Color', 'r'); % z in red
    hold off;
    
    % Add title, labels, legend, and grid
    title('Convergence of h(u) and z');
    xlabel('Iteration');
    ylabel('Objective Values');
    legend('h(u)', 'z', 'Location', 'best'); % Add legend
    grid on; % Enable grid
    
    % Output results
    pi_opt = pi;  % Optimized dual variables
    fprintf('Optimization completed in %d iterations.\n', iter);
    fprintf('Best dual objective value: %.4f\n', h_best);
end
