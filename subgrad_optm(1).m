function [h_best, pi_opt, iter, okcom, newnl, best_primal] = subgrad_optm(dimX, dimY, k, com, max_iter, theta_init)
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
    %   best_primal  - Best primal objective value found

    % Initialize variables
    pi = zeros(dimX * dimY * 2, 1);  % Dual variables
    h_best = inf;  % Best dual objective value found (minimize h)
    theta = theta_init;  % Initial step size multiplier
    h_lbd = 0;
    s_k = 10; % for s^k rule
    best_primal = inf;

    % Initialize solution tracking variables
    okcom = zeros(k, 1);  % Track which commodities are routed
    newnl = [];   % Store routes for each commodity

    % To track progress over iterations
    h_values = zeros(max_iter, 1);  % Store h(pi) values for each iteration
    primal_values = zeros(max_iter, 1);
    
    % Start of ergodic sequence
    x_ergodic = zeros(2 * dimX * dimY, 2* dimX * dimY, k);
    x_old = zeros(2 * dimX * dimY, 2* dimX * dimY, k);
    best_x = zeros(2 * dimX * dimY, 2* dimX * dimY, k);

    % Precompute adjacency matrix for path finding
    adj_matrix = ones(2 * dimX * dimY, 2 * dimX * dimY);

    % Iterate for subgradient optimization
    for iter = 1:max_iter
        % Initialize x_{ijl} and tracker of feasible routes
        x_new = zeros(2 * dimX * dimY, 2* dimX * dimY, k);
        feasible_routes = cell(1, k);
        feasible_routes(:) = {0};

        % Step 1: Solve Lagrangean subproblem using gsp
        nl = gsp(dimX, dimY, pi, k, com);
        
        % Step 2: Calculate h(pi) and extract feasible routes
        h_pi = sum(pi);  % Start with sum of dual variables
        iter_okcom = zeros(k, 1);  % Initialize okcom for this iteration
        last = 0;        % Initialize route traversal
        
        for i = 1:k
            % Find route for contact pair
            first = last + 1;
            slask = find(nl(last+1:length(nl)) == com(i,1));
            if isempty(slask)
                continue;  % Skip if no route found
            end
            last = slask(1) + first - 1;
            route = nl(first:last);  % Extract route nodes
            
            % Calculate cost of the route
            cost = sum(pi(route));
            
            % If route is feasible, update h(pi) and store details
            if cost < 1
                % Add route to solution
                for idx = 1:length(route)-1
                    x_new(route(idx), route(idx+1), i) = 1;
                end
                h_pi = h_pi + (1 - cost);  % Reward feasibility
                iter_okcom(i) = 1;
            end
        end
        
        % Update okcom with current iteration's results
        okcom = iter_okcom;
        
        % Update best dual value
        if h_pi < h_best
            h_best = h_pi;
            x_ergodic = x_new;
        end

        %HEURISTIC AND ERGODIC SEQUENCE STUFF
        %Calculate s^k rule stuff
        if iter ~= 1
            % Calculate denominators and numerators for s^k-rule
            denom1 = 0;
            for s = 0:iter-2
                denom1 = denom1 + (s+1)^s_k;
            end

            num1 = 0;
            for s = 0:iter-1
                num1 = num1 + (s+1)^s_k;
            end
    
            % Update x_ergodic
            x_ergodic = (denom1 / num1) * x_ergodic + (iter^(s_k) / num1) * x_old;
        else
            % First iteration: Initialize x_ergodic to x_new (x_0)
            x_ergodic = x_new; %For first iteration
        end
        
        % Step 3: Compute subgradients
        gamma = zeros(size(pi));
        for i = 1:length(pi)
            gamma(i) = 1 - sum(nl == i);
        end
        
        % Step 4: Update dual variables
        step = theta * (h_pi - h_lbd) / (norm(gamma)^2 + 1e-10);
        pi = max(0, pi - step * gamma);
        
        % Step 5: Update parameters
        if mod(iter, 5) == 0
            theta = theta * 0.95;
        end
        
        % Store values for plotting
        h_values(iter) = h_pi;
        primal_values(iter) = best_primal;
        x_old = x_new;
    end
    
    % Return results
    pi_opt = pi;  % Optimized dual variables
end
