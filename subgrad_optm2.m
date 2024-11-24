function varargout = subgrad_optm4(dimX, dimY, k, com, max_iter, theta_init)
    % 
    % === How This Algorithm Works ===
    % 1. Start by allowing wires to overlap (relaxing the no-sharing rule)
    % 2. Add penalties (dual variables) for sharing nodes
    % 3. Repeatedly:
    %    - Find shortest paths considering these penalties
    %    - Update penalties based on how many wires share nodes
    %    - Every 100 iterations, try to find a real solution where no wires overlap
    %
    % === The Grid Structure ===
    % - We have two identical grids of size dimX × dimY (like two floors)
    % - First layer: nodes 1 to dimX*dimY
    % - Second layer: nodes (dimX*dimY + 1) to 2*dimX*dimY
    % - Wires can move within each layer or switch between layers
    %
    % === Inputs ===
    % dimX, dimY : Size of each grid layer (width and height)
    % k          : Number of pairs of points to connect (called commodities)
    % com        : k×2 matrix where each row is [start_point, end_point] to connect
    % max_iter   : Maximum number of times to try improving the solution
    % theta_init : Initial step size for updating penalties (like how aggressive we are)
    %
    % === Outputs ===
    % h_best      : Best lower bound on solution cost we found
    % pi          : Final penalties for using each node
    % iter        : How many iterations we actually ran
    % okcom       : Which pairs of points we successfully connected
    % nl          : The paths we found
    % best_primal : Total length of all wires in our best solution

    % === Initialize Variables ===
    % These variables keep track of our solution and help us improve it
    pi = zeros(dimX * dimY * 2, 1);     % Dual variables (Lagrange multipliers)
    h_best = inf;                        % Best dual bound
    theta = theta_init;                  % Step size multiplier
    h_lbd = 0;                          % Lower bound for step size calculation
    best_primal = inf;                  % Best primal objective value found
    okcom = zeros(k, 1);                % Tracks which commodities are routed

    % Initialize arrays for solution tracking
    h_values = zeros(max_iter, 1);       % History of dual values
    primal_values = zeros(max_iter, 1);  % History of primal values
    x_ergodic = zeros(2 * dimX * dimY, 2* dimX * dimY, k);  % Ergodic average of solutions
    x_old = zeros(2 * dimX * dimY, 2* dimX * dimY, k);      % Previous iteration's solution
    best_x = zeros(2 * dimX * dimY, 2* dimX * dimY, k);     % Best feasible solution found

    % Precompute adjacency matrix for path finding
    % This matrix will be modified during primal heuristic to avoid node conflicts
    adj_matrix = ones(2 * dimX * dimY, 2 * dimX * dimY);

    % === Main Optimization Loop ===
    % Each iteration tries to improve our solution
    for iter = 1:max_iter
        % === Step 1: Find Shortest Paths ===
        % Try to connect each pair of points, considering penalties
        % Initialize solution matrix for current iteration
        % x_new(i,j,l) = 1 if commodity l uses edge (i,j)
        x_new = zeros(2 * dimX * dimY, 2* dimX * dimY, k);

        % Solve Lagrangean Subproblem
        % gsp finds shortest paths for each commodity considering dual variables as costs
        nl = gsp(dimX, dimY, pi, k, com);

        % === Step 2: Calculate Solution Quality ===
        % Figure out how good our current solution is
        h_pi = sum(pi);  % Start with sum of all dual variables
        iter_okcom = zeros(k, 1);  % Track feasible commodities in this iteration
        last = 0;  % Index tracker for reading paths from nl

        % Process each commodity's path from gsp solution
        for i = 1:k
            first = last + 1;
            % Find next occurrence of starting point
            slask = find(nl(last+1:length(nl)) == com(i,1));
            if isempty(slask)
                continue;  % Skip if no path found
            end
            last = slask(1) + first - 1;
            route = nl(first:last);  % Extract path nodes

            % Calculate path cost using dual variables
            cost = sum(pi(route));
            
            % If path cost < 1, it's beneficial to use this path
            if cost < 1
                % Add path to solution
                for idx = 1:length(route)-1
                    x_new(route(idx), route(idx+1), i) = 1;
                end
                h_pi = h_pi + (1 - cost);  % Update dual objective
                iter_okcom(i) = 1;  % Mark commodity as routed
            end
        end

        % Update which commodities are currently routed
        okcom = iter_okcom;

        % Update best dual solution if improved
        if h_pi < h_best
            h_best = h_pi;
            % After calculating x_new and before the primal heuristic
s_k = iter;  % Using iteration number as s^k
if iter == 1
    x_ergodic = x_new;
else
    x_ergodic = ((s_k - 1) * x_ergodic + x_new) / s_k;
end
        end

        % === Step 3: Try to Find a Real Solution (Every 100 iterations) ===
        % Attempt to find paths that don't share any nodes
        if mod(iter, 100) == 0
            % Initialize variables for primal solution attempt
            current_x = zeros(2 * dimX * dimY, 2* dimX * dimY, k);  % Solution matrix
            node_usage = zeros(2 * dimX * dimY, 1);  % Track used nodes
            current_primal = 0;  % Current solution cost
            current_okcom = zeros(k, 1);  % Track routed commodities

            % Calculate priority for each commodity based on frequency in ergodic solution
            % More frequently used paths get higher priority
            priorities = zeros(k, 1);
            for l = 1:k
                priorities(l) = sum(sum(x_ergodic(:,:,l)));
            end
            [~, order] = sort(priorities, 'descend');  % Sort commodities by priority

            % Try to route each commodity in priority order
            for idx = 1:k
                l = order(idx);
                % Create weighted graph avoiding used nodes
                % Used nodes get higher weights to discourage their use
                weights = adj_matrix + 5 * node_usage;
                % Find path avoiding used nodes where possible
                path = fast_path(weights, com(l,1), com(l,2), dimX, dimY);

                if ~isempty(path)
                    % Check if path nodes are available (not used by other routes)
                    path_nodes = unique(path);
                    if all(node_usage(path_nodes) == 0)
                        % Add path to solution
                        path_length = length(path) - 1;
                        for i = 1:path_length
                            current_x(path(i), path(i+1), l) = 1;
                            node_usage(path(i)) = 1;  % Mark nodes as used
                        end
                        node_usage(path(end)) = 1;
                        current_primal = current_primal + path_length;  % Update solution cost
                        current_okcom(l) = 1;  % Mark commodity as routed
                    end
                end
            end

            % Update best primal solution if more commodities are routed
            if sum(current_okcom) > sum(okcom)
                best_primal = current_primal;
                best_x = current_x;
                okcom = current_okcom;
            end
        end

        % === Step 4: Update Penalties ===
        % Change penalties based on how many wires share nodes
        % Compute subgradients for dual update
        % Subgradient for each node = 1 - number of times node is used
        gamma = zeros(size(pi));
        for i = 1:length(pi)
            gamma(i) = 1 - sum(nl == i);
        end

        % Update dual variables using subgradient step
        step = theta * (h_pi - h_lbd) / (norm(gamma)^2 + 1e-10);
        pi = max(0, pi - step * gamma);  % Project to non-negative orthant

        % === Step 5: Update Parameters ===
        % Reduce step size periodically to help convergence
        if mod(iter, 5) == 0
            theta = theta * 0.95;
        end

        % Store values for analysis
        h_values(iter) = h_pi;
        primal_values(iter) = best_primal;
        x_old = x_new;
    end

    % Return results
    varargout{1} = h_best;
    varargout{2} = pi;
    varargout{3} = iter;
    varargout{4} = okcom;
    varargout{5} = nl;
    varargout{6} = best_primal;
end

function path = fast_path(weights, start, target, dimX, dimY)
    % FAST_PATH - Finds a good path between two points on our two-layer grid
    %
    % === How It Works ===
    % 1. Start at the beginning point
    % 2. Each step:
    %    - Look at all neighboring spots we haven't visited
    %    - Score each neighbor based on:
    %      * How "expensive" it is to go there (from weights)
    %      * How far it is from our target
    %      * How many times we need to switch floors
    %    - Move to the best-scoring neighbor
    % 3. Stop when we:
    %    - Reach the target
    %    - Can't find any unvisited neighbors
    %    - Take too many steps (to prevent infinite loops)
    %
    % === Inputs ===
    % weights : Matrix showing how "expensive" each move is
    % start   : Starting point number
    % target  : Ending point number
    % dimX    : Width of each layer
    % dimY    : Height of each layer
    %
    % === Output ===
    % path : List of points making up our route (empty if no path found)
    
    n = size(weights, 1);
    path = [start];
    current = start;
    visited = false(n, 1);
    visited(start) = true;
    max_steps = min(100, 8*max(dimX, dimY));  % Limit path length to prevent infinite loops
    steps = 0;
    grid_size = dimX * dimY;

    while current ~= target && steps < max_steps
        steps = steps + 1;
        neighbors = find(weights(current,:) > 0 & ~visited');
        if isempty(neighbors)
            % If we can't find any unvisited neighbors, give up and return empty path
            path = [];
            return;
        end

        % === Score Each Possible Move ===
        % For each neighbor, calculate how good it would be to go there
        scores = zeros(length(neighbors), 1);
        for i = 1:length(neighbors)
            neighbor = neighbors(i);
            
            % Convert node numbers to actual coordinates
            [curr_x, curr_y, curr_layer] = node_to_coords(current, dimX, dimY);
            [neigh_x, neigh_y, neigh_layer] = node_to_coords(neighbor, dimX, dimY);
            [targ_x, targ_y, targ_layer] = node_to_coords(target, dimX, dimY);
            
            % Calculate how far this neighbor is from our target
            dist_to_target = abs(neigh_x - targ_x) + abs(neigh_y - targ_y);
            
            % Count how many times we'll need to switch layers
            layer_switches = 0;
            if curr_layer ~= neigh_layer  % Penalty for switching layers now
                layer_switches = layer_switches + 1;
            end
            if neigh_layer ~= targ_layer  % Penalty if we'll need to switch later
                layer_switches = layer_switches + 1;
            end
            
            % Final score combines:
            % 1. Cost of using this spot (from weights)
            % 2. How far we'll still be from target
            % 3. Penalties for switching layers
            scores(i) = weights(current, neighbor) + dist_to_target + 2*layer_switches;
        end

        % Move to the neighbor with the best (lowest) score
        [~, best_idx] = min(scores);
        current = neighbors(best_idx);
        path = [path current];
        visited(current) = true;
    end

    % If we didn't reach the target, return empty path
    if current ~= target
        path = [];
    end
end

function [x, y, layer] = node_to_coords(node, dimX, dimY)
    % NODE_TO_COORDS - Converts a node number to its actual position in the grid
    %
    % === How It Works ===
    % - First dimX*dimY nodes are on layer 1
    % - Remaining nodes are on layer 2
    % - Within each layer, nodes are numbered left-to-right, then top-to-bottom
    %
    % === Example ===
    % For a 3×2 grid:
    % Layer 1:        Layer 2:
    % 4 5 6          8 10 12
    % 1 2 3          7  9 11
    %
    % === Inputs ===
    % node  : Node number to convert
    % dimX  : Width of each layer
    % dimY  : Height of each layer
    %
    % === Outputs ===
    % x     : Horizontal position (1 to dimX)
    % y     : Vertical position (1 to dimY)
    % layer : Which layer (1 or 2)
    
    grid_size = dimX * dimY;
    
    if node <= grid_size
        % Node is in first layer
        layer = 1;
        [x, y] = ind2sub([dimX, dimY], node);
    else
        % Node is in second layer
        layer = 2;
        node_adj = node - grid_size;
        [x, y] = ind2sub([dimX, dimY], node_adj);
    end
end
