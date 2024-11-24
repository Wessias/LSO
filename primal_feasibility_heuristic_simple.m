function [x_feasible, feasible_routes] = primal_feasibility_heuristic_simple(x_ergodic, com, dimX, dimY, delta)
    % PRIMAL_FEASIBILITY_HEURISTIC
    % Constructs a feasible primal solution based on the ergodic solution.
    % Logical paths are excluded from routing.
    %
    % INPUT:
    %   x_ergodic      - Ergodic solution (dimX * dimY * k matrix)
    %   com            - Contact pairs [k x 2 matrix]
    %   dimX, dimY     - Graph dimensions
    %   delta          - Threshold for binary conversion
    %
    % OUTPUT:
    %   x_feasible     - Feasible primal solution (binary matrix)
    %   feasible_routes - Feasible routes for each contact pair

    % Initialize outputs
    k = size(com, 1);  % Number of contact pairs
    x_feasible = zeros(size(x_ergodic));  % Feasible primal solution
    feasible_routes = cell(k, 1);        % Store routes for each pair

    % Step 1: Threshold x_ergodic to create binary solution
    x_binary = x_ergodic > delta;

    % Debugging: Display thresholded matrix
    disp('Thresholded x_binary for each contact pair:');
    for l = 1:k
        fprintf('Contact pair %d:\n', l);
        disp(x_binary(:, :, l));
    end

    % Step 2: Validate Paths for Each Contact Pair
    for l = 1:k
        start_node = com(l, 1);
        end_node = com(l, 2);

        % Debugging: Print nodes
        fprintf('Processing contact pair %d: start_node=%d, end_node=%d\n', l, start_node, end_node);

        % Extract path from binary solution
        route = extract_path(x_binary(:, :, l), start_node, end_node);

        % Debugging: Print extracted route
        fprintf('Extracted route for pair %d: ', l);
        disp(route);

        % Validate route
        if is_valid_route(route, start_node, end_node)
            % Debugging: Valid route
            fprintf('Valid route found for pair %d: ', l);
            disp(route);

            % Valid route; add to feasible solution
            feasible_routes{l} = route;
            for idx = 1:(length(route)-1)
                i = route(idx);
                j = route(idx+1);
                x_feasible(i, j, l) = 1;
            end
        else
            % Debugging: Invalid route
            fprintf('No valid route for pair %d\n', l);
        end
    end

    % Step 3: Prune Paths to Ensure Feasibility
    x_feasible = prune_paths(x_feasible, feasible_routes, dimX, dimY);

    % Debugging: Final feasible routes
    fprintf('Final feasible routes:\n');
    disp(feasible_routes);

    % Nested helper function to extract a path
    function route = extract_path(x_binary_l, start_node, end_node)
        % EXTRACT_PATH
        % Extracts a route from binary x_binary for a single contact pair.
        %
        % INPUT:
        %   x_binary_l  - Binary solution for a single contact pair (matrix)
        %   start_node  - Start node of the contact pair
        %   end_node    - End node of the contact pair
        % OUTPUT:
        %   route       - Extracted route as a list of nodes (physical only)

        route = []; % Initialize route
        current_node = start_node;

        while current_node ~= end_node
            % Find outgoing arcs (excluding logical arcs)
            next_nodes = find(x_binary_l(current_node, :) > 0);

            % Debugging: Print potential next nodes
            fprintf('Current node: %d, Next nodes: ', current_node);
            disp(next_nodes);

            if isempty(next_nodes)
                route = []; % No valid path
                return;
            end

            % Assume first valid node
            route = [route, current_node];
            current_node = next_nodes(1);
        end

        route = [route, end_node]; % Complete the route
    end

    % Nested helper function to validate a route
    function is_valid = is_valid_route(route, start_node, end_node)
        % IS_VALID_ROUTE
        % Checks if the route is valid and uses only physical arcs.
        %
        % INPUT:
        %   route       - List of nodes representing the route
        %   start_node  - Start node of the contact pair
        %   end_node    - End node of the contact pair
        % OUTPUT:
        %   is_valid    - Boolean indicating if the route is valid

        is_valid = ~isempty(route) && route(1) == start_node && route(end) == end_node;
    end

    % Nested helper function to prune conflicting paths
    function x_pruned = prune_paths(x_feasible, feasible_routes, dimX, dimY)
        % PRUNE_PATHS
        % Ensures no node is used by more than one path.
        %
        % INPUT:
        %   x_feasible      - Current feasible primal solution
        %   feasible_routes - Routes for each contact pair
        %   dimX, dimY      - Graph dimensions
        % OUTPUT:
        %   x_pruned        - Pruned feasible primal solution

        x_pruned = x_feasible;
        used_nodes = zeros(dimX * dimY * 2, 1); % Track node usage

        for l = 1:length(feasible_routes)
            if isempty(feasible_routes{l})
                continue; % Skip if no route
            end

            route = feasible_routes{l};
            for idx = 1:length(route)
                node = route(idx);
                if used_nodes(node) > 0
                    % Conflict detected; prune route
                    x_pruned(:, :, l) = 0;
                    feasible_routes{l} = []; % Mark as invalid
                    break;
                end
                used_nodes(node) = used_nodes(node) + 1;
            end
        end
    end
end
