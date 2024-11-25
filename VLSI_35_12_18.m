dimX= 35;
dimY = 12;
k=18;
com =[
     2   394;
     4   398;
     5   389;
     6   410;
     7   414;
     9   413;
    10   387;
    13   416;
    16   412;
    20   399;
    22   392;
    23   411;
    24   391;
    26   415;
    28   407;
    31   390;
    34   393;
    35   386];


[h_best, pi_opt, iter, okcom, newnl, best_primal, route] = subgrad_optm_og(dimX, dimY, ...
    k, com, 1000, 2);


% Initialize an empty list to store all nodes
all_nodes_used = [];
empty_routes = []; % To store indices of empty or invalid routes

% Iterate through feasible_routes
for i = 1:length(route)
    curRoute = route{i}; % Get the current route
    if isempty(curRoute) || isequal(curRoute, 0)
        % Record the index of the empty or invalid route
        empty_routes = [empty_routes; i];
    else
        % Concatenate nodes from the current route as a column vector
        all_nodes_used = [all_nodes_used; curRoute(:)]; % Ensure curRoute is a column
    end
end



% Remove elements of com corresponding to empty_routes
com_updated = com; % Create a copy of com
com_updated(empty_routes, :) = []; % Delete rows at indices in empty_routes



%visagrid(dimX, dimY, all_nodes_used, com_updated, pi_opt, 25);
