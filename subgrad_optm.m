function [h_best, pi_opt, iter, okcom, newnl] = subgrad_optm(dimX, dimY, k, com, max_iter, theta_init)
    % SUBGRADIENT_OPTIMIZATION
    % Solves the Lagrangean dual for VLSI problem using subgradient optimization.
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
    pi = zeros(dimX * dimY * 2, 1);  % Initialize dual variables
    h_best = inf;  % Best dual value found
    h_star = 0;     % Lower bound of the dual problem (set to 0)
    theta = theta_init;  % Initialize step size multiplier
    

    % Iterate for subgradient optimization
    for iter = 1:max_iter
        % Step 1: Solve Lagrangean subproblem using gsp
        nl = gsp(dimX, dimY, pi, k, com);
        

        
        % Calculate cost per route; remove route with
        % cost > 1 (required routes stored in nl and pairs in com)
        h_pi = 0; % Init 
        okcom = [];
        newnl = [];

        last = 0;
        for i = 1 : k;
            first = last+1;
            slask = find(nl(last+1:length(nl)) == com(i,1));
            last = slask(1)+first-1;
            cost = sum(pi(nl(first:last)));
            if (cost < 1)
                okcom = [okcom i]; newnl = [newnl; nl(first:last)]; h_pi = h_pi + 1 - cost; %Wrong?
            end
        end

        %STILL HAVE TO FINISH CALC OF h_pi
        h_pi = h_pi + sum(pi);
        
        
        % Update best dual value
        h_best = min(h_best, h_pi);
        
        % Step 2: Compute subgradients
        gamma = zeros(size(pi));  % Initialize subgradient
        for i = 1:length(pi)
            gamma(i) = 1 - sum(newnl == i);  % Count routes passing through node i
            %CHECK IF THIS MAKES SENSE (think it does)
        end
        
        % Step 3: Calculate step length
        step_length = max(1e-4, theta * (h_pi - h_best) / norm(gamma)^2);

        fprintf("Step length %d \n", step_length);
        
        % Step 4: Update dual variables
        pi = max(0, pi + step_length * gamma);  % Projection ensures pi >= 0
        

        % Step 5: Update theta periodically
        if mod(iter, 10) == 0
            theta = theta * 0.95;
        end
        
        % Print progress
        fprintf('Iteration %d: h(pi) = %.4f\n', iter, h_pi);
        
        
    end
    
    % Output optimized dual variables
    pi_opt = pi;
    fprintf('Optimization completed in %d iterations.\n', iter);
    fprintf('Best dual objective value: %.4f\n', h_best);
end
