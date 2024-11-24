% Example inputs
dimX = 8; dimY = 6; k = 7; delta = 0.9;
x_ergodic = rand(2 * dimX * dimY, 2 * dimX * dimY, k); % Example ergodic solution
com = [1 48; 2 42; 3 43; 4 44; 5 45; 6 46; 7 47]; % Example contact pairs

% Run heuristic
[x_feasible, feas_routes] = primal_feasibility_heuristic_simple(x_ergodic, com, dimX, dimY, delta);

% Display results
disp('Feasible Routes:');
disp(feas_routes);
