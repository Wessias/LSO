function compare_performance()
    % COMPARE_PERFORMANCE - Compares the performance of subgrad_optm and subgrad_optm2
    % Tests both implementations on different problem sizes and measures:
    % 1. Execution time
    % 2. Solution quality
    % 3. Convergence speed
    
    % Test cases
    test_cases = {
        'VLSI_8_6_7.m',      % Small problem
        'VLSI_30_10_15.m',   % Medium problem
        'VLSI_30_11_15.m',   % Medium problem
        'VLSI_35_12_18.m',   % Large problem
        'VLSI_35_13_17.m',   % Large problem
        'VLSI_35_14_16.m'    % Large problem
    };
    
    % Parameters for optimization
    max_iter = 1000;
    theta_init = 2;
    
    % Initialize results storage
    n_cases = length(test_cases);
    results = struct();
    
    % Create figure for plotting
    figure('Position', [100, 100, 1200, 800]);
    
    % Process each test case
    for i = 1:n_cases
        test_case = test_cases{i};
        fprintf('\nProcessing test case: %s\n', test_case);
        
        % Load problem data
        run(test_case);
        fprintf('Problem dimensions: %d x %d, %d pairs\n', dimX, dimY, k);
        
        % Test original implementation
        tic;
        [h_best1, pi_opt1, iter1, okcom1, newnl1] = subgrad_optm(dimX, dimY, k, com, max_iter, theta_init);
        time1 = toc;
        
        % Test new implementation
        tic;
        [h_best2, pi_opt2, iter2, okcom2, newnl2, best_primal] = subgrad_optm2(dimX, dimY, k, com, max_iter, theta_init);
        time2 = toc;
        
        % Store results
        results(i).test_case = test_case;
        results(i).dimensions = sprintf('%dx%d', dimX, dimY);
        results(i).pairs = k;
        results(i).original_time = time1;
        results(i).new_time = time2;
        results(i).original_obj = h_best1;
        results(i).new_obj = h_best2;
        results(i).original_pairs = length(okcom1);
        results(i).new_pairs = best_primal;
        results(i).speedup = time1/time2;
        
        % Print results for this case
        fprintf('\nResults for %s:\n', test_case);
        fprintf('Original Implementation:\n');
        fprintf('  Time: %.2f seconds\n', time1);
        fprintf('  Objective: %.4f\n', h_best1);
        fprintf('  Feasible pairs: %d\n', length(okcom1));
        
        fprintf('\nNew Implementation:\n');
        fprintf('  Time: %.2f seconds\n', time2);
        fprintf('  Objective: %.4f\n', h_best2);
        fprintf('  Feasible pairs: %d\n', best_primal);
        fprintf('  Speedup: %.2fx\n', time1/time2);
    end
    
    % Create summary plots
    subplot(2,2,1);
    bar([results.original_time; results.new_time]');
    title('Execution Time Comparison');
    xlabel('Test Case');
    ylabel('Time (seconds)');
    legend('Original', 'New');
    set(gca, 'XTickLabel', {results.dimensions});
    grid on;
    
    subplot(2,2,2);
    bar([results.speedup]);
    title('Speedup Factor (Original/New)');
    xlabel('Test Case');
    ylabel('Speedup');
    set(gca, 'XTickLabel', {results.dimensions});
    grid on;
    
    subplot(2,2,3);
    bar([results.original_obj; results.new_obj]');
    title('Objective Value Comparison');
    xlabel('Test Case');
    ylabel('Objective Value');
    legend('Original', 'New');
    set(gca, 'XTickLabel', {results.dimensions});
    grid on;
    
    subplot(2,2,4);
    bar([results.original_pairs; results.new_pairs]');
    title('Feasible Pairs Comparison');
    xlabel('Test Case');
    ylabel('Number of Feasible Pairs');
    legend('Original', 'New');
    set(gca, 'XTickLabel', {results.dimensions});
    grid on;
    
    % Save results to file
    fid = fopen('performance_comparison.txt', 'w');
    fprintf(fid, 'Performance Comparison Results\n');
    fprintf(fid, '============================\n\n');
    
    for i = 1:n_cases
        fprintf(fid, 'Test Case: %s\n', results(i).test_case);
        fprintf(fid, 'Dimensions: %s\n', results(i).dimensions);
        fprintf(fid, 'Number of pairs: %d\n', results(i).pairs);
        fprintf(fid, '\nOriginal Implementation:\n');
        fprintf(fid, '  Time: %.2f seconds\n', results(i).original_time);
        fprintf(fid, '  Objective: %.4f\n', results(i).original_obj);
        fprintf(fid, '  Feasible pairs: %d\n', results(i).original_pairs);
        fprintf(fid, '\nNew Implementation:\n');
        fprintf(fid, '  Time: %.2f seconds\n', results(i).new_time);
        fprintf(fid, '  Objective: %.4f\n', results(i).new_obj);
        fprintf(fid, '  Feasible pairs: %d\n', results(i).new_pairs);
        fprintf(fid, '  Speedup: %.2fx\n', results(i).speedup);
        fprintf(fid, '\n----------------------------\n\n');
    end
    
    % Calculate and write summary statistics
    avg_speedup = mean([results.speedup]);
    avg_obj_improvement = mean(([results.new_obj] - [results.original_obj]) ./ [results.original_obj] * 100);
    avg_pairs_improvement = mean(([results.new_pairs] - [results.original_pairs]) ./ [results.original_pairs] * 100);
    
    fprintf(fid, '\nSummary Statistics:\n');
    fprintf(fid, '==================\n');
    fprintf(fid, 'Average speedup: %.2fx\n', avg_speedup);
    fprintf(fid, 'Average objective improvement: %.2f%%\n', avg_obj_improvement);
    fprintf(fid, 'Average feasible pairs improvement: %.2f%%\n', avg_pairs_improvement);
    
    fclose(fid);
    
    % Display summary in console
    fprintf('\nPerformance comparison completed. Results saved to performance_comparison.txt\n');
    fprintf('Average speedup across all test cases: %.2fx\n', avg_speedup);
    fprintf('Average objective improvement: %.2f%%\n', avg_obj_improvement);
    fprintf('Average feasible pairs improvement: %.2f%%\n', avg_pairs_improvement);
end
