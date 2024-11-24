dimX=30;
dimY=11;
k=15;

com = [24   301; 
        8   312; 
       21   325; 
       25   320; 
       15   321; 
       27   309; 
       30   323; 
       26   319;  
        2   327;
       29   314; 
       22   306; 
       16   318; 
       19   303; 
       10   305;  
        7   316];

[h_best, pi_opt, iter, okcom, newnl, best_primal, route] = subgrad_optm_og(dimX, dimY, ...
    k, com, 2000, 2);
