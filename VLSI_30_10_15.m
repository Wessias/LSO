dimX=30;
dimY=10;
k=15;
com = [24   271; 
        8   282; 
       21   295; 
       25   290; 
       15   291; 
       27   279; 
       30   293; 
       26   289;  
        2   297;
       29   284; 
       22   276; 
       16   288; 
       19   273; 
       10   275;  
        7   286];

[h_best, pi_opt, iter, okcom, newnl, best_primal, route] = subgrad_optm_og(dimX, dimY, ...
    k, com, 1000, 2);