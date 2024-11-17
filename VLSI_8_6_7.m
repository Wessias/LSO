dimX=8;
dimY=6;
k=7;
com = [1 45; 
       2 43; 
       3 44; 
       4 42; 
       5 46; 
       6 47; 
       7 48];


[h_best, pi_opt, iter, okcom, newnl] = subgrad_optm(dimX, dimY, ...
    k, com, 1000, 2);
