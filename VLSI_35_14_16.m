dimX= 35;
dimY = 14;
k=16;

com = [ 3   470;
     5   475;
     7   462;
     8   487;
     9   464;
    12   485;
    13   465;
    16   467;
    19   459;
    23   490;
    25   468;
    26   477;
    27   466;
    29   488;
    32   461;
    35   460];


[h_best, pi_opt, iter, okcom, newnl, best_primal, route] = subgrad_optm_og(dimX, dimY, ...
    k, com, 5000, 2);
