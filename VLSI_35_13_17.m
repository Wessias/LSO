dimX= 35;
dimY = 13;
k=17;
com =[
     3   431;
     5   435;
     6   422;
     7   448;
     8   427;
    10   453;
    11   428;
    14   423;
    17   451;
    21   429;
    23   437;
    24   426;
    25   450;
    27   425;
    30   421;
    33   445;
    35   452];

[h_best, pi_opt, iter, okcom, newnl, best_primal, route] = subgrad_optm_og(dimX, dimY, ...
    k, com, 2000, 2);