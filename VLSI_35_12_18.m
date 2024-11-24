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
