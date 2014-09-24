function paper_add_cse_perf()

FAST423_20_WO_NO_CSE_1 = [2000 2000 2000 1 678.37 ;  2500 2500 2500 1 1382.23 ;  3000 3000 3000 1 2224.71 ;  3500 3500 3500 1 3418.86 ;  4000 4000 4000 1 4991.19 ;  4500 4500 4500 1 7038.42 ;  5000 5000 5000 1 9603.15 ;  5500 5500 5500 1 13145 ;  6000 6000 6000 1 16617.2 ;  6500 6500 6500 1 21235.8 ; ];
FAST423_20_WO_YES_CSE_1 = [2000 2000 2000 1 677.131 ;  2500 2500 2500 1 1298.87 ;  3000 3000 3000 1 2211.32 ;  3500 3500 3500 1 3433.93 ;  4000 4000 4000 1 5005.58 ;  4500 4500 4500 1 7206.57 ;  5000 5000 5000 1 9656.99 ;  5500 5500 5500 1 12945.4 ;  6000 6000 6000 1 16708.8 ;  6500 6500 6500 1 21193 ; ];
FAST423_20_ST_NO_CSE_1 = [2000 2000 2000 1 698.274 ;  2500 2500 2500 1 1361.02 ;  3000 3000 3000 1 2252.99 ;  3500 3500 3500 1 3543.69 ;  4000 4000 4000 1 5033.19 ;  4500 4500 4500 1 7260.45 ;  5000 5000 5000 1 9897.89 ;  5500 5500 5500 1 13246.6 ;  6000 6000 6000 1 16869.1 ;  6500 6500 6500 1 21760.7 ; ];
FAST423_20_ST_YES_CSE_1 = [2000 2000 2000 1 681.314 ;  2500 2500 2500 1 1336.34 ;  3000 3000 3000 1 2318.26 ;  3500 3500 3500 1 3491.52 ;  4000 4000 4000 1 4946.46 ;  4500 4500 4500 1 7258.42 ;  5000 5000 5000 1 9635.18 ;  5500 5500 5500 1 13083.1 ;  6000 6000 6000 1 16682.8 ;  6500 6500 6500 1 21552.7 ; ];
FAST423_20_PW_NO_CSE_1 = [2000 2000 2000 1 702.407 ;  2500 2500 2500 1 1325.51 ;  3000 3000 3000 1 2247.61 ;  3500 3500 3500 1 3546.19 ;  4000 4000 4000 1 4985.58 ;  4500 4500 4500 1 7285.62 ;  5000 5000 5000 1 9798.67 ;  5500 5500 5500 1 13409.3 ;  6000 6000 6000 1 16986.7 ;  6500 6500 6500 1 21774.5 ; ];
FAST423_20_PW_YES_CSE_1 = [2000 2000 2000 1 687.902 ;  2500 2500 2500 1 1326.59 ;  3000 3000 3000 1 2224.28 ;  3500 3500 3500 1 3532.14 ;  4000 4000 4000 1 4977.95 ;  4500 4500 4500 1 7364.75 ;  5000 5000 5000 1 9836.04 ;  5500 5500 5500 1 13239.7 ;  6000 6000 6000 1 16925.7 ;  6500 6500 6500 1 21677.9 ; ];

FAST423_20_WO_NO_CSE_2 = [2000 2000 2000 2 868.512 ;  2500 2500 2500 2 1399.36 ;  3000 3000 3000 2 2379.56 ;  3500 3500 3500 2 3651.47 ;  4000 4000 4000 2 4811.77 ;  4500 4500 4500 2 7164.1 ;  5000 5000 5000 2 9507.45 ;  5500 5500 5500 2 13107.4 ;  6000 6000 6000 2 15941.2 ;  6500 6500 6500 2 20589.2 ; ];
FAST423_20_WO_YES_CSE_2 = [2000 2000 2000 2 700.978 ;  2500 2500 2500 2 1364.22 ;  3000 3000 3000 2 2359.28 ;  3500 3500 3500 2 3662.55 ;  4000 4000 4000 2 5115.62 ;  4500 4500 4500 2 7389.86 ;  5000 5000 5000 2 9762.05 ;  5500 5500 5500 2 13515.6 ;  6000 6000 6000 2 16215.3 ;  6500 6500 6500 2 20804.1 ; ];
FAST423_20_ST_NO_CSE_2 = [2000 2000 2000 2 686.128 ;  2500 2500 2500 2 1446.66 ;  3000 3000 3000 2 2640.29 ;  3500 3500 3500 2 3798.94 ;  4000 4000 4000 2 5420.79 ;  4500 4500 4500 2 7698.72 ;  5000 5000 5000 2 10497.5 ;  5500 5500 5500 2 14054 ;  6000 6000 6000 2 16834.3 ;  6500 6500 6500 2 21892.2 ; ];
FAST423_20_ST_YES_CSE_2 = [2000 2000 2000 2 672.876 ;  2500 2500 2500 2 1370.41 ;  3000 3000 3000 2 2389.39 ;  3500 3500 3500 2 3679.21 ;  4000 4000 4000 2 5317.85 ;  4500 4500 4500 2 7635.78 ;  5000 5000 5000 2 9906.77 ;  5500 5500 5500 2 13938.1 ;  6000 6000 6000 2 16662.9 ;  6500 6500 6500 2 21723.6 ; ];
FAST423_20_PW_NO_CSE_2 = [2000 2000 2000 2 779.054 ;  2500 2500 2500 2 1489.09 ;  3000 3000 3000 2 2511.32 ;  3500 3500 3500 2 3894.35 ;  4000 4000 4000 2 5210.9 ;  4500 4500 4500 2 7613.58 ;  5000 5000 5000 2 10164.4 ;  5500 5500 5500 2 13848.5 ;  6000 6000 6000 2 16913.7 ;  6500 6500 6500 2 21639.2 ; ];
FAST423_20_PW_YES_CSE_2 = [ 2000 2000 2000 2 746.895 ;  2500 2500 2500 2 1432.23 ;  3000 3000 3000 2 2457.66 ;  3500 3500 3500 2 3861.9 ;  4000 4000 4000 2 5259.71 ;  4500 4500 4500 2 7682.77 ;  5000 5000 5000 2 10142.4 ;  5500 5500 5500 2 13991.6 ;  6000 6000 6000 2 16677.5 ;  6500 6500 6500 2 21680.2 ; ];

FAST424_26_WO_NO_CSE_1 = [2000 1600 2000 1 604.009 ;  2500 1600 2500 1 865.436 ;  3000 1600 3000 1 1199.57 ;  3500 1600 3500 1 1629.08 ;  4000 1600 4000 1 2047.28 ;  4500 1600 4500 1 2641.32 ;  5000 1600 5000 1 3170.31 ;  5500 1600 5500 1 3978.41 ;  6000 1600 6000 1 4578.85 ;  6500 1600 6500 1 5448.43 ;  7000 1600 7000 1 6276.99 ;  7500 1600 7500 1 7336.91 ;  8000 1600 8000 1 7954.83 ;  8500 1600 8500 1 9347.94 ;  9000 1600 9000 1 10562.8 ; ];
FAST424_26_WO_YES_CSE_1 = [2000 1600 2000 1 570.221 ;  2500 1600 2500 1 877.696 ;  3000 1600 3000 1 1248.12 ;  3500 1600 3500 1 1727.09 ;  4000 1600 4000 1 2112.59 ;  4500 1600 4500 1 2773.55 ;  5000 1600 5000 1 3371.21 ;  5500 1600 5500 1 4080.9 ;  6000 1600 6000 1 4865.19 ;  6500 1600 6500 1 5720.37 ;  7000 1600 7000 1 6489.87 ;  7500 1600 7500 1 7565.03 ;  8000 1600 8000 1 8290.21 ;  8500 1600 8500 1 9586.17 ;  9000 1600 9000 1 10620.2 ; ];
FAST424_26_ST_NO_CSE_1 = [2000 1600 2000 1 589.652 ;  2500 1600 2500 1 985.629 ;  3000 1600 3000 1 1309.07 ;  3500 1600 3500 1 1928.78 ;  4000 1600 4000 1 2177.2 ;  4500 1600 4500 1 2876.49 ;  5000 1600 5000 1 3343.82 ;  5500 1600 5500 1 4230.26 ;  6000 1600 6000 1 4740.78 ;  6500 1600 6500 1 5850.58 ;  7000 1600 7000 1 6408.17 ;  7500 1600 7500 1 7796.94 ;  8000 1600 8000 1 8147.35 ;  8500 1600 8500 1 9520.9 ;  9000 1600 9000 1 10969.9 ; ];
FAST424_26_ST_YES_CSE_1 = [2000 1600 2000 1 582.397 ;  2500 1600 2500 1 969.107 ;  3000 1600 3000 1 1309.56 ;  3500 1600 3500 1 1923.24 ;  4000 1600 4000 1 2186.45 ;  4500 1600 4500 1 2862.65 ;  5000 1600 5000 1 3320.99 ;  5500 1600 5500 1 4210.78 ;  6000 1600 6000 1 4714.97 ;  6500 1600 6500 1 5798.03 ;  7000 1600 7000 1 6386.26 ;  7500 1600 7500 1 7780.19 ;  8000 1600 8000 1 8094.15 ;  8500 1600 8500 1 9383.32 ;  9000 1600 9000 1 10789.6 ; ];
FAST424_26_PW_NO_CSE_1 = [2000 1600 2000 1 613.723 ;  2500 1600 2500 1 939.396 ;  3000 1600 3000 1 1311.66 ;  3500 1600 3500 1 1754.32 ;  4000 1600 4000 1 2204.29 ;  4500 1600 4500 1 2931.82 ;  5000 1600 5000 1 3538.5 ;  5500 1600 5500 1 4232.88 ;  6000 1600 6000 1 5011.4 ;  6500 1600 6500 1 5940.95 ;  7000 1600 7000 1 6799 ;  7500 1600 7500 1 7974.47 ;  8000 1600 8000 1 8599.29 ;  8500 1600 8500 1 9984.27 ;  9000 1600 9000 1 11055 ; ];
FAST424_26_PW_YES_CSE_1 = [2000 1600 2000 1 598.03 ;  2500 1600 2500 1 914.23 ;  3000 1600 3000 1 1288.37 ;  3500 1600 3500 1 1761.32 ;  4000 1600 4000 1 2185.33 ;  4500 1600 4500 1 2875.43 ;  5000 1600 5000 1 3496.84 ;  5500 1600 5500 1 4233.65 ;  6000 1600 6000 1 5081.56 ;  6500 1600 6500 1 5914.07 ;  7000 1600 7000 1 6753.16 ;  7500 1600 7500 1 7916.57 ;  8000 1600 8000 1 8570.9 ;  8500 1600 8500 1 10015.5 ;  9000 1600 9000 1 11057.2 ; ];

FAST424_26_WO_NO_CSE_2 = [2000 1600 2000 2 721.981 ;  2500 1600 2500 2 1005.3 ;  3000 1600 3000 2 1407.9 ;  3500 1600 3500 2 1809.28 ;  4000 1600 4000 2 2162.66 ;  4500 1600 4500 2 2764.06 ;  5000 1600 5000 2 3238.16 ;  5500 1600 5500 2 4072.11 ;  6000 1600 6000 2 4573.51 ;  6500 1600 6500 2 5446.92 ;  7000 1600 7000 2 6362.4 ;  7500 1600 7500 2 7165.41 ;  8000 1600 8000 2 7852.82 ;  8500 1600 8500 2 9269.53 ;  9000 1600 9000 2 10678.4 ; ];
FAST424_26_WO_YES_CSE_2 = [2000 1600 2000 2 615.982 ;  2500 1600 2500 2 1051.31 ;  3000 1600 3000 2 1555.03 ;  3500 1600 3500 2 2029.64 ;  4000 1600 4000 2 2392.27 ;  4500 1600 4500 2 3019.78 ;  5000 1600 5000 2 3654.26 ;  5500 1600 5500 2 4425.31 ;  6000 1600 6000 2 5013.05 ;  6500 1600 6500 2 6251.37 ;  7000 1600 7000 2 7255.92 ;  7500 1600 7500 2 8252.89 ;  8000 1600 8000 2 9043.68 ;  8500 1600 8500 2 10486.8 ;  9000 1600 9000 2 11750.1 ; ];
FAST424_26_ST_NO_CSE_2 = [2000 1600 2000 2 680.275 ;  2500 1600 2500 2 1308.72 ;  3000 1600 3000 2 1618.97 ;  3500 1600 3500 2 2293.46 ;  4000 1600 4000 2 2444.12 ;  4500 1600 4500 2 3301.08 ;  5000 1600 5000 2 3609.57 ;  5500 1600 5500 2 4660.9 ;  6000 1600 6000 2 5014.44 ;  6500 1600 6500 2 6680.63 ;  7000 1600 7000 2 7430.83 ;  7500 1600 7500 2 8535.65 ;  8000 1600 8000 2 8902.19 ;  8500 1600 8500 2 10294.6 ;  9000 1600 9000 2 11919.8 ; ];
FAST424_26_ST_YES_CSE_2 = [2000 1600 2000 2 667.81 ;  2500 1600 2500 2 1224.76 ;  3000 1600 3000 2 1593.06 ;  3500 1600 3500 2 2261.03 ;  4000 1600 4000 2 2397.11 ;  4500 1600 4500 2 3271.36 ;  5000 1600 5000 2 3553.64 ;  5500 1600 5500 2 4621.84 ;  6000 1600 6000 2 4996.93 ;  6500 1600 6500 2 6573.5 ;  7000 1600 7000 2 7335.67 ;  7500 1600 7500 2 8445.64 ;  8000 1600 8000 2 8779.12 ;  8500 1600 8500 2 10118.1 ;  9000 1600 9000 2 11575.9 ; ];
FAST424_26_PW_NO_CSE_2 = [2000 1600 2000 2 810.816 ;  2500 1600 2500 2 1226.92 ;  3000 1600 3000 2 1701.55 ;  3500 1600 3500 2 2235.59 ;  4000 1600 4000 2 2668.2 ;  4500 1600 4500 2 3471.58 ;  5000 1600 5000 2 4065.79 ;  5500 1600 5500 2 4883.73 ;  6000 1600 6000 2 5607.35 ;  6500 1600 6500 2 6645.19 ;  7000 1600 7000 2 7636.22 ;  7500 1600 7500 2 8327.57 ;  8000 1600 8000 2 9639.73 ;  8500 1600 8500 2 11185.9 ;  9000 1600 9000 2 12457.9 ; ];
FAST424_26_PW_YES_CSE_2 = [2000 1600 2000 2 780.048 ;  2500 1600 2500 2 1197.83 ;  3000 1600 3000 2 1698.77 ;  3500 1600 3500 2 2232.17 ;  4000 1600 4000 2 2645.78 ;  4500 1600 4500 2 3309.35 ;  5000 1600 5000 2 4025.99 ;  5500 1600 5500 2 4819.72 ;  6000 1600 6000 2 5490.71 ;  6500 1600 6500 2 6808.92 ;  7000 1600 7000 2 7906.84 ;  7500 1600 7500 2 8984.9 ;  8000 1600 8000 2 9925.24 ;  8500 1600 8500 2 11411.4 ;  9000 1600 9000 2 12735.5 ; ];

close all;
add_cse_perf_plot(FAST424_26_WO_NO_CSE_1, FAST424_26_WO_YES_CSE_1, ...
                  FAST424_26_ST_NO_CSE_1, FAST424_26_ST_YES_CSE_1, ...
                  FAST424_26_PW_NO_CSE_1, FAST424_26_PW_YES_CSE_1, ...
                  '<4,2,4>, one recursive step', 'adds_cse_424_1', 0);

add_cse_perf_plot(FAST424_26_WO_NO_CSE_2, FAST424_26_WO_YES_CSE_2, ...
                  FAST424_26_ST_NO_CSE_2, FAST424_26_ST_YES_CSE_2, ...
                  FAST424_26_PW_NO_CSE_2, FAST424_26_PW_YES_CSE_2, ...
                  '<4,2,4>, two recursive steps', 'adds_cse_424_2', 0);
              
add_cse_perf_plot(FAST423_20_WO_NO_CSE_1, FAST423_20_WO_YES_CSE_1, ...
                  FAST423_20_ST_NO_CSE_1, FAST423_20_ST_YES_CSE_1, ...
                  FAST423_20_PW_NO_CSE_1, FAST423_20_PW_YES_CSE_1, ...
                  '<4,2,3>, one recursive step', 'adds_cse_423_1', 0);

add_cse_perf_plot(FAST423_20_WO_NO_CSE_2, FAST423_20_WO_YES_CSE_2, ...
                  FAST423_20_ST_NO_CSE_2, FAST423_20_ST_YES_CSE_2, ...
                  FAST423_20_PW_NO_CSE_2, FAST423_20_PW_YES_CSE_2, ...
                  '<4,2,3>, two recursive steps', 'adds_cse_423_2', 1);

              