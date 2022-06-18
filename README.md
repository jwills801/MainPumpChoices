This code uses dynamic programming to minimize the size of accumulators on the HHEA. This code takes in the rail flows which resulted from a lagrange multiplier optimization on the pressure rail selection.


The up to data files are M_and_P_DP_five as well as BestFlowRateMP_v3.m

Use BestFlowRateMP_v3.m to perform a grid search and then use M_and_P_DP_five.m to run the simulation with the optimal flow rate - which was obtained through the grid search

There is code at the bottom of BestFlowRateMP_v3.m which plots results for a few different sampling frequencies


## NOTE!
Now I am tying to generalize the algorithm. The code to do this is in main_v1.m

This code is still a work in progress
