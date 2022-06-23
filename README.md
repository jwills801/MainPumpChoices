This code uses dynamic programming to minimize the size of accumulators on the HHEA. This code takes in the rail flows which resulted from a lagrange multiplier optimization on the pressure rail selection.

main_v1.m is in working condition. This is the code to use. Only use others for constrained optimization. main_v1 does not depend on any files except a drive cycle .mat file.

This code is general. Any number of rails can be used with any net flow sign.

## OLD WAY OF DOING THINGS
This way constrains the options so that each rail recieves uni-directional flow from the pump. The old way is faster but is suboptimal. This way also only works for 2 rails.

The up to data files are M_and_P_DP_five as well as BestFlowRateMP_v3.m

Use BestFlowRateMP_v3.m to perform a grid search and then use M_and_P_DP_five.m to run the simulation with the optimal flow rate - which was obtained through the grid search

There is code at the bottom of BestFlowRateMP_v3.m which plots results for a few different sampling frequencies
