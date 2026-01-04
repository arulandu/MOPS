# MOPS
This repository implements https://arxiv.org/pdf/math-ph/0409066 in Julia (see `./julia`). This port was unit tested against outputs of the original Maple implementation, in `./maple`. 

Note that `GBC_cont` is the ported original implementation of the contiguous generalized binomial coefficient calculation. This has been shown to be incorrect. Thus, the port uses `GBC_cont_explicit` instead which is faithful to the paper's algorithm. 

Despite this, there are errors in the actual algorithm of the paper itself that prevent this implementation from generating valid Jacobi polynomials. 
