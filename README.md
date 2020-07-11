# LALP

Matlab code for inference on impulse responses using local projection or VAR, with or without lag augmentation

**Reference:**
Montiel Olea, José Luis, and Mikkel Plagborg-Møller (2020), "Local Projection Inference is Simpler and More Robust Than You Think", https://scholar.princeton.edu/mikkelpm/publications/lp_inference (paper + supplemental appendix)

Tested in Matlab R2020a on Windows 10 PC, 64-bit

## Contents

**[functions](functions):** Matlab routines
- [functions/ir_estim.m](ir_estim.m): main function for impulse response inference

**[examples](examples):** Empirical examples
- [examples/gk.m](gk.m): example based on Gertler & Karadi (2015)

**[simulations](simulations):** Simulation studies
- [simulations/run_sim_ar1.m](run_sim_ar1.m): AR(1) study
- [simulations/run_sim_var.m](run_sim_var.m): VAR(p) study
- [simulations/create_figs.m](create_figs.m): create figures of results
