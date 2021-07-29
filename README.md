# Lag-augmented Local Projections

Matlab code for inference on impulse responses using local projection or VAR, with or without lag augmentation

**Reference:**
Montiel Olea, José Luis, and Mikkel Plagborg-Møller (2020), "Local Projection Inference is Simpler and More Robust Than You Think", https://scholar.princeton.edu/mikkelpm/publications/lp_inference (paper + supplemental appendix)

Tested in: Matlab R2020a on Windows 10 PC (64-bit) and Matlab R2019a on Macbook Pro 2019

## Contents

**[functions](functions):** Matlab routines
- [ir_estim.m](functions/ir_estim.m): main function for impulse response inference

**[examples](examples):** Empirical examples
- [gk.m](examples/gk.m): example based on Gertler & Karadi (2015)

**[simulations](simulations):** Simulation studies
- [run_sim_ar1.m](simulations/run_sim_ar1.m): AR(1) study
- [run_sim_var.m](simulations/run_sim_var.m): bivariate VAR study
- [run_sim_var_calib.m](simulations/run_sim_var_calib.m): empirically calibrated VAR study
- [create_figs.m](simulations/create_figs.m): create figures of AR(1) and bivariate VAR results
- [create_table.m](simulations/create_table.m): create table of AR(1) and bivariate VAR results
- [create_fig_var_calib.m](simulations/create_fig_var_calib.m): create figure of empirically calibrated VAR results

## Acknowledgements

Plagborg-Møller acknowledges that this material is based upon work supported by the National Science Foundation under Grant #1851665.
