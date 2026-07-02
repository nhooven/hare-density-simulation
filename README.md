## Movement simulation for unmarked density estimation

The aim of this project is to evaluate a continuous-time movement modeling approach to estimating animal travel speed (Noonan et al. 2019) in deriving a key input for the
random encounter model (REM), a popular unmarked density estimator for camera trap detections (Rowcliffe et al. 2008). Our general approach is as follows:

- Fit continuous-time stochastic process (CTSP) movement models to real GPS collar data from snowshoe hares (*Lepus americanus*)
- Extract model parameters
- Initialize simulations with extracted parameters
- Resample simulated trajectories to several fix sampling rates, fix acquisition success rates, and track durations
- Fit CTSP models to resampled data
- Extract speed estimates
- Use speeds as inputs to REMs, and evaluate performance compared to "true" animal density

This workflow heavily depends on the "ctmm" package in R (Calabrese et al. 2016).

A preprint describing the results of this work is available at: <https://doi.org/10.2139/ssrn.6293639>, and you can find a complementary repository including code
to conduct this approach on real camera trap data at <https://github.com/nhooven/hare-density-cameras>.


## References
Calabrese, J. M., C. H. Fleming, and E. Gurarie. 2016. ctmm: an R package for analyzing animal relocation data as a continuous‐time stochastic process. 
Methods in Ecology and Evolution 7:1124–1132.

Noonan, M. J., C. H. Fleming, T. S. Akre, J. Drescher-Lehman, E. Gurarie, A.-L. Harrison, R. Kays, and J. M. Calabrese. 2019. 
Scale-insensitive estimation of speed and distance traveled from animal tracking data. Movement Ecology 7:35.

Rowcliffe, J. M., J. Field, S. T. Turvey, and C. Carbone. 2008. Estimating animal density using camera traps without the need for individual recognition. 
Journal of Applied Ecology 45:1228–1236.



