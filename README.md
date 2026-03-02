# mp\_modelcnfctls

Replication Files for "Evaluating Monetary Policy Counterfactuals: (When) Do We Need Structural Models?" by Caravello, McKay \& Wolf

Tested in Matlab R2022b on a Dell Inspiron 15.

In order to produce most Figures and Tables, the posterior mode and posterior draws for policy causal effect matrices are required.
Due to its size (~7.5 GB), these are not included in this repository. They must be downloaded from [here](https://www.dropbox.com/scl/fo/zi78j833q1py3w312dm0x/ADee6QczKcb7dqTgQXxRTwY?rlkey=ryb5ilyx0ywjftk2j1lmsz882&e=1&dl=0). Then, the /suff\_stats folder in this repository must be replaced with the one downloaded above.

To ensure that all codes run, the variable "path"---located near the top of the various m-files---needs to be changed to reflect the correct path in the user's machine.

# Contents

[readme.pdf](https://github.com/tcaravello/varplus/blob/main/varplus_readme.pdf): contains detailed instructions for using this repository. The main points are summarized below.

[main\_varplus.m](https://github.com/tcaravello/varplus/blob/main/main_varplus.m): produces all figures and tables in the paper.

[applications](https://github.com/tcaravello/varplus/tree/main/applications): creates Table 5.1 and Figures 6-9, D.1-D.2

* [second\_moments](https://github.com/tcaravello/varplus/tree/main/applications/second_moments): creates Table 5.1 and Figures 6, 7, 8 and D.1-D.2
* [hist\_evol](https://github.com/tcaravello/varplus/tree/main/applications/hist_evol): creates Figure 7
* [hist\_scenario](https://github.com/tcaravello/varplus/tree/main/applications/hist_scenario): creates Figures 8 and 9

[invertibility](https://github.com/tcaravello/varplus/tree/main/invertibility): functions used to create Figure 1 of the paper.

[model\_estim](https://github.com/tcaravello/varplus/tree/main/model_estim): obtains the posterior for the causal effect matrices for all models, as well as the posterior model probabilities.

* [get\_param\_post.m](https://github.com/tcaravello/varplus/tree/main/model_estim/get_param_post.m): obtain descriptive statistics for posterior of strucutral parameters.
* [get\_posterior\_mode.m](https://github.com/tcaravello/varplus/tree/main/model_estim/get_posterior_mode.m): obtains the posterior mode.
* [get\_posterior\_probs.m](https://github.com/tcaravello/varplus/tree/main/model_estim/get_posterior_probs.m): get posterior model probabilities, and create posterior draws across different subsets of models.
* [plot\_model\_irfs.m](https://github.com/tcaravello/varplus/tree/main/model_estim/plot_model_irfs.m): plot fitted impulse responses, as well as impulse responses for different models and horizons.
* [sample\_posterior.m](https://github.com/tcaravello/varplus/tree/main/model_estim/sample_posterior.m): obtain posterior draws for each model.

[var\_inputs](https://github.com/tcaravello/varplus/tree/main/var_inputs): estimate 1) Wold Decomposition, 2) IRF to monetary policy shock using VARs

* [run\_var\_fcst\_evol.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_fcst_evol.m): creates VAR forecasts for each date, used when constructing the historical evolution.
* [run\_var\_fcst\_scenario.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_fcst_scenario.m): creates VAR forecasts at a given point in time, used when constructing the historical scenario.
* [run\_var\_mbc.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_mbc.m): estimates the impulse response to the Main Business Cycle shock by Angeletos, Collard \& Dellas (2020, AER)
* [run\_var\_mp\_adrr.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_mp_adrr.m): estimates the impulse response to monetary policy shocks, using the shock series constructed by Aruoba \& Drechsel (2024) and Romer \& Romer (2004)
* [run\_var\_spf\_fcst\_compare.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_spf_fcst_compare.m): compares the forecasting performance of VARs and the Survey of Profesional Forecasters.
* [run\_var\_swfactors.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_swfactors.m): compares the forecasting performance of VARs with and without the Stock \& Watson (2016) factors.
* [run\_var\_wold.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_wold.m): estimates Wold IRFs for the sample 1960Q1 - 2019Q4
* [run\_var\_wold\_early.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_wold_early.m): estimates Wold IRFs for the sample 1960Q1 - 2007Q1

[\_auxiliary\_functions](https://github.com/tcaravello/varplus/tree/main/_auxiliary_functions): auxiliary functions and routines.

# Acknowledgements

* Our solution of the HANK model uses: first, several files from the replication codes of the article Ahn et al. (2017); second, the CompEcon toolbox of Miranda and Fackler, available here: www4.ncsu.edu/~pfackler/compecon; and third, the ergodicdist.m
  function, written by Marco Maffezzoli. The codes closely build on those used by one of the authors in Wolf (2024).
* The files jbfill.m and winsorize.m are available on Mathworks file exchange.
* The file regcyc.m is taken from the replication codes for the article Hamilton (2018).
