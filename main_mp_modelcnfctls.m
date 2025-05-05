%% REPLICATION FILES: Evaluating Monetary Policy Counterfactuals: (When) Do We Need Structural Models?
% Tomas Caravello, Alisdair McKay, and Christian Wolf
% this version: 05/05/2025

%% HOUSEKEEPING
 
clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')
warning('off','MATLAB:declareGlobalBeforeUse')

path = '/Users/tomyc/Dropbox (MIT)/mp_modelcnfctls/code/github_public/mp_modelcnfctls';
vintage = '';

save_fig = 0; 

addpath([path vintage '/_auxiliary_functions'])
addpath([path vintage '/var_inputs'])
addpath([path vintage '/var_inputs/_results'])
addpath([path vintage '/model_estim'])
addpath([path vintage '/applications/second_moments'])
addpath([path vintage '/applications/hist_scenario'])
addpath([path vintage '/applications/hist_evol'])
addpath([path vintage '/invertibility/sw_sequence'])
addpath([path vintage '/suff_stats/ratex']);
addpath([path vintage '/suff_stats/behavioral']);
addpath([path vintage '/suff_stats/mix']);

%% THEORETICAL INVERTIBILITY ANALYSIS (FIGURE 1)

% counterfactuals (true and approximate)

get_2mom_ninvwold_all;

% forecasting

get_2mom_fcstvars;

%% MONETARY POLICY IRFs

%----------------------------------------------------------------
% Figure 2
%----------------------------------------------------------------

run_var_mp_adrr;
%%
%----------------------------------------------------------------
% Figures 3-5
%----------------------------------------------------------------

plot_model_irfs;

%% APPLICATIONS

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

% sample

indic_early  = 0; % 1: early sample, 0 standard sample.

% Purely empirical, or use models?

indic_emp = 0; % 1 purely empirical, 0 use model-based IRFs

indic_RE    = 0; % only RE models?
indic_behav = 0; % only behavioral models?
indic_joint = 1; % all models?

% target counterfactual rule

cnfctl_0y       = 0; % output gap targeting
cnfctl_0pi      = 0; % inflation targeting
cnfctl_0ib      = 0; % nominal rate peg
cnfctl_tylr     = 0; % Taylor rule
cnfctl_ngdp     = 0; % NGDP targeting
cnfctl_ibtarget = 0; % rate target
cnfctl_optpol   = 1; % dual mandate

%----------------------------------------------------------------
% Figures 6-7
%----------------------------------------------------------------
 
indic_emp = 1;

get_cnfctl_stats;

%----------------------------------------------------------------
% Figure 8
%----------------------------------------------------------------

indic_emp = 0;

get_cnfctl_mbc;

%----------------------------------------------------------------
% Figure 9
%----------------------------------------------------------------

get_historical_evol;

%----------------------------------------------------------------
% Figure 10 and Figure D.3
%----------------------------------------------------------------

% RE

indic_RE    = 1; % only RE models?
indic_behav = 0; % only behavioral models?
indic_joint = 0; % all models?

get_historical_scenario;

% behavioral

indic_RE    = 0; % only RE models?
indic_behav = 1; % only behavioral models?
indic_joint = 0; % all models?

get_historical_scenario;

% joint

indic_RE     = 0; % only RE models?
indic_behav  = 0; % only behavioral models?
indic_joint  = 1; % all models?
indic_models = 1; % show individual models?

get_historical_scenario;

%----------------------------------------------------------------
% Figure 11
%----------------------------------------------------------------

decompose_realrates_brank;

%----------------------------------------------------------------
% Figure D.1 
%----------------------------------------------------------------

% second moments, early sample

indic_emp   = 0;
indic_early = 1;

get_cnfctl_stats;

%----------------------------------------------------------------
% Figure D.2 
%----------------------------------------------------------------

% historical evolution, only empirical IRFs

indic_emp   = 1;
indic_early = 0;

get_historical_evol;

%----------------------------------------------------------------
% Table 4.1
%----------------------------------------------------------------

get_posterior_probs;

%----------------------------------------------------------------
% Table C.2
%----------------------------------------------------------------

get_param_post;

%----------------------------------------------------------------
% Table D.1
%----------------------------------------------------------------

run_var_spf_fcst_compare;

%----------------------------------------------------------------
% Table D.2
%----------------------------------------------------------------

run_var_swfactors;