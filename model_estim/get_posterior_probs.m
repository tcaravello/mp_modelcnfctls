%% COMPUTE THE POSTERIOR PROBABILITIES FOR EACH MODEL, BOTH BEHAVIORAL AND NON-BEHAVIORAL.
% Tomas Caravello, Alisdair McKay and Christian Wolf

% This code uses the formula from Geweke (1999) as explained in Herbst and
% Schorfheide (2016) ch. 4.6

%% HOUSEKEEPING
 
task = '/model_estim';

models_use = [1; 1];
model_code = {'rank','hank'};

addpath([path vintage task '/_results'])

cd([path vintage task  ]);

compute_model_draws = 0;

tau = 0.5; %truncation: only use the best tau fraction of the draws.

%% PRIOR MODEL PROBABILITIES

% among non-behavioral models
priors_non_behav = [1/2;1/2];
% among behavioral models
priors_behav = [1/2;1/2];
% among all models

behav_total_prior = 1/2;
priors_all_models = [(1-behav_total_prior)*priors_non_behav; behav_total_prior * priors_behav];

%% LOAD PARAMETER DRAWS 

n_models_total = length(models_use);

% ratex
param_all = cell(n_models_total,1);
log_posterior_all = cell(n_models_total,1);
settings_all = cell(n_models_total,1);

% behav
param_all_behav = cell(n_models_total,1);
log_posterior_all_behav = cell(n_models_total,1);
settings_all_behav = cell(n_models_total,1);

for i_model = 1:n_models_total

if models_use(i_model) == 1

%ratex model
eval(['load',' ',model_code{i_model},...
    '_draws_other param_collector log_posterior_collector settings_MH'])    
param_all{i_model}= param_collector(:,settings_MH.N_burn+1:end);
log_posterior_all{i_model} = log_posterior_collector(:,settings_MH.N_burn+1:end);
settings_all{i_model} = settings_MH;

clear param_collector log_posterior_collector settings_MH
% behavioral model 
eval(['load',' ',model_code{i_model},...
    '_draws_other_behav param_collector log_posterior_collector settings_MH']) 
param_all_behav{i_model} = param_collector(:,settings_MH.N_burn+1:end);
log_posterior_all_behav{i_model} = log_posterior_collector(:,settings_MH.N_burn+1:end);
settings_all_behav{i_model} = settings_MH;
clear param_collector log_posterior_collector settings_MH
 
end
end

cd([path vintage task]);

global T_use

T_use = settings_all{1}.T_use;
cov_mat = 0; % 0 is non-diagonal

%% GET EMPIRICAL ESTIMATES

% These are needed because we need to add the log determinant to the
% posterior.

get_empirical_targets;

%% GET PRIORS

% Priors are needed to evaluate the posterior.

get_priors_general;

%% CONSTRUCT GEWEKE'S FUNCTION

marginal_lik_all  = zeros(n_models_total,1);
quantiles_use_all = cell(n_models_total,1);

marginal_lik_all_behav  = zeros(n_models_total,1);
quantiles_use_all_behav = cell(n_models_total,1);

for i_model = 1:n_models_total
if models_use(i_model) == 1
var_init = diag(settings_all{i_model}.var_mat_draws_init);
val_init = settings_all{i_model}; 
index_use = find(var_init~=0);

settings_all{i_model}.index_use  = index_use;
settings_all{i_model}.tau = tau;
settings_all{i_model}.target_Sigma_inv = target_Sigma_inv;
settings_all{i_model}.N_IRFS = length(Psi_i_mat_mean);
settings_all{i_model}.N_IRFS = length(Psi_i_mat_mean);


[marginal_lik_all(i_model), quantiles_use_all{i_model}] = ...
    harmonic_mean_geweke(param_all{i_model}, log_posterior_all{i_model}, settings_all{i_model},model_code{i_model});

[marginal_lik_all_behav(i_model), quantiles_use_all_behav{i_model}] = ...
    harmonic_mean_geweke(param_all_behav{i_model}, log_posterior_all_behav{i_model}, settings_all{i_model},model_code{i_model});

end
end

% Parameter means
param_mean_all = cell(n_models_total,1);
param_mean_all_behav = cell(n_models_total,1);

for i_model = 1:n_models_total
if models_use(i_model) == 1
param_mean_all{i_model} = mean(param_all{i_model},2);
param_mean_all_behav{i_model} = mean(param_all_behav{i_model},2);

end
end

%% CREATE TABLE

fprintf(' \n \n \n')
fprintf('\\begin{table}[tp!]  \n')
fprintf('     \\centering \n')
fprintf('     \\bigskip \n')
fprintf('\\begin{tabular}{c c cc c c} \\toprule \\toprule \n')
fprintf('Model & \\phantom{xxxx} & Baseline & Behavioral & \\phantom{xxxx} & Total \\\\ \n')
fprintf('     \\midrule \n')
for i_model = 1:n_models_total
if models_use(i_model) == 1
    numbs_aux = zeros(3,length(model_posteriors_all_models));
    numbs_aux(1,i_model) = 1;
    numbs_aux(2,i_model+n_models_total) = 1;
    numbs_aux(3,i_model) = 1;
    numbs_aux(3,i_model+n_models_total) = 1;
fprintf('%s && %1.4f & %1.4f && %1.4f \\\\ \n', ...
    upper(model_code{i_model}),numbs_aux*model_posteriors_all_models);
end
end
fprintf('& Total   && %1.4f  & %1.4f && %1.4f \\\\ \\bottomrule \n', [sum(model_posteriors_all_models(1:n_models_total)), sum(model_posteriors_all_models(n_models_total+1:end)), sum(model_posteriors_all_models(:))]);
fprintf('\\end{tabular} \n')
fprintf('\\caption{Posterior probabilities across the four models, assuming a uniform prior. The posterior model probabilities are computed as in \\eqref{eq:model_probs}.} \n')
fprintf('\\label{tab:mode_post} \n')
fprintf('\\end{table} \n')

% models for both versions

model_posteriors_behav_non_behav = cell(n_models_total,1);
for i_model = 1:n_models_total
if models_use(i_model) == 1
    model_posteriors_behav_non_behav{i_model} = model_posteriors_all_models([i_model (i_model+n_models_total)])/sum(model_posteriors_all_models([i_model (i_model+n_models_total)]));
end
end

%% CREATE DRAWS FROM MODELS

if compute_model_draws == 1

M_draws = settings_all{1}.N_keep;
% create uniform draws
unif_all_models = rand(M_draws,1);
% create discrete unform draws for each model. This is better than
% taking draws in order due to autocorrelation.
model_uni_draws = unidrnd(M_draws,M_draws,1); 

%create collectors
T_save = settings_all{1}.T_save;
T_use  = settings_all{1}.T_use;
N_keep = settings_all{1}.N_keep;
N_models = length(model_posteriors_all_models);

Pi_m_col = zeros([T_save T_save N_keep N_models]);
Y_m_col = zeros([T_save T_save N_keep N_models]);
R_n_m_col = zeros([T_save T_save N_keep N_models]);
m_fit_col = zeros([T_use 2 N_keep N_models]);

% non-behavioral 

% rank
load rank_draws_main
Pi_m_col(:,:,:,1) = Pi_m_collector;
R_n_m_col(:,:,:,1) = R_n_m_collector;
Y_m_col(:,:,:,1) = Y_m_collector;
m_fit_col(:,:,:,1) = m_fit_collector;

clear Pi_m_collector R_n_m_collector Y_m_collector m_fit_collector

% hank
load hank_draws_main
Pi_m_col(:,:,:,2) = Pi_m_collector;
R_n_m_col(:,:,:,2) = R_n_m_collector;
Y_m_col(:,:,:,2) = Y_m_collector;
m_fit_col(:,:,:,2) = m_fit_collector;

clear Pi_m_collector R_n_m_collector Y_m_collector m_fit_collector

% behavioral 

% rank
load rank_draws_main_behav
Pi_m_col(:,:,:,3) = Pi_m_collector;
R_n_m_col(:,:,:,3) = R_n_m_collector;
Y_m_col(:,:,:,3) = Y_m_collector;
m_fit_col(:,:,:,3) = m_fit_collector;

clear Pi_m_collector R_n_m_collector Y_m_collector m_fit_collector

% HANK
load hank_draws_main_behav
Pi_m_col(:,:,:,4) = Pi_m_collector;
R_n_m_col(:,:,:,4) = R_n_m_collector;
Y_m_col(:,:,:,4) = Y_m_collector;
m_fit_col(:,:,:,4) = m_fit_collector;

clear Pi_m_collector R_n_m_collector Y_m_collector m_fit_collector

%% SAMPLING

% non-behavioral

[Pi_m_collector_non_behav, Y_m_collector_non_behav, R_n_m_collector_non_behav, m_fit_collector_non_behav] = ...
    sample_from_models(M_draws, model_posteriors_non_behav, Pi_m_col(:,:,:,[1 2]), Y_m_col(:,:,:,[1 2]), R_n_m_col(:,:,:,[1 2]),...
     m_fit_col(:,:,:,[1 2]));

% behavioral 

[Pi_m_collector_behav, Y_m_collector_behav, R_n_m_collector_behav, m_fit_collector_behav] = ...
    sample_from_models(M_draws, model_posteriors_behav, Pi_m_col(:,:,:,[3 4]), Y_m_col(:,:,:,[3 4]), R_n_m_col(:,:,:,[3 4]),...
     m_fit_col(:,:,:,[3 4]));

% only RANK

[Pi_m_collector_only_rank, Y_m_collector_only_rank, R_n_m_collector_only_rank, m_fit_collector_only_rank] = ...
    sample_from_models(M_draws, model_posteriors_rank, Pi_m_col(:,:,:,[1 3]), Y_m_col(:,:,:,[1 3]), R_n_m_col(:,:,:,[1 3]),...
     m_fit_col(:,:,:,[1 3]));

% only HANK

[Pi_m_collector_only_hank, Y_m_collector_only_hank, R_n_m_collector_only_hank, m_fit_collector_only_hank] = ...
    sample_from_models(M_draws, model_posteriors_hank, Pi_m_col(:,:,:,[2 4]), Y_m_col(:,:,:,[2 4]), R_n_m_col(:,:,:,[2 4]),...
     m_fit_col(:,:,:,[2 4]));

% all models

[Pi_m_collector, Y_m_collector, R_n_m_collector, m_fit_collector] = ...
    sample_from_models(M_draws, model_posteriors_all_models, Pi_m_col, Y_m_col, R_n_m_col,...
     m_fit_col);

%% SAVE RESULTS

if save_results == 1

aux_path = [path vintage '/suff_stats'];

cd([aux_path '/joint'])

save all_models_draws  Pi_m_collector Y_m_collector R_n_m_collector m_fit_collector model_posteriors_all_models
save rank_behav_and_non_behav_draws Pi_m_collector_only_rank Y_m_collector_only_rank R_n_m_collector_only_rank m_fit_collector_only_rank model_posteriors_rank
save hank_behav_and_non_behav_draws Pi_m_collector_only_hank Y_m_collector_only_hank R_n_m_collector_only_hank m_fit_collector_only_hank model_posteriors_hank


cd([aux_path '/ratex'])
save non_behav_models_draws  Pi_m_collector_non_behav Y_m_collector_non_behav R_n_m_collector_non_behav m_fit_collector_non_behav model_posteriors_non_behav

cd([aux_path '/behavioral_fixed'])
save behav_all_models_draws  Pi_m_collector_behav Y_m_collector_behav R_n_m_collector_behav m_fit_collector_behav model_posteriors_behav
end

end