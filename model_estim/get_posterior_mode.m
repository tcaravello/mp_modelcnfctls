%% FIND POSTERIOR MODE AND GET IRFS FOR ANY OF THE FOUR STRUCTURAL MODELS.
% Tomas Caravello, Alisdair McKay and Christian Wolf

%% HOUSEKEEPING
 
clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

path = '/Users/tomyc/Dropbox (MIT)/mp_modelcnfctls/code/github_public/mp_modelcnfctls';
vintage = '';
task = '/model_estim';

model = '/hank'; %type '/rank' for RANK and'/hank' for HANK, 

behavioral = 1; %1 if you want to estimate the behavioral parameters, 0 if you wanna do ratex.

%% SETTINGS

global T T_use n_shock save_fig shocks_match

shocks_match = 2; %1 if you want to match only the AD shock, 2 if you want to match both the AD and RR shocks. (Using 1 will give an error).

save_results = 0; %1 if you want to save the results.

save_fig = 1; %1 if you want to save the IRF plots.

% only applies if model = hank. 1: load Jacobians from previous run.
% 0: compute steady state and jacobians from scratch, takes a few seconds.
load_hank = 1;

%1 if you want to use diagonal, 0 if you want to use non-diagonal.
% in the paper we use non-diagonal because it better represents the
% informativeness of the data to do model selection.
cov_mat = 0; 

n_shock = 25;       % Number of news shocks to use
T_use = 25;         % Horizon to match IRFs
T = T_use + 275;    % Truncation horizon for model solution
T_save = 200;        % Horizon used to save IRF matrices.


addpath([path vintage '/_auxiliary_functions'])
addpath([path vintage task '/_results'])
addpath([path vintage '/var_inputs/_results'])

cd([path vintage task]);

%% EMPIRICAL TARGETS

get_empirical_targets;

%% CALIBRATION

get_calibration_general;

%% JACOBIANS

get_baseline_jacobians;

%% PRIORS

get_priors_general;

%% OBTAIN POSTERIOR MODE

if behavioral == 1

if strcmp(model, '/rank')
     x0 = [0.71    0.86    0.99999    0.946    0.99999    5.3    0.47    m_d_mean     m_f_mean];
elseif strcmp(model, '/hank')
     x0 = [0.9598    0.85    0.99999    0.95     0.99999    5.3219    0.47    m_d_mean    m_f_mean];
end


lower_bound = [0,0,0,0,0,0,0,0,0]; 
%upper_bound = [1,1,1,1,1,inf,1,1,1];
upper_bound = [1,1,1,1,1,1,1,1,1]-(1e-13);
upper_bound(6) = inf;

else
    if strcmp(model, '/rank')
        x0 = [0.73    0.85   0.99999    0.74   0.99999    5.3    0.47];
    elseif strcmp(model, '/hank')
        x0 = [0.9654    0.9530    0.999999    0.8763    0.999999    5.7686    0.4618];
    end
lower_bound = [0,0,0,0,0,0,0]; 
%upper_bound = [1,1,1,1,1,inf,1];
upper_bound = [1,1,1,1,1,1,1]-(1e-13);
upper_bound(6) = inf;

end


% use low levels of tolerance for illustration purpose only.

options_min = optimset('Display','iter','MaxFunEvals',1000,'tolFun',1e-4, 'TolX', 1e-4);

% get posterior mode
[param_sol,posterior_mode,exit_flag,output] = fminsearchcon(@(x) solve_model(x,target,target_Sigma_inv,model),x0, lower_bound,upper_bound, [],[],[],options_min);

if strcmp(model, '/rank')

fprintf([' h: %1.4f \n calvo_p: %1.5f \n zeta_p: %1.3f \n' ...
    'calvo_w: %1.5f \n zeta_w: %1.3f \n' ...
    ' inv adj cost: %2.2f \n psi_uti: %1.3f \n m_d %1.4f \n m_f %1.4f \n'],param_sol)
else

fprintf([' theta: %1.4f \n calvo_p: %1.5f \n zeta_p: %1.3f \n' ...
    'calvo_w: %1.5f \n zeta_w: %1.3f \n' ...
    ' inv adj cost: %2.2f \n psi_uti: %1.3f \n m_d %1.4f \n m_f %1.4f '],param_sol)
end

kappa_p = (1-beta*param_sol(2))*(1-param_sol(2))/param_sol(2);
kappa_w = (1-beta*param_sol(4))*(1-param_sol(4))/param_sol(4);

if length(param_sol)>7
beta_p = beta*param_sol(2)*param_sol(9)*(1+kappa_p/(1-beta*param_sol(2)*param_sol(9)));
beta_w = beta*param_sol(4)*param_sol(9)*(1+kappa_w/(1-beta*param_sol(4)*param_sol(9)));
else
beta_p = beta;
beta_w = beta;
end
fprintf([' \n kappa_p %1.4f \n kappa_w: %1.4f \n beta_p: %1.4f \n' ...
    'beta_w: %1.4f'],[kappa_p kappa_w beta_p beta_w])

%% FINAL IRFs

[error_fit, Pi_m_model, Y_m_model, R_n_m_model, m_fit_model,...
    Pi_m_fit, Y_m_fit, R_n_m_fit, C_m_model, Inv_m_model, C_m_fit, Inv_m_fit] = solve_model(param_sol,target,target_Sigma_inv,model);
 
if shocks_match ==2 

    Pi_m_fit = zeros(T_use,2);
    Y_m_fit  = zeros(T_use,2);
    R_n_m_fit  = zeros(T_use,2);

    for j=1:2
    Pi_m_fit(:,j) = Pi_m_model(1:T_use,1:T_use)*m_fit_model(:,j);
    Y_m_fit(:,j)  = Y_m_model(1:T_use,1:T_use)*m_fit_model(:,j);
    R_n_m_fit(:,j)  = R_n_m_model(1:T_use,1:T_use)*m_fit_model(:,j);
    end

else

m_fit_model_2 = [m_fit_model(1:end-1);0]; 
Pi_m_fit_2 = Pi_m_model(1:T_use,1:T_use) * m_fit_model_2;
Y_m_fit_2 = Y_m_model(1:T_use,1:T_use) * m_fit_model_2;
R_n_m_fit_2 = R_n_m_model(1:T_use,1:T_use) * m_fit_model_2;

end

%% SAVE POSTERIOR MODE

if save_results == 1

if behavioral == 1

    if strcmp(model, '/rank')

Pi_m_rank = Pi_m_model;
Y_m_rank = Y_m_model;
R_n_m_rank = R_n_m_model; 
m_fit_rank = m_fit_model;

Pi_m_base_rank = Pi_m_rank(1:T_save,1:T_save); 
Y_m_base_rank = Y_m_rank(1:T_save,1:T_save); 
I_m_base_rank = R_n_m_rank(1:T_save,1:T_save);

cd([path vintage task '/_results']);

if shocks_match == 2
    save params_rank_ad_behav_2shocks T T_use target_Sigma_inv param_sol posterior_mode Pi_m_rank Y_m_rank R_n_m_rank m_fit_rank
    save rank_base_ad_behav_2shocks Pi_m_base_rank Y_m_base_rank I_m_base_rank
else
save params_rank_ad_behav T T_use target_Sigma_inv param_sol posterior_mode Pi_m_rank Y_m_rank R_n_m_rank m_fit_rank
save rank_base_ad_behav Pi_m_base_rank Y_m_base_rank I_m_base_rank
end

cd([path vintage task]);

elseif strcmp(model, '/hank')

Pi_m_hank = Pi_m_model;
Y_m_hank = Y_m_model;
R_n_m_hank = R_n_m_model; 
m_fit_hank = m_fit_model;

Pi_m_base_hank = Pi_m_hank(1:T_save,1:T_save); 
Y_m_base_hank = Y_m_hank(1:T_save,1:T_save); 
I_m_base_hank = R_n_m_hank(1:T_save,1:T_save);

cd([path vintage task '/_results']);

if shocks_match == 2
    save params_hank_ad_behav_2shocks T T_use target_Sigma_inv param_sol posterior_mode Pi_m_hank Y_m_hank R_n_m_hank m_fit_hank
    save hank_base_ad_behav_2shocks Pi_m_base_hank Y_m_base_hank I_m_base_hank
else
save params_hank_ad_behav T T_use target_Sigma_inv param_sol posterior_mode Pi_m_hank Y_m_hank R_n_m_hank m_fit_hank
save hank_base_ad_behav Pi_m_base_hank Y_m_base_hank I_m_base_hank
end

cd([path vintage task]);

end
else %no behavioral

if strcmp(model, '/rank')

Pi_m_rank = Pi_m_model;
Y_m_rank = Y_m_model;
R_n_m_rank = R_n_m_model; 
m_fit_rank = m_fit_model;

Pi_m_base_rank = Pi_m_rank(1:T_save,1:T_save); 
Y_m_base_rank = Y_m_rank(1:T_save,1:T_save); 
I_m_base_rank = R_n_m_rank(1:T_save,1:T_save);

cd([path vintage task '/_results']);

if shocks_match == 2
    save params_rank_ad_2shocks T T_use target_Sigma_inv param_sol posterior_mode Pi_m_rank Y_m_rank R_n_m_rank m_fit_rank
    save rank_base_ad_2shocks Pi_m_base_rank Y_m_base_rank I_m_base_rank
else
save params_rank_ad T T_use target_Sigma_inv param_sol posterior_mode Pi_m_rank Y_m_rank R_n_m_rank m_fit_rank
save rank_base_ad Pi_m_base_rank Y_m_base_rank I_m_base_rank
end

cd([path vintage task]);

elseif strcmp(model, '/hank')

Pi_m_hank = Pi_m_model;
Y_m_hank = Y_m_model;
R_n_m_hank = R_n_m_model; 
m_fit_hank = m_fit_model;

Pi_m_base_hank = Pi_m_hank(1:T_save,1:T_save); 
Y_m_base_hank = Y_m_hank(1:T_save,1:T_save); 
I_m_base_hank = R_n_m_hank(1:T_save,1:T_save);

cd([path vintage task '/_results']);

if shocks_match == 2
    save params_hank_ad_2shocks T T_use target_Sigma_inv param_sol posterior_mode Pi_m_hank Y_m_hank R_n_m_hank m_fit_hank
    save hank_base_ad_2shocks Pi_m_base_hank Y_m_base_hank I_m_base_hank
else
    save params_hank_ad T T_use target_Sigma_inv param_sol posterior_mode Pi_m_hank Y_m_hank R_n_m_hank m_fit_hank
    save hank_base_ad Pi_m_base_hank Y_m_base_hank I_m_base_hank
end

cd([path vintage task]);

end
end
end

%% PLOT RESULTS

%----------------------------------------------------------------
% General Settings
%----------------------------------------------------------------

% plot horizon

IRF_hor_plot = T_use;

% results scaling

scale = 1;

% go to results folder

cd([path vintage task '/_results']);

if strcmp(model, '/rank')
    if behavioral == 0
        model_legend = 'RANK';
    elseif behavioral == 1
        model_legend = 'B-RANK';
    end
elseif strcmp(model, '/hank')
    if behavioral == 0
        model_legend = 'HANK';
    elseif behavioral == 1
        model_legend = 'B-HANK';
    end
end

%----------------------------------------------------------------
% Figure
%----------------------------------------------------------------

settings.colors.black  = [0 0 0];
settings.colors.grey   = [230/255 230/255 230/255];
settings.colors.dgrey  = [130/255 130/255 130/255];
settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.lblue  = 0.25 * settings.colors.blue + 0.75 * [1 1 1];

plotwidth = 0.27;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * gapsize + 2 * plotwidth];

if shocks_match == 2
shocks_legend    = {'AD','RR'};
% AD shock

for shock_use = 1:2
figure(shock_use)
set(gcf, 'units','normalized','outerposition',[0.0 0.0 0.9 0.9]);

subplot(1,3,1)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
jbfill(1:1:IRF_hor_plot,scale*squeeze(Y_m_lb_emp(1:IRF_hor_plot,1,shock_use))',...
       scale*squeeze(Y_m_ub_emp(1:IRF_hor_plot,1,shock_use))',settings.colors.grey,settings.colors.grey,0,1);
hold on
plot(1:1:IRF_hor_plot,scale*squeeze(Y_m_emp(1:IRF_hor_plot,1,shock_use)),':','Color',settings.colors.dgrey,'LineWidth',4)
hold on
plot(1:1:IRF_hor_plot,scale*Y_m_fit(1:IRF_hor_plot,shock_use),'Color',settings.colors.blue,'LineWidth',4)
hold on
xlim([1 IRF_hor_plot])
%ylim([-2 2])
%yticks([-2 -1 0 1 2])
set(gcf,'color','w')
title('Output','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
ylabel('\% Deviation','interpreter','latex','FontSize',20)
grid on
hold off

subplot(1,3,2)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
jbfill(1:1:IRF_hor_plot,scale*squeeze(Pi_m_lb_emp(1:IRF_hor_plot,1,shock_use))',...
       scale*squeeze(Pi_m_ub_emp(1:IRF_hor_plot,1,shock_use))',settings.colors.grey,settings.colors.grey,0,1);
hold on
plot(1:1:IRF_hor_plot,scale*squeeze(Pi_m_emp(1:IRF_hor_plot,1,shock_use)),':','Color',settings.colors.dgrey,'LineWidth',4)
hold on
plot(1:1:IRF_hor_plot,scale*4*Pi_m_fit(1:IRF_hor_plot,shock_use),'Color',settings.colors.blue,'LineWidth',4)
hold on
xlim([1 IRF_hor_plot])
set(gcf,'color','w')
title('Inflation','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
grid on
hold off

subplot(1,3,3)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(3);
pos(3) = plotwidth;
set(gca,'Position', pos)
jbfill(1:1:IRF_hor_plot,scale*squeeze(R_n_m_lb_emp(1:IRF_hor_plot,1,shock_use))',...
       scale*squeeze(R_n_m_ub_emp(1:IRF_hor_plot,1,shock_use))',settings.colors.grey,settings.colors.grey,0,1);
hold on
plot(1:1:IRF_hor_plot,scale*squeeze(R_n_m_emp(1:IRF_hor_plot,1,shock_use)),':','Color',settings.colors.dgrey,'LineWidth',4)
hold on
plot(1:1:IRF_hor_plot,scale*4*R_n_m_fit(1:IRF_hor_plot,shock_use),'Color',settings.colors.blue,'LineWidth',4)
hold on
xlim([1 IRF_hor_plot])
set(gcf,'color','w')
title('Interest Rate','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
legend({'',['Empirical -' shocks_legend{shock_use}],model_legend},'Location','Southeast','fontsize',18,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.25*pos(3) 1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
if save_fig == 1
    print([model_legend '_IRFs_implementation_' num2str(n_shock) '_2shocks_' shocks_legend{shock_use}], '-dpng')
end

end
else

close all

figure(1)
set(gcf, 'units','normalized','outerposition',[0.0 0.0 0.9 0.9]);

subplot(1,3,1)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
jbfill(1:1:IRF_hor_plot,scale*Y_m_lb_emp(1:IRF_hor_plot)',...
       scale*Y_m_ub_emp(1:IRF_hor_plot)',settings.colors.grey,settings.colors.grey,0,1);
hold on
plot(1:1:IRF_hor_plot,scale*Y_m_emp(1:IRF_hor_plot),':','Color',settings.colors.dgrey,'LineWidth',4)
hold on
plot(1:1:IRF_hor_plot,scale*Y_m_fit(1:IRF_hor_plot),'Color',settings.colors.blue,'LineWidth',4)
hold on
xlim([1 IRF_hor_plot])
%ylim([-2 2])
%yticks([-2 -1 0 1 2])
set(gcf,'color','w')
title('Output','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
ylabel('\% Deviation','interpreter','latex','FontSize',20)
grid on
hold off

subplot(1,3,2)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
jbfill(1:1:IRF_hor_plot,scale*Pi_m_lb_emp(1:IRF_hor_plot)',...
       scale*Pi_m_ub_emp(1:IRF_hor_plot)',settings.colors.grey,settings.colors.grey,0,1);
hold on
plot(1:1:IRF_hor_plot,scale*Pi_m_emp(1:IRF_hor_plot),':','Color',settings.colors.dgrey,'LineWidth',4)
hold on
plot(1:1:IRF_hor_plot,scale*4*Pi_m_fit(1:IRF_hor_plot),'Color',settings.colors.blue,'LineWidth',4)
hold on
xlim([1 IRF_hor_plot])
set(gcf,'color','w')
title('Inflation','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
grid on
hold off

subplot(1,3,3)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(3);
pos(3) = plotwidth;
set(gca,'Position', pos)
jbfill(1:1:IRF_hor_plot,scale*R_n_m_lb_emp(1:IRF_hor_plot)',...
       scale*R_n_m_ub_emp(1:IRF_hor_plot)',settings.colors.grey,settings.colors.grey,0,1);
hold on
plot(1:1:IRF_hor_plot,scale*R_n_m_emp(1:IRF_hor_plot),':','Color',settings.colors.dgrey,'LineWidth',4)
hold on
plot(1:1:IRF_hor_plot,scale*4*R_n_m_fit(1:IRF_hor_plot),'Color',settings.colors.blue,'LineWidth',4)
hold on
xlim([1 IRF_hor_plot])
set(gcf,'color','w')
title('Interest Rate','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
legend({'','Empirical',model_legend},'Location','Southeast','fontsize',18,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.25*pos(3) 1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
if save_fig == 1
    print([model_legend '_IRFs_implementation_' num2str(n_shock)], '-dpng')
end

end

cd([path vintage task]);