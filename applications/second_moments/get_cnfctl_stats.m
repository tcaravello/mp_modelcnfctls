%% UNCONDITIONAL SECOND MOMENTS UNDER COUNTERFACTUAL POLICY RULE
% Tomas Caravello, Alisdair McKay, Christian Wolf

%% HOUSEKEEPING

experiment = '/applications/second_moments';

save_fig = 1;

addpath([path vintage '/_auxiliary_functions'])
addpath([path vintage '/var_inputs/_results']);
cd([path vintage experiment]);

%% IMPORTS & SETTINGS

%----------------------------------------------------------------
% Experiment
%----------------------------------------------------------------

% show individual models?
 
indic_models = 1;

%----------------------------------------------------------------
% Policy Shock Sufficient Statistics
%----------------------------------------------------------------

% import

import_suffstats

% sizes

T       = size(Pi_m_draws,1);
n_draws = size(Pi_m_draws,3);

shock_max = T;

%----------------------------------------------------------------
% Wold IRFs
%----------------------------------------------------------------

if indic_early == 0

    load wold_results
    
elseif indic_early == 1
    
    load wold_results_early
    
end

wold_base = IS_wold.Theta_OLS([9 2 10],:,1:T); % variables: (pi, y, i)

series_names = series_names([9 2 10]);

n_y      = size(wold_base,1);
n_shocks = size(wold_base,2);

clear IS_wold

%----------------------------------------------------------------
% Specify Counterfactual Rule
%----------------------------------------------------------------

cnfctl_0y       = 0; % output gap targeting
cnfctl_0pi      = 0; % inflation targeting
cnfctl_0ib      = 0; % nominal rate peg
cnfctl_tylr     = 0; % Taylor rule
cnfctl_ngdp     = 0; % NGDP targeting
cnfctl_ibtarget = 0; % rate target
cnfctl_optpol   = 1; % optimal dual mandate

set_cnfctl_rule

%% CONSTRUCT COUNTERFACTUAL WOLD IRFs

%----------------------------------------------------------------
% Counterfactual Posterior Draws
%----------------------------------------------------------------

wold_cnfctl = NaN(n_y,n_shocks,T,n_draws);

for i_draw = 1:n_draws
    
% get policy shock IRFs

Pi_m = Pi_m_draws(:,1:shock_max,i_draw);
Y_m  = Y_m_draws(:,1:shock_max,i_draw);
I_m  = I_m_draws(:,1:shock_max,i_draw);

for i_shock = 1:n_shocks

% get the corresponding base sequences

pi_z = squeeze(wold_base(1,i_shock,:));
y_z  = squeeze(wold_base(2,i_shock,:));
i_z  = squeeze(wold_base(3,i_shock,:));

% find best fit to counterfactual rule

if cnfctl_optpol == 0

[pi_z_cnfctl,y_z_cnfctl,i_z_cnfctl] = cnfctl_fn(A_pi,A_y,A_i,wedge,...
    Pi_m,Y_m,I_m,pi_z,y_z,i_z);

elseif cnfctl_optpol == 1

[pi_z_cnfctl,y_z_cnfctl,i_z_cnfctl] = optpol_fn(W_pi,W_y,W_i,...
    Pi_m,Y_m,I_m,pi_z,y_z,i_z);

end

wold_cnfctl(1,i_shock,:,i_draw) = pi_z_cnfctl;
wold_cnfctl(2,i_shock,:,i_draw) = y_z_cnfctl;
wold_cnfctl(3,i_shock,:,i_draw) = i_z_cnfctl;

end

end

%----------------------------------------------------------------
% Individual Models
%----------------------------------------------------------------

if indic_models == 1

n_models           = 4;
wold_cnfctl_models = NaN(n_y,n_shocks,T,n_models,n_draws);

for i_model = 1:n_models
    
for i_draw = 1:n_draws
    
% get policy shock IRFs

if i_model == 1

    Pi_m = Pi_m_rank_draws(:,1:shock_max,i_draw);
    Y_m  = Y_m_rank_draws(:,1:shock_max,i_draw);
    I_m  = I_m_rank_draws(:,1:shock_max,i_draw);

elseif i_model == 2

    Pi_m = Pi_m_hank_draws(:,1:shock_max,i_draw);
    Y_m  = Y_m_hank_draws(:,1:shock_max,i_draw);
    I_m  = I_m_hank_draws(:,1:shock_max,i_draw);

elseif i_model == 3

    Pi_m = Pi_m_brank_draws(:,1:shock_max,i_draw);
    Y_m  = Y_m_brank_draws(:,1:shock_max,i_draw);
    I_m  = I_m_brank_draws(:,1:shock_max,i_draw);

elseif i_model == 4

    Pi_m = Pi_m_bhank_draws(:,1:shock_max,i_draw);
    Y_m  = Y_m_bhank_draws(:,1:shock_max,i_draw);
    I_m  = I_m_bhank_draws(:,1:shock_max,i_draw);

end

for i_shock = 1:n_shocks

% get the corresponding base sequences

pi_z = squeeze(wold_base(1,i_shock,:));
y_z  = squeeze(wold_base(2,i_shock,:));
i_z  = squeeze(wold_base(3,i_shock,:));

% find best fit to counterfactual rule

if cnfctl_optpol == 0

[pi_z_cnfctl,y_z_cnfctl,i_z_cnfctl] = cnfctl_fn(A_pi,A_y,A_i,wedge,...
    Pi_m,Y_m,I_m,pi_z,y_z,i_z);

elseif cnfctl_optpol == 1

[pi_z_cnfctl,y_z_cnfctl,i_z_cnfctl] = optpol_fn(W_pi,W_y,W_i,...
    Pi_m,Y_m,I_m,pi_z,y_z,i_z);

end

wold_cnfctl_models(1,i_shock,:,i_model,i_draw) = pi_z_cnfctl;
wold_cnfctl_models(2,i_shock,:,i_model,i_draw) = y_z_cnfctl;
wold_cnfctl_models(3,i_shock,:,i_model,i_draw) = i_z_cnfctl;

end

end

end

end

%----------------------------------------------------------------
% Empirics
%----------------------------------------------------------------

if indic_emp == 1

% imports

clear Pi_m_draws
clear Y_m_draws
clear I_m_draws

import_suffstats_emp

% sizes

T       = size(Pi_m_draws,1);
n_draws = size(Pi_m_draws,3);

% shock space

shock_max = size(Pi_m_draws,2);

% counterfactuals

wold_cnfctl_emp = NaN(n_y,n_shocks,T,n_draws);

for i_draw = 1:n_draws
    
% get policy shock IRFs

Pi_m = Pi_m_draws(:,1:shock_max,i_draw);
Y_m  = Y_m_draws(:,1:shock_max,i_draw);
I_m  = I_m_draws(:,1:shock_max,i_draw);

for i_shock = 1:n_shocks

% get the corresponding base sequences

pi_z = squeeze(wold_base(1,i_shock,:));
y_z  = squeeze(wold_base(2,i_shock,:));
i_z  = squeeze(wold_base(3,i_shock,:));

% find best fit to counterfactual rule

if cnfctl_optpol == 0

[pi_z_cnfctl,y_z_cnfctl,i_z_cnfctl] = cnfctl_fn(A_pi,A_y,A_i,wedge,...
    Pi_m,Y_m,I_m,pi_z,y_z,i_z);

elseif cnfctl_optpol == 1

[pi_z_cnfctl,y_z_cnfctl,i_z_cnfctl] = optpol_fn(W_pi,W_y,W_i,...
    Pi_m,Y_m,I_m,pi_z,y_z,i_z);

end

wold_cnfctl_emp(1,i_shock,:,i_draw) = pi_z_cnfctl;
wold_cnfctl_emp(2,i_shock,:,i_draw) = y_z_cnfctl;
wold_cnfctl_emp(3,i_shock,:,i_draw) = i_z_cnfctl;

end

end

end

%% GET COUNTERFACTUAL BUSINESS-CYCLE STATISTICS

%----------------------------------------------------------------
% Baseline
%----------------------------------------------------------------

% VMA-implied variance-covariance matrix

base.cov = zeros(n_y,n_y);
for i_hor = 1:T
    base.cov = base.cov + wold_base(:,:,i_hor) * wold_base(:,:,i_hor)';
end

% correlations

base.corr = zeros(n_y,n_y);
for i_y = 1:n_y
    for ii_y = 1:n_y
        base.corr(i_y,ii_y) = base.cov(i_y,ii_y)/sqrt(base.cov(i_y,i_y) * base.cov(ii_y,ii_y));
    end
end

%----------------------------------------------------------------
% Counterfactual (draws)
%----------------------------------------------------------------

% VMA-implied variance-covariance matrix

cnfctl.cov = zeros(n_y,n_y,n_draws);
for i_draw = 1:n_draws
    for i_hor = 1:T
        cnfctl.cov(:,:,i_draw) = cnfctl.cov(:,:,i_draw) + wold_cnfctl(:,:,i_hor,i_draw) * wold_cnfctl(:,:,i_hor,i_draw)';
    end
end

cnfctl.cov_lb  = quantile(cnfctl.cov,0.16,3);
cnfctl.cov_med = quantile(cnfctl.cov,0.5,3);
cnfctl.cov_ub  = quantile(cnfctl.cov,0.84,3);

% correlations

cnfctl.corr = zeros(n_y,n_y,n_draws);
for i_draw = 1:n_draws
    for i_y = 1:n_y
        for ii_y = 1:n_y
            cnfctl.corr(i_y,ii_y,i_draw) = cnfctl.cov(i_y,ii_y,i_draw)/sqrt(cnfctl.cov(i_y,i_y,i_draw) * cnfctl.cov(ii_y,ii_y,i_draw));
        end
    end
end

cnfctl.corr_lb  = quantile(cnfctl.corr,0.16,3);
cnfctl.corr_med = quantile(cnfctl.corr,0.5,3);
cnfctl.corr_ub  = quantile(cnfctl.corr,0.84,3);

%----------------------------------------------------------------
% Counterfactual (models)
%----------------------------------------------------------------

if indic_models == 1

% VMA-implied variance-covariance matrix

cnfctl_models.cov = zeros(n_y,n_y,n_models,n_draws);
for i_model = 1:n_models
    for i_draw = 1:n_draws
        for i_hor = 1:T
            cnfctl_models.cov(:,:,i_model,i_draw) = cnfctl_models.cov(:,:,i_model,i_draw) ...
                + wold_cnfctl_models(:,:,i_hor,i_model,i_draw) * wold_cnfctl_models(:,:,i_hor,i_model,i_draw)';
        end
    end
end

cnfctl_models.cov = quantile(cnfctl_models.cov,0.5,4);

end

%----------------------------------------------------------------
% Counterfactual (empirics)
%----------------------------------------------------------------

if indic_emp == 1

% VMA-implied variance-covariance matrix

cnfctl_emp.cov = zeros(n_y,n_y,n_draws);
for i_draw = 1:n_draws
    for i_hor = 1:T
        cnfctl_emp.cov(:,:,i_draw) = cnfctl_emp.cov(:,:,i_draw) + wold_cnfctl_emp(:,:,i_hor,i_draw) * wold_cnfctl_emp(:,:,i_hor,i_draw)';
    end
end

cnfctl_emp.cov_lb  = quantile(cnfctl_emp.cov,0.16,3);
cnfctl_emp.cov_med = quantile(cnfctl_emp.cov,0.5,3);
cnfctl_emp.cov_ub  = quantile(cnfctl_emp.cov,0.84,3);

end

%% PLOT COUNTERFACTUALS

%----------------------------------------------------------------
% General Settings
%----------------------------------------------------------------

% plot horizon

IRF_hor_plot = 30;

% color settings

settings.colors.black  = [0 0 0];
settings.colors.dgrey  = [130/255 130/255 130/255];
settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.lblue  = 0.25 * settings.colors.blue + 0.75 * [1 1 1];
settings.colors.orange = [255/255 153/255 51/255];
settings.colors.lorange= 0.25 * settings.colors.orange + 0.75 * [1 1 1];
settings.colors.models = [196/255 174/255 120/255; ... % beige
                            204/255 0/255 0/255; ... % red
                            102/255 178/255 255/255]; % blue

% plot size

plotwidth = 0.27;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * gapsize + 2 * plotwidth];

% variable order

var_order = [2 1 3];

%----------------------------------------------------------------
% Posterior Standard Deviations
%----------------------------------------------------------------

if indic_RE == 1
    cd([path vintage experiment '/_results/re']);
elseif indic_behav == 1
    cd([path vintage experiment '/_results/behav']);
elseif indic_joint == 1
    cd([path vintage experiment '/_results/joint']);
end

n_kernel = 1001;
n_gap    = 20;

% baseline

figure

for ii_y = 1:n_y

i_y = var_order(ii_y);

grid_lb       = 0;
grid_ub       = 2 * max(sqrt(cnfctl.cov_ub(i_y,i_y)),sqrt(base.cov(i_y,i_y)));
grid_plot     = linspace(grid_lb,grid_ub,n_kernel)';
dist_plot     = kernel(grid_plot,winsorize(squeeze(sqrt(cnfctl.cov(i_y,i_y,:))),95),0.4);
dist_plot     = dist_plot ./ sum(dist_plot * grid_plot(2));

plot_lb_indx = find(dist_plot ~= 0, 1 );
plot_ub_indx = find(dist_plot ~= 0, 1, 'last' );

[~,base_pos] = min(abs(grid_plot - sqrt(base.cov(i_y,i_y))));

plot_lb_indx = max(min(plot_lb_indx - n_gap, base_pos - n_gap), 1);
plot_ub_indx = min(max(plot_ub_indx + n_gap, base_pos + n_gap), n_kernel);

subplot(1,3,ii_y)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(ii_y);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
plot(0,0,':','Color',settings.colors.black,'LineWidth',4)
hold on
jbfill(0,0,0,...
    settings.colors.lblue,settings.colors.blue,0,0.5);
hold on
if indic_models == 1
    for i_model = 1:n_models
        if i_model <= 2
            plot(0,0,'-','Color',settings.colors.models(i_model,:),'LineWidth',4)
        else
            plot(0,0,'-.','Color',settings.colors.models(i_model-2,:),'LineWidth',4)
        end
        hold on
    end
    hold on
end
plot(grid_plot,dist_plot,'Color',settings.colors.blue,'LineWidth',4)
hold on
jbfill(grid_plot,0*dist_plot,dist_plot,...
    settings.colors.lblue,settings.colors.lblue,0,0.5);
hold on
plot([sqrt(base.cov(i_y,i_y)) sqrt(base.cov(i_y,i_y))],[0 100],':','Color',settings.colors.black,'LineWidth',4)
hold on
if indic_models == 1
    for i_model = 1:n_models
        if i_model <= 2
            plot([sqrt(cnfctl_models.cov(i_y,i_y,i_model)) sqrt(cnfctl_models.cov(i_y,i_y,i_model))],[0 100],'-','Color',settings.colors.models(i_model,:),'LineWidth',4)
        else
            plot([sqrt(cnfctl_models.cov(i_y,i_y,i_model)) sqrt(cnfctl_models.cov(i_y,i_y,i_model))],[0 100],'-.','Color',settings.colors.models(i_model-2,:),'LineWidth',3)
        end
        hold on
    end
end
if ii_y < 3
    xlim([grid_plot(plot_lb_indx) grid_plot(plot_ub_indx)])
else
    xlim([1.95 grid_plot(plot_ub_indx)])
end
ylim([0 1.2 * max(dist_plot)])
% yticks([])
xlabel('St. Dev.','interpreter','latex','FontSize',20)
set(gcf,'color','w')
title(series_names(i_y),'interpreter','latex','fontsize',24)
grid on

if ii_y == 3
    if indic_models == 1
        legend({'Data','Counterfct''l','RANK','HANK','B-RANK','B-HANK'},'Location','Northeast','fontsize',18,'interpreter','latex',...
            'NumColumns',2,'Orientation','horizontal')
    else
        legend({'Data','Counterfct''l'},'Location','Northwest','fontsize',18,'interpreter','latex','NumColumns',2,'Orientation','horizontal')
    end
end

hold off

end

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.25*pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');

if save_fig == 1 

if indic_early == 0
    print('cnfctl_histograms_optpol','-dpng');
elseif indic_early == 1
    print('cnfctl_histograms_optpol_early','-dpng');
end
    
end

% with empirical shocks (density)

if indic_emp == 1

figure

for ii_y = 1:n_y

i_y = var_order(ii_y);

grid_lb_1       = 0;
grid_ub_1       = 2 * max(sqrt(cnfctl.cov_ub(i_y,i_y)),sqrt(base.cov(i_y,i_y)));
grid_plot_1     = linspace(grid_lb_1,grid_ub_1,n_kernel)';
dist_plot_1     = kernel(grid_plot_1,winsorize(squeeze(sqrt(cnfctl.cov(i_y,i_y,:))),95),0.4);
dist_plot_1     = dist_plot_1 ./ sum(dist_plot_1 * grid_plot_1(2));

plot_lb_indx_1 = find(dist_plot_1 ~= 0, 1 );
plot_ub_indx_1 = find(dist_plot_1 ~= 0, 1, 'last' );

[~,base_pos_1] = min(abs(grid_plot_1 - sqrt(base.cov(i_y,i_y))));

plot_lb_indx_1 = max(min(plot_lb_indx_1 - n_gap, base_pos_1 - n_gap), 1);
plot_ub_indx_1 = min(max(plot_ub_indx_1 + n_gap, base_pos_1 + n_gap), n_kernel);

grid_lb_2       = 0;
grid_ub_2       = 2 * max(sqrt(cnfctl_emp.cov_ub(i_y,i_y)),sqrt(base.cov(i_y,i_y)));
grid_plot_2     = linspace(grid_lb_2,grid_ub_2,n_kernel)';
dist_plot_2     = kernel(grid_plot_2,winsorize(squeeze(sqrt(cnfctl_emp.cov(i_y,i_y,:))),95),0.4);
dist_plot_2     = dist_plot_2 ./ sum(dist_plot_2 * grid_plot_2(2));

plot_lb_indx_2 = find(dist_plot_2 ~= 0, 1 );
plot_ub_indx_2 = find(dist_plot_2 ~= 0, 1, 'last' );

[~,base_pos_2] = min(abs(grid_plot_2 - sqrt(base.cov(i_y,i_y))));

plot_lb_indx_2 = max(min(plot_lb_indx_2 - n_gap, base_pos_2 - n_gap), 1);
plot_ub_indx_2 = min(max(plot_ub_indx_2 + n_gap, base_pos_2 + n_gap), n_kernel);

subplot(1,3,ii_y)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(ii_y);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
plot(0,0,':','Color',settings.colors.black,'LineWidth',4)
hold on
jbfill(0,0,0,settings.colors.lblue,settings.colors.blue,0,0.5);
hold on
jbfill(0,0,0,settings.colors.lorange,settings.colors.orange,0,0.5);
hold on
plot(grid_plot_1,dist_plot_1,'Color',settings.colors.blue,'LineWidth',4)
hold on
plot(grid_plot_2,dist_plot_2,'Color',settings.colors.orange,'LineWidth',4)
jbfill(grid_plot_1,0*dist_plot_1,dist_plot_1,...
    settings.colors.lblue,settings.colors.lblue,0,0.5);
hold on
jbfill(grid_plot_2,0*dist_plot_2,dist_plot_2,...
    settings.colors.lorange,settings.colors.lorange,0,0.5);
hold on
plot([sqrt(base.cov(i_y,i_y)) sqrt(base.cov(i_y,i_y))],[0 100],':','Color',settings.colors.black,'LineWidth',4)
hold on
if ii_y < 3
    xlim([grid_plot_1(plot_lb_indx_1) grid_plot_1(plot_ub_indx_1)])
else
    xlim([1.95 grid_plot_1(plot_ub_indx_1)])
end
ylim([0 1.2 * max(dist_plot_1)])
% yticks([])
xlabel('St. Dev.','interpreter','latex','FontSize',20)
set(gcf,'color','w')
title(series_names(i_y),'interpreter','latex','fontsize',24)
grid on

if ii_y == 3
    legend({'Data','w/ Extrapolation','Empirics Only'},'Location','Northeast','fontsize',18,'interpreter','latex','NumColumns',1)
end

hold off

end

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.25*pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');

if save_fig == 1 

if indic_early == 0
    print('cnfctl_histograms_optpol_emp','-dpng');
elseif indic_early == 1
    print('cnfctl_histograms_optpol_early_emp','-dpng');
end
    
end

end
