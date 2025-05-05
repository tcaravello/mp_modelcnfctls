%% PLOT IRFs FOR ESTIMATED MODELS
% Tomas Caravello, Alisdair McKay, Christian Wolf

%% HOUSEKEEPING

task = '/model_estim';

save_fig = 1;

addpath([path vintage '/_auxiliary_functions'])
addpath([path vintage '/model_estim/_results'])
addpath([path vintage '/var_inputs/_results'])
addpath([path vintage '/suff_stats/ratex']);
addpath([path vintage '/suff_stats/behavioral']);

cd([path vintage task]);

%% IMPORTS & SETTINGS

%----------------------------------------------------------------
% News Shock Horizons
%----------------------------------------------------------------

T      = 200;
n_news = 1;

m_target_all = zeros(T,n_news);

m_target_all(16,1) = 1;

% n_news = 2;
% 
% m_target_all = zeros(T,n_news);
% 
% m_target_all(5,1)  = 1;
% m_target_all(16,2) = 1;

%----------------------------------------------------------------
% VAR Targets
%----------------------------------------------------------------

load IRFs_adrr_results

clear Pi_m_draws_emp Pi_m_var_emp R_n_m_draws_emp R_n_m_var_emp Y_m_draws_emp Y_m_var_emp

%----------------------------------------------------------------
% Matched Model IRFs
%----------------------------------------------------------------

% HANK

load params_hank_mode

T = size(Pi_m_hank,1);

m_fit_hank = [m_fit_hank;zeros(T-length(m_fit_hank),2)];

Pi_m_hank_match  = Pi_m_hank * m_fit_hank;
Y_m_hank_match   = Y_m_hank * m_fit_hank;
R_n_m_hank_match = R_n_m_hank * m_fit_hank;

clear param_sol posterior_mode m_fit_hank Pi_m_hank Y_m_hank R_n_m_hank T_use target_Sigma_inv

load hank_draws_main

T       = size(Pi_m_collector,1);
n_draws = size(Pi_m_collector,3);

Pi_m_hank_news_draws  = NaN(T,n_news,n_draws);
Y_m_hank_news_draws   = NaN(T,n_news,n_draws);
R_n_m_hank_news_draws = NaN(T,n_news,n_draws);

for i_news = 1:n_news
    for i_draw = 1:n_draws
        m_target = m_target_all(:,i_news);
        Pi_m_hank_news_draws(:,i_news,i_draw)  = Pi_m_collector(:,:,i_draw) ...
                            * m_target;
        Y_m_hank_news_draws(:,i_news,i_draw)   = Y_m_collector(:,:,i_draw) ...
                            * m_target;
        R_n_m_hank_news_draws(:,i_news,i_draw) = R_n_m_collector(:,:,i_draw) ...
                            * m_target;
    end
end

Pi_m_hank_news_med  = quantile(Pi_m_hank_news_draws,0.5,3);
Y_m_hank_news_med   = quantile(Y_m_hank_news_draws,0.5,3);
R_n_m_hank_news_med = quantile(R_n_m_hank_news_draws,0.5,3);

Pi_m_hank_news_lb   = quantile(Pi_m_hank_news_draws,0.16,3);
Y_m_hank_news_lb    = quantile(Y_m_hank_news_draws,0.16,3);
R_n_m_hank_news_lb  = quantile(R_n_m_hank_news_draws,0.16,3);

Pi_m_hank_news_ub   = quantile(Pi_m_hank_news_draws,0.84,3);
Y_m_hank_news_ub    = quantile(Y_m_hank_news_draws,0.84,3);
R_n_m_hank_news_ub  = quantile(R_n_m_hank_news_draws,0.84,3);

clear m_fit_collector Pi_m_hank_match_draws Y_m_hank_match_draws R_n_m_hank_match_draws ...
    Pi_m_hank_news_draws Y_m_hank_news_draws R_n_m_hank_news_draws ...
    Pi_m_collector Y_m_collector R_n_m_collector i_draw i_news

% RANK

load params_rank_mode

T = size(Pi_m_rank,1);

m_fit_rank = [m_fit_rank;zeros(T-length(m_fit_rank),2)];

Pi_m_rank_match  = Pi_m_rank * m_fit_rank;
Y_m_rank_match   = Y_m_rank * m_fit_rank;
R_n_m_rank_match = R_n_m_rank * m_fit_rank;

clear param_sol posterior_mode m_fit_rank Pi_m_rank Y_m_rank R_n_m_rank T_use target_Sigma_inv

load rank_draws_main

T       = size(Pi_m_collector,1);
n_draws = size(Pi_m_collector,3);

Pi_m_rank_news_draws  = NaN(T,n_news,n_draws);
Y_m_rank_news_draws   = NaN(T,n_news,n_draws);
R_n_m_rank_news_draws = NaN(T,n_news,n_draws);

for i_news = 1:n_news
    for i_draw = 1:n_draws
        m_target = m_target_all(:,i_news);
        Pi_m_rank_news_draws(:,i_news,i_draw)  = Pi_m_collector(:,:,i_draw) ...
                            * m_target;
        Y_m_rank_news_draws(:,i_news,i_draw)   = Y_m_collector(:,:,i_draw) ...
                            * m_target;
        R_n_m_rank_news_draws(:,i_news,i_draw) = R_n_m_collector(:,:,i_draw) ...
                            * m_target;
    end
end

Pi_m_rank_news_med  = quantile(Pi_m_rank_news_draws,0.5,3);
Y_m_rank_news_med   = quantile(Y_m_rank_news_draws,0.5,3);
R_n_m_rank_news_med = quantile(R_n_m_rank_news_draws,0.5,3);

Pi_m_rank_news_lb   = quantile(Pi_m_rank_news_draws,0.16,3);
Y_m_rank_news_lb    = quantile(Y_m_rank_news_draws,0.16,3);
R_n_m_rank_news_lb  = quantile(R_n_m_rank_news_draws,0.16,3);

Pi_m_rank_news_ub   = quantile(Pi_m_rank_news_draws,0.84,3);
Y_m_rank_news_ub    = quantile(Y_m_rank_news_draws,0.84,3);
R_n_m_rank_news_ub  = quantile(R_n_m_rank_news_draws,0.84,3);

clear m_fit_collector Pi_m_rank_news_draws Y_m_rank_news_draws R_n_m_rank_news_draws ...
    Pi_m_collector Y_m_collector R_n_m_collector i_draw i_news

% B-HANK

load params_hank_mode_behav

T = size(Pi_m_hank,1);

m_fit_hank = [m_fit_hank;zeros(T-length(m_fit_hank),2)];

Pi_m_bhank_match  = Pi_m_hank * m_fit_hank;
Y_m_bhank_match   = Y_m_hank * m_fit_hank;
R_n_m_bhank_match = R_n_m_hank * m_fit_hank;

clear param_sol posterior_mode m_fit_hank Pi_m_hank Y_m_hank R_n_m_hank T_use target_Sigma_inv

load hank_draws_main_behav

T       = size(Pi_m_collector,1);
n_draws = size(Pi_m_collector,3);

Pi_m_bhank_news_draws  = NaN(T,n_news,n_draws);
Y_m_bhank_news_draws   = NaN(T,n_news,n_draws);
R_n_m_bhank_news_draws = NaN(T,n_news,n_draws);

for i_news = 1:n_news
    for i_draw = 1:n_draws
        m_target = m_target_all(:,i_news);
        Pi_m_bhank_news_draws(:,i_news,i_draw)  = Pi_m_collector(:,:,i_draw) ...
                            * m_target;
        Y_m_bhank_news_draws(:,i_news,i_draw)   = Y_m_collector(:,:,i_draw) ...
                            * m_target;
        R_n_m_bhank_news_draws(:,i_news,i_draw) = R_n_m_collector(:,:,i_draw) ...
                            * m_target;
    end
end

Pi_m_bhank_news_med  = quantile(Pi_m_bhank_news_draws,0.5,3);
Y_m_bhank_news_med   = quantile(Y_m_bhank_news_draws,0.5,3);
R_n_m_bhank_news_med = quantile(R_n_m_bhank_news_draws,0.5,3);

Pi_m_bhank_news_lb   = quantile(Pi_m_bhank_news_draws,0.16,3);
Y_m_bhank_news_lb    = quantile(Y_m_bhank_news_draws,0.16,3);
R_n_m_bhank_news_lb  = quantile(R_n_m_bhank_news_draws,0.16,3);

Pi_m_bhank_news_ub   = quantile(Pi_m_bhank_news_draws,0.84,3);
Y_m_bhank_news_ub    = quantile(Y_m_bhank_news_draws,0.84,3);
R_n_m_bhank_news_ub  = quantile(R_n_m_bhank_news_draws,0.84,3);

clear m_fit_collector Pi_m_bhank_news_draws Y_m_bhank_news_draws R_n_m_bhank_news_draws ...
    Pi_m_collector Y_m_collector R_n_m_collector i_draw i_news

% B-RANK

load params_rank_mode_behav

T = size(Pi_m_rank,1);

m_fit_rank = [m_fit_rank;zeros(T-length(m_fit_rank),2)];

Pi_m_brank_match  = Pi_m_rank * m_fit_rank;
Y_m_brank_match   = Y_m_rank * m_fit_rank;
R_n_m_brank_match = R_n_m_rank * m_fit_rank;

clear param_sol posterior_mode m_fit_rank Pi_m_rank Y_m_rank R_n_m_rank T_use target_Sigma_inv

load rank_draws_main_behav

T       = size(Pi_m_collector,1);
n_draws = size(Pi_m_collector,3);

Pi_m_brank_news_draws  = NaN(T,n_news,n_draws);
Y_m_brank_news_draws   = NaN(T,n_news,n_draws);
R_n_m_brank_news_draws = NaN(T,n_news,n_draws);

for i_news = 1:n_news
    for i_draw = 1:n_draws
        m_target = m_target_all(:,i_news);
        Pi_m_brank_news_draws(:,i_news,i_draw)  = Pi_m_collector(:,:,i_draw) ...
                            * m_target;
        Y_m_brank_news_draws(:,i_news,i_draw)   = Y_m_collector(:,:,i_draw) ...
                            * m_target;
        R_n_m_brank_news_draws(:,i_news,i_draw) = R_n_m_collector(:,:,i_draw) ...
                            * m_target;
    end
end

Pi_m_brank_news_med  = quantile(Pi_m_brank_news_draws,0.5,3);
Y_m_brank_news_med   = quantile(Y_m_brank_news_draws,0.5,3);
R_n_m_brank_news_med = quantile(R_n_m_brank_news_draws,0.5,3);

Pi_m_brank_news_lb   = quantile(Pi_m_brank_news_draws,0.16,3);
Y_m_brank_news_lb    = quantile(Y_m_brank_news_draws,0.16,3);
R_n_m_brank_news_lb  = quantile(R_n_m_brank_news_draws,0.16,3);

Pi_m_brank_news_ub   = quantile(Pi_m_brank_news_draws,0.84,3);
Y_m_brank_news_ub    = quantile(Y_m_brank_news_draws,0.84,3);
R_n_m_brank_news_ub  = quantile(R_n_m_brank_news_draws,0.84,3);

clear m_fit_collector Pi_m_brank_news_draws Y_m_brank_news_draws R_n_m_brank_news_draws ...
    Pi_m_collector Y_m_collector R_n_m_collector i_draw i_news

%% PLOT RESULTS

%----------------------------------------------------------------
% General Settings
%----------------------------------------------------------------

% plot horizon

IRF_hor_plot = 25;

% colors

settings.colors.black  = [0 0 0];
settings.colors.grey   = [200/255 200/255 200/255];
settings.colors.dgrey  = [130/255 130/255 130/255];
settings.colors.orange = [204/255 102/255 0/255];
settings.colors.green = [37/255 152/255 14/255];
settings.colors.navyblue = [0/255 0/255 50/255];
settings.colors.purple = [160/255 32/255 240/255];
settings.colors.models = [196/255 174/255 120/255; ... % beige
                            204/255 0/255 0/255; ... % red
                            102/255 178/255 255/255]; % blue
settings.colors.lmodels = 0.25 * settings.colors.models + 0.75 * [1 1 1];

% size

plotwidth = 0.27;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2+0.02;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * gapsize + 2 * plotwidth];

% go to results folder

cd([path vintage task '/_results']);

%----------------------------------------------------------------
% Matched IRFs
%----------------------------------------------------------------

scale(1) = 1/max(abs(R_n_m_emp(1:IRF_hor_plot,1,1)));
scale(2) = 1/max(abs(R_n_m_emp(1:IRF_hor_plot,1,2)));

for i_shock = 1:2

figure

subplot(1,3,1)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
if i_shock == 1
    jbfill(0:1:IRF_hor_plot-1,scale(i_shock)*Y_m_lb_emp(1:IRF_hor_plot,1,i_shock)',...
           scale(i_shock)*Y_m_ub_emp(1:IRF_hor_plot,1,i_shock)',settings.colors.grey,settings.colors.grey,0,0.5);
    hold on
    plot(0:1:IRF_hor_plot-1,scale(i_shock)*Y_m_emp(1:IRF_hor_plot,1,i_shock),'Color',settings.colors.dgrey,'LineWidth',4)
    hold on
elseif i_shock == 2
    jbfill(0:1:IRF_hor_plot-1,scale(i_shock)*Y_m_lb_emp(1:IRF_hor_plot,1,i_shock)',...
           scale(i_shock)*Y_m_ub_emp(1:IRF_hor_plot,1,i_shock)',settings.colors.purple,settings.colors.purple,0,0.25);
    hold on
    plot(0:1:IRF_hor_plot-1,scale(i_shock)*Y_m_emp(1:IRF_hor_plot,1,i_shock),'Color',settings.colors.purple,'LineWidth',4)
    hold on
end
plot(0:1:IRF_hor_plot-1,scale(i_shock)*Y_m_rank_match(1:IRF_hor_plot,i_shock),'-','Color',settings.colors.models(1,:),'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot-1,scale(i_shock)*Y_m_hank_match(1:IRF_hor_plot,i_shock),'-','Color',settings.colors.models(2,:),'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot-1,scale(i_shock)*Y_m_brank_match(1:IRF_hor_plot,i_shock),':','Color',settings.colors.models(1,:),'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot-1,scale(i_shock)*Y_m_bhank_match(1:IRF_hor_plot,i_shock),':','Color',settings.colors.models(2,:),'LineWidth',4)
hold on
xlim([0 IRF_hor_plot-1])
ylim([-1 0.5])
yticks([-1:0.25:0.5])
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
hold on
if i_shock == 1
    jbfill(0:1:IRF_hor_plot-1,scale(i_shock)*Pi_m_lb_emp(1:IRF_hor_plot,1,i_shock)',...
           scale(i_shock)*Pi_m_ub_emp(1:IRF_hor_plot,1,i_shock)',settings.colors.grey,settings.colors.grey,0,0.5);
    hold on
    plot(0:1:IRF_hor_plot-1,scale(i_shock)*Pi_m_emp(1:IRF_hor_plot,1,i_shock),'Color',settings.colors.dgrey,'LineWidth',4)
    hold on
elseif i_shock == 2
    jbfill(0:1:IRF_hor_plot-1,scale(i_shock)*Pi_m_lb_emp(1:IRF_hor_plot,1,i_shock)',...
       scale(i_shock)*Pi_m_ub_emp(1:IRF_hor_plot,1,i_shock)',settings.colors.purple,settings.colors.purple,0,0.25);
    hold on
    plot(0:1:IRF_hor_plot-1,scale(i_shock)*Pi_m_emp(1:IRF_hor_plot,1,i_shock),'Color',settings.colors.purple,'LineWidth',4)
    hold on
end
plot(0:1:IRF_hor_plot-1,scale(i_shock)*4*Pi_m_rank_match(1:IRF_hor_plot,i_shock),'-','Color',settings.colors.models(1,:),'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot-1,scale(i_shock)*4*Pi_m_hank_match(1:IRF_hor_plot,i_shock),'-','Color',settings.colors.models(2,:),'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot-1,scale(i_shock)*4*Pi_m_brank_match(1:IRF_hor_plot,i_shock),':','Color',settings.colors.models(1,:),'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot-1,scale(i_shock)*4*Pi_m_bhank_match(1:IRF_hor_plot,i_shock),':','Color',settings.colors.models(2,:),'LineWidth',4)
hold on
xlim([0 IRF_hor_plot-1])
ylim([-0.5 0.5])
yticks([-0.5:0.25:0.5])
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
hold on
if i_shock == 1
    jbfill(0:1:IRF_hor_plot-1,scale(i_shock)*R_n_m_lb_emp(1:IRF_hor_plot,1,i_shock)',...
           scale(i_shock)*R_n_m_ub_emp(1:IRF_hor_plot,1,i_shock)',settings.colors.grey,settings.colors.grey,0,0.5);
    hold on
    plot(0:1:IRF_hor_plot-1,scale(i_shock)*R_n_m_emp(1:IRF_hor_plot,1,i_shock),'Color',settings.colors.dgrey,'LineWidth',4)
    hold on
elseif i_shock == 2
    jbfill(0:1:IRF_hor_plot-1,scale(i_shock)*R_n_m_lb_emp(1:IRF_hor_plot,1,i_shock)',...
       scale(i_shock)*R_n_m_ub_emp(1:IRF_hor_plot,1,i_shock)',settings.colors.purple,settings.colors.purple,0,0.25);
    hold on
    plot(0:1:IRF_hor_plot-1,scale(i_shock)*R_n_m_emp(1:IRF_hor_plot,1,i_shock),'Color',settings.colors.purple,'LineWidth',4)
    hold on
end
plot(0:1:IRF_hor_plot-1,scale(i_shock)*4*R_n_m_rank_match(1:IRF_hor_plot,i_shock),'-','Color',settings.colors.models(1,:),'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot-1,scale(i_shock)*4*R_n_m_hank_match(1:IRF_hor_plot,i_shock),'-','Color',settings.colors.models(2,:),'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot-1,scale(i_shock)*4*R_n_m_brank_match(1:IRF_hor_plot,i_shock),':','Color',settings.colors.models(1,:),'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot-1,scale(i_shock)*4*R_n_m_bhank_match(1:IRF_hor_plot,i_shock),':','Color',settings.colors.models(2,:),'LineWidth',4)
hold on
xlim([0 IRF_hor_plot-1])
ylim([-0.5 1.5])
yticks([-0.5:0.25:1.5])
set(gcf,'color','w')
title('Interest Rate','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
legend({'','','RANK','HANK','B-RANK','B-HANK'},'NumColumns',2,'Location','Northeast','fontsize',18,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.2*2.25*pos(3) 1.1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print(['estim_models_irf_matching_' num2str(i_shock)], '-dpng')

end

%----------------------------------------------------------------
% IRF Extrapolation
%----------------------------------------------------------------

% plot horizon

IRF_hor_plot = 41;

indx_fig = 5;

% RANK vs. HANK

for i_news = 1:n_news

indx_fig = indx_fig + 1;

figure(indx_fig)

subplot(1,3,1)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
jbfill(0:1:IRF_hor_plot-1,Y_m_rank_news_lb(1:IRF_hor_plot,i_news)',...
       Y_m_rank_news_ub(1:IRF_hor_plot,i_news)',settings.colors.lmodels(1,:),settings.colors.lmodels(1,:),0,0.5);
hold on
jbfill(0:1:IRF_hor_plot-1,Y_m_hank_news_lb(1:IRF_hor_plot,i_news)',...
       Y_m_hank_news_ub(1:IRF_hor_plot,i_news)',settings.colors.lmodels(2,:),settings.colors.lmodels(2,:),0,0.5);
hold on
plot(0:1:IRF_hor_plot-1,Y_m_rank_news_med(1:IRF_hor_plot,i_news),'-','Color',settings.colors.models(1,:),'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot-1,Y_m_hank_news_med(1:IRF_hor_plot,i_news),'-','Color',settings.colors.models(2,:),'LineWidth',4)
hold on
xlim([0 IRF_hor_plot-1])
ylim([-0.3 0.2])
yticks([-0.3:0.1:0.2])
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
hold on
jbfill(0:1:IRF_hor_plot-1,4*Pi_m_rank_news_lb(1:IRF_hor_plot,i_news)',...
       4*Pi_m_rank_news_ub(1:IRF_hor_plot,i_news)',settings.colors.lmodels(1,:),settings.colors.lmodels(1,:),0,0.5);
hold on
jbfill(0:1:IRF_hor_plot-1,4*Pi_m_hank_news_lb(1:IRF_hor_plot,i_news)',...
       4*Pi_m_hank_news_ub(1:IRF_hor_plot,i_news)',settings.colors.lmodels(2,:),settings.colors.lmodels(2,:),0,0.5);
hold on
plot(0:1:IRF_hor_plot-1,4*Pi_m_rank_news_med(1:IRF_hor_plot,i_news),'-','Color',settings.colors.models(1,:),'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot-1,4*Pi_m_hank_news_med(1:IRF_hor_plot,i_news),'-','Color',settings.colors.models(2,:),'LineWidth',4)
hold on
xlim([0 IRF_hor_plot-1])
ylim([-0.3 0.1])
yticks([-0.3:0.1:0.1])
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
hold on
plot(0,0,'-','Color',settings.colors.models(1,:),'LineWidth',4)
hold on
plot(0,0,'-','Color',settings.colors.models(2,:),'LineWidth',4)
hold on
jbfill(0:1:IRF_hor_plot-1,4*R_n_m_rank_news_lb(1:IRF_hor_plot,i_news)',...
       4*R_n_m_rank_news_ub(1:IRF_hor_plot,i_news)',settings.colors.lmodels(1,:),settings.colors.lmodels(1,:),0,0.5);
hold on
jbfill(0:1:IRF_hor_plot-1,4*R_n_m_hank_news_lb(1:IRF_hor_plot,i_news)',...
       4*R_n_m_hank_news_ub(1:IRF_hor_plot,i_news)',settings.colors.lmodels(2,:),settings.colors.lmodels(2,:),0,0.5);
hold on
plot(0:1:IRF_hor_plot-1,4*R_n_m_rank_news_med(1:IRF_hor_plot,i_news),'-','Color',settings.colors.models(1,:),'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot-1,4*R_n_m_hank_news_med(1:IRF_hor_plot,i_news),'-','Color',settings.colors.models(2,:),'LineWidth',4)
hold on
xlim([0 IRF_hor_plot-1])
ylim([-0.4 0.3])
yticks([-0.4:0.1:0.3])
set(gcf,'color','w')
title('Interest Rate','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
legend({'RANK','HANK'},'Location','Northeast','fontsize',18,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.2*2.25*pos(3) 1.2*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
if save_fig == 1
    print(['estim_models_rankhank_news_' num2str(i_news)], '-dpng')
end

end

% RANK vs. BRANK

for i_news = 1:n_news

indx_fig = indx_fig + 1;

figure(indx_fig)

subplot(1,3,1)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
jbfill(0:1:IRF_hor_plot-1,Y_m_rank_news_lb(1:IRF_hor_plot,i_news)',...
       Y_m_rank_news_ub(1:IRF_hor_plot,i_news)',settings.colors.lmodels(1,:),settings.colors.lmodels(1,:),0,0.5);
hold on
jbfill(0:1:IRF_hor_plot-1,Y_m_brank_news_lb(1:IRF_hor_plot,i_news)',...
       Y_m_brank_news_ub(1:IRF_hor_plot,i_news)',settings.colors.lmodels(3,:),settings.colors.lmodels(3,:),0,0.5);
hold on
plot(0:1:IRF_hor_plot-1,Y_m_rank_news_med(1:IRF_hor_plot,i_news),'-','Color',settings.colors.models(1,:),'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot-1,Y_m_brank_news_med(1:IRF_hor_plot,i_news),'-','Color',settings.colors.models(3,:),'LineWidth',4)
hold on
xlim([0 IRF_hor_plot-1])
ylim([-0.3 0.2])
yticks([-0.3:0.1:0.2])
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
hold on
jbfill(0:1:IRF_hor_plot-1,4*Pi_m_rank_news_lb(1:IRF_hor_plot,i_news)',...
       4*Pi_m_rank_news_ub(1:IRF_hor_plot,i_news)',settings.colors.lmodels(1,:),settings.colors.lmodels(1,:),0,0.5);
hold on
jbfill(0:1:IRF_hor_plot-1,4*Pi_m_brank_news_lb(1:IRF_hor_plot,i_news)',...
       4*Pi_m_brank_news_ub(1:IRF_hor_plot,i_news)',settings.colors.lmodels(3,:),settings.colors.lmodels(3,:),0,0.5);
hold on
plot(0:1:IRF_hor_plot-1,4*Pi_m_rank_news_med(1:IRF_hor_plot,i_news),'-','Color',settings.colors.models(1,:),'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot-1,4*Pi_m_brank_news_med(1:IRF_hor_plot,i_news),'-','Color',settings.colors.models(3,:),'LineWidth',4)
hold on
xlim([0 IRF_hor_plot-1])
ylim([-0.3 0.1])
yticks([-0.3:0.1:0.1])
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
hold on
plot(0,0,'-','Color',settings.colors.models(1,:),'LineWidth',4)
hold on
plot(0,0,'-','Color',settings.colors.models(3,:),'LineWidth',4)
hold on
jbfill(0:1:IRF_hor_plot-1,4*R_n_m_rank_news_lb(1:IRF_hor_plot,i_news)',...
       4*R_n_m_rank_news_ub(1:IRF_hor_plot,i_news)',settings.colors.lmodels(1,:),settings.colors.lmodels(1,:),0,0.5);
hold on
jbfill(0:1:IRF_hor_plot-1,4*R_n_m_brank_news_lb(1:IRF_hor_plot,i_news)',...
       4*R_n_m_brank_news_ub(1:IRF_hor_plot,i_news)',settings.colors.lmodels(3,:),settings.colors.lmodels(3,:),0,0.5);
hold on
plot(0:1:IRF_hor_plot-1,4*R_n_m_rank_news_med(1:IRF_hor_plot,i_news),'-','Color',settings.colors.models(1,:),'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot-1,4*R_n_m_brank_news_med(1:IRF_hor_plot,i_news),'-','Color',settings.colors.models(3,:),'LineWidth',4)
hold on
xlim([0 IRF_hor_plot-1])
ylim([-0.4 0.3])
yticks([-0.4:0.1:0.3])
set(gcf,'color','w')
title('Interest Rate','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
legend({'RANK','B-RANK'},'Location','Northeast','fontsize',18,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.2*2.25*pos(3) 1.2*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
if save_fig == 1
    print(['estim_models_rankbrank_news_' num2str(i_news)], '-dpng')
end

end

cd([path vintage task]);