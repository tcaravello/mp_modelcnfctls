function [Pi_m_collector, Y_m_collector, R_n_m_collector, m_fit_collector]...
    = sample_from_models(M_draws, model_probs, Pi_m_col, Y_m_col, R_n_m_col, m_fit_col)

global shocks_match
% create uniform draws
unif_all_models = randsample(length(model_probs), M_draws, true,model_probs);
model_uni_draws = randsample(M_draws, M_draws, true);

% placeholders 

Pi_m_collector = zeros(size(Pi_m_col, [1 2 3]));
Y_m_collector = zeros(size(Y_m_col, [1 2 3]));
R_n_m_collector = zeros(size(R_n_m_col, [1 2 3]));
if shocks_match == 2
    m_fit_collector = zeros(size(m_fit_col, [1 2 3]));
else
    m_fit_collector = zeros(size(m_fit_col, [1 2]));
end


% draw from each model
for i = 1: M_draws
Pi_m_collector(:,:,i)  = Pi_m_col(:,:,model_uni_draws(i),unif_all_models(i));
Y_m_collector(:,:,i)  = Y_m_col(:,:,model_uni_draws(i),unif_all_models(i));
R_n_m_collector(:,:,i) = R_n_m_col(:,:,model_uni_draws(i),unif_all_models(i));
if shocks_match == 2
    m_fit_collector(:,:,i) = m_fit_col(:,:,model_uni_draws(i), unif_all_models(i));
else
    m_fit_collector(:,i) = m_fit_col(:,model_uni_draws(i), unif_all_models(i));
end
end

end