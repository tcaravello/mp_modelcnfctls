%----------------------------------------------------------------
% Empirical Shock
%----------------------------------------------------------------

load IRFs_adrr_results

T = size(Pi_m_emp,1);
N_draws = size(Pi_m_draws_emp,3);

clear Pi_m_draws Y_m_draws I_m_draws
% No need to multiply or divide by 4 since empirical IRFs are all measured
% in annualized terms.
Pi_m_draws(:,:,:) = squeeze(Pi_m_draws_emp);
Y_m_draws(:,:,:) = squeeze(Y_m_draws_emp);
I_m_draws(:,:,:) = squeeze(R_n_m_draws_emp);

Pi_m_base  = Pi_m_emp;
I_m_base   = R_n_m_emp;
Y_m_base   = Y_m_emp;

% % temporary: only keep first shock
% 
% Pi_m_draws = Pi_m_draws(:,1,:);
% Y_m_draws  = Y_m_draws(:,1,:);
% I_m_draws  = I_m_draws(:,1,:);
% 
% Pi_m_base  = Pi_m_base(:,:,1);
% I_m_base   = I_m_base(:,:,1);
% Y_m_base   = Y_m_base(:,:,1);