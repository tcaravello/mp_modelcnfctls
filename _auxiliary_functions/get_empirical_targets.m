
% Tuning parameters if covariance matrix is non-diagonal

beta_cov_mat = 1;
theta_1_cov_mat = 1;
theta_2_cov_mat = 1;
bandwidth_cov_mat = 8;

global Y_m_emp Y_m_var_emp Pi_m_emp Pi_m_var_emp R_n_m_emp R_n_m_var_emp Y_m_lb_emp Y_m_ub_emp Pi_m_lb_emp Pi_m_ub_emp R_n_m_lb_emp R_n_m_ub_emp

load IRFs_adrr_results
% add the first element, then just put a high variance on it for the RR
% shock
target = [squeeze(Pi_m_emp(1:T_use,1,1))/4; squeeze(Y_m_emp(1:T_use,1,1)); squeeze(R_n_m_emp(1:T_use,1,1))/4;...
    squeeze(Pi_m_emp(1:T_use,1,2))/4; squeeze(Y_m_emp(1:T_use,1,2)); squeeze(R_n_m_emp(1:T_use,1,2))/4;];

% RR is the second shock.

% compute variance-covariance matrix.
%Psi_i_mat = [squeeze(Pi_m_draws_emp(2:T_use,:))/4;Y_m_draws_emp(1:T_use,:);R_n_m_draws_emp(1:T_use,:)/4];
Psi_i_mat = [squeeze(Pi_m_draws_emp(1:T_use,1,1,:))/4; squeeze(Y_m_draws_emp(1:T_use,1,1,:)); squeeze(R_n_m_draws_emp(1:T_use,1,1,:))/4;...
    squeeze(Pi_m_draws_emp(1:T_use,1,2,:))/4; squeeze(Y_m_draws_emp(1:T_use,1,2,:)); squeeze(R_n_m_draws_emp(1:T_use,1,2,:))/4;];

Psi_i_mat_mean = mean(Psi_i_mat,2);
M_draws = size(Psi_i_mat,2);

Psi_i_var = zeros(6*T_use,6*T_use);
for n_draw = 1:M_draws
    Psi_i_var = Psi_i_var + (1/M_draws)*((Psi_i_mat(:,n_draw) -Psi_i_mat_mean)  * (Psi_i_mat(:,n_draw) - Psi_i_mat_mean)') ; 
end

%build triangular kernel

kernel_aux_1 = 1:-1/bandwidth_cov_mat:0;
kernel_aux_2 = [flip(kernel_aux_1(2:end)),kernel_aux_1];
mat_kernel_1 = zeros(T_use,3*T_use);
for t_mat = 1:T_use
mat_kernel_1(t_mat,t_mat:t_mat+2*bandwidth_cov_mat) = kernel_aux_2;
end
mat_kernel_same_var = mat_kernel_1(1:T_use,bandwidth_cov_mat+1:bandwidth_cov_mat+T_use).^theta_1_cov_mat;
mat_kernel_other_var = beta_cov_mat*mat_kernel_1(1:T_use,bandwidth_cov_mat+1:bandwidth_cov_mat+T_use).^theta_2_cov_mat;

mat_aux_1 = eye(6); %on diagonal elements
mat_aux_2 = ones(6,6)-eye(6); %off diagonal elements.

adj_mat_var_diag = kron(eye(6),mat_kernel_same_var);
adj_mat_var = adj_mat_var_diag +kron(ones(6,6)-eye(6),mat_kernel_other_var);


Sigma_var_hat = Psi_i_var.*adj_mat_var;
% put very high variance in impact response for RR shock
Sigma_var_hat(3*T_use+1,3*T_use+1) = 10^(6);
Sigma_var_hat(4*T_use+1,4*T_use+1) = 10^(6);
target_Sigma_inv_non_diag = Sigma_var_hat\eye(size(Sigma_var_hat,1));

Sigma_var_hat_diag = Psi_i_var.*adj_mat_var_diag;
Sigma_var_hat_diag(3*T_use+1,3*T_use+1) = 10^(6);
Sigma_var_hat_diag(4*T_use+1,4*T_use+1) = 10^(6);
target_Sigma_inv_diag = Sigma_var_hat_diag\eye(size(Sigma_var_hat_diag,1));
if cov_mat == 1

target_Sigma_inv = target_Sigma_inv_diag;

else

target_Sigma_inv = target_Sigma_inv_non_diag;
end

