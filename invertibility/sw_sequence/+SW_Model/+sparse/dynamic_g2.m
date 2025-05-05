function [g2_v, T_order, T] = dynamic_g2(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(18, 1);
end
[T_order, T] = SW_Model.sparse.dynamic_g2_tt(y, x, params, steady_state, T_order, T);
g2_v = NaN(0, 1);
end
