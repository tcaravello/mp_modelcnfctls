function [T_order, T] = static_g1_tt(y, x, params, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = SW_Model.sparse.static_resid_tt(y, x, params, T_order, T);
T_order = 1;
if size(T, 1) < 18
    T = [T; NaN(18 - size(T, 1), 1)];
end
T(17) = 1-(T(12)/(1+T(12))+1/(1+T(12)));
T(18) = 1/(1-T(12))-T(12)/(1-T(12));
end
