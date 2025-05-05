function [T_order, T] = static_resid_tt(y, x, params, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 16
    T = [T; NaN(16 - size(T, 1), 1)];
end
T(1) = 1+params(35)/100;
T(2) = 1/(1+params(6)/100);
T(3) = T(2)^(-1)*T(1)^params(13);
T(4) = 1-(1-params(12))/T(1);
T(5) = params(15)*((1-params(9))/params(9)*(T(3)-(1-params(12)))/(params(9)^params(9)*(1-params(9))^(1-params(9))/(params(15)*(T(3)-(1-params(12)))^params(9)))^(1/(1-params(9))))^(params(9)-1);
T(6) = T(1)*T(4)*T(5);
T(7) = 1/(params(10)/(1-params(10)));
T(8) = T(1)*T(2)*T(1)^(-params(13));
T(9) = 1/(1+T(8));
T(10) = T(1)^2;
T(11) = params(11)*T(10);
T(12) = params(14)/T(1);
T(13) = (1-T(12))/(params(13)*(1+T(12)));
T(14) = (1-params(19))*(1-T(8)*params(19))/params(19)/(1+(params(15)-1)*params(3));
T(15) = T(8)/(1+T(8));
T(16) = 1/(1+(params(21)-1)*params(1))*(1-params(17))*(1-T(8)*params(17))/((1+T(8))*params(17));
end
