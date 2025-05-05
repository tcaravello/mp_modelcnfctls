function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = SW_Model.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(33, 1);
    residual(1) = (y(26)) - (params(9)*y(43)+(1-params(9))*y(47));
    residual(2) = (y(42)) - (y(43)*T(3));
    residual(3) = (y(43)) - (y(47)+y(46)-y(44));
    residual(4) = (y(44)) - (y(42)+y(13));
    residual(5) = (y(36)) - (y(29)+T(14)*(y(45)*1/T(5)+y(16)+y(55)*T(13)));
    residual(6) = (y(45)) - (y(27)*1/T(7)-y(48)+y(61)*(T(10)-(1-params(12)))/T(10)+y(62)*(1-params(12))/T(10));
    residual(7) = (y(35)) - (y(27)+y(15)*T(6)/(1+T(6))+y(54)*1/(1+T(6))+(y(46)-y(63))*T(15)-y(48)*T(7));
    residual(8) = (y(23)) - (y(28)+y(35)*(1-params(36)-T(1)*T(11)*T(12))+y(36)*T(1)*T(11)*T(12)+y(42)*(T(10)-(1-params(12)))*T(12));
    residual(9) = (y(23)) - (params(15)*(y(26)+params(9)*y(44)+(1-params(9))*y(46)));
    residual(10) = (y(47)) - (y(46)*params(20)+y(35)*T(8)-y(15)*T(9));
    residual(11) = (y(33)) - (y(13)*(1-T(11))+y(36)*T(11)+y(29)*T(5)*T(11));
    residual(12) = (y(49)) - (params(9)*y(51)+(1-params(9))*y(40)-y(26));
    residual(13) = (y(50)) - (T(3)*y(51));
    residual(14) = (y(51)) - (y(40)+y(41)-y(52));
    residual(15) = (y(52)) - (y(50)+y(14));
    residual(16) = (y(38)) - (y(29)+T(14)*(y(53)*1/T(5)+y(18)+y(57)*T(13)));
    residual(17) = (y(53)) - (y(58)-y(25)+y(27)*1/T(7)+y(64)*(T(10)-(1-params(12)))/T(10)+y(65)*(1-params(12))/T(10));
    residual(18) = (y(37)) - (y(27)+y(17)*T(6)/(1+T(6))+y(56)*1/(1+T(6))+(y(41)-y(60))*T(15)-(y(25)-y(58))*T(7));
    residual(19) = (y(24)) - (y(28)+y(37)*(1-params(36)-T(1)*T(11)*T(12))+y(38)*T(1)*T(11)*T(12)+y(50)*(T(10)-(1-params(12)))*T(12));
    residual(20) = (y(24)) - (params(15)*(y(26)+params(9)*y(52)+(1-params(9))*y(41)));
    residual(21) = (y(39)) - (y(31)+T(16)*(params(18)*y(19)+y(58)*T(13)+y(49)*T(17)));
    residual(22) = (y(40)) - (y(32)+y(20)*T(14)+y(59)*T(13)/(1+T(13))+y(19)*params(16)/(1+T(13))-y(39)*(1+params(16)*T(13))/(1+T(13))+y(58)*T(13)/(1+T(13))+(params(20)*y(41)+y(37)*T(8)-y(17)*T(9)-y(40))*T(18));
    residual(23) = (y(25)) - (y(39)*params(22)*(1-params(25))+(1-params(25))*params(24)*(y(24)-y(23))+(1-params(25))*params(23)*(y(24)-y(23)-y(4)+y(3))+params(25)*y(5)+y(30));
    residual(24) = (y(26)) - (params(26)*y(6)+params(37)*x(it_, 2));
    residual(25) = (y(27)) - (params(28)*y(7)+params(38)*x(it_, 3));
    residual(26) = (y(28)) - (params(29)*y(8)+params(39)*x(it_, 4)+x(it_, 2)*params(37)*params(2));
    residual(27) = (y(29)) - (params(31)*y(9)+params(40)*x(it_, 5));
    residual(28) = (y(30)) - (params(32)*y(10)+params(41)*x(it_, 1));
    residual(29) = (y(31)) - (params(33)*y(11)+y(22)-params(8)*y(2));
    residual(30) = (y(22)) - (params(42)*x(it_, 6));
    residual(31) = (y(32)) - (params(34)*y(12)+y(21)-params(7)*y(1));
    residual(32) = (y(21)) - (params(43)*x(it_, 7));
    residual(33) = (y(34)) - (y(14)*(1-T(11))+y(38)*T(11)+y(29)*params(11)*T(4)*T(11));

end
