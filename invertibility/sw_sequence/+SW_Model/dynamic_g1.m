function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = SW_Model.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(33, 72);
g1(1,26)=1;
g1(1,43)=(-params(9));
g1(1,47)=(-(1-params(9)));
g1(2,42)=1;
g1(2,43)=(-T(3));
g1(3,43)=1;
g1(3,44)=1;
g1(3,46)=(-1);
g1(3,47)=(-1);
g1(4,13)=(-1);
g1(4,42)=(-1);
g1(4,44)=1;
g1(5,29)=(-1);
g1(5,16)=(-T(14));
g1(5,36)=1;
g1(5,55)=(-(T(13)*T(14)));
g1(5,45)=(-(1/T(5)*T(14)));
g1(6,27)=(-(1/T(7)));
g1(6,61)=(-((T(10)-(1-params(12)))/T(10)));
g1(6,45)=1;
g1(6,62)=(-((1-params(12))/T(10)));
g1(6,48)=1;
g1(7,27)=(-1);
g1(7,15)=(-(T(6)/(1+T(6))));
g1(7,35)=1;
g1(7,54)=(-(1/(1+T(6))));
g1(7,46)=(-T(15));
g1(7,63)=T(15);
g1(7,48)=T(7);
g1(8,23)=1;
g1(8,28)=(-1);
g1(8,35)=(-(1-params(36)-T(1)*T(11)*T(12)));
g1(8,36)=(-(T(1)*T(11)*T(12)));
g1(8,42)=(-((T(10)-(1-params(12)))*T(12)));
g1(9,23)=1;
g1(9,26)=(-params(15));
g1(9,44)=(-(params(15)*params(9)));
g1(9,46)=(-(params(15)*(1-params(9))));
g1(10,15)=T(9);
g1(10,35)=(-T(8));
g1(10,46)=(-params(20));
g1(10,47)=1;
g1(11,29)=(-(T(5)*T(11)));
g1(11,13)=(-(1-T(11)));
g1(11,33)=1;
g1(11,36)=(-T(11));
g1(12,26)=1;
g1(12,40)=(-(1-params(9)));
g1(12,49)=1;
g1(12,51)=(-params(9));
g1(13,50)=1;
g1(13,51)=(-T(3));
g1(14,40)=(-1);
g1(14,41)=(-1);
g1(14,51)=1;
g1(14,52)=1;
g1(15,14)=(-1);
g1(15,50)=(-1);
g1(15,52)=1;
g1(16,29)=(-1);
g1(16,18)=(-T(14));
g1(16,38)=1;
g1(16,57)=(-(T(13)*T(14)));
g1(16,53)=(-(1/T(5)*T(14)));
g1(17,25)=1;
g1(17,27)=(-(1/T(7)));
g1(17,58)=(-1);
g1(17,64)=(-((T(10)-(1-params(12)))/T(10)));
g1(17,53)=1;
g1(17,65)=(-((1-params(12))/T(10)));
g1(18,25)=T(7);
g1(18,27)=(-1);
g1(18,17)=(-(T(6)/(1+T(6))));
g1(18,37)=1;
g1(18,56)=(-(1/(1+T(6))));
g1(18,58)=(-T(7));
g1(18,41)=(-T(15));
g1(18,60)=T(15);
g1(19,24)=1;
g1(19,28)=(-1);
g1(19,37)=(-(1-params(36)-T(1)*T(11)*T(12)));
g1(19,38)=(-(T(1)*T(11)*T(12)));
g1(19,50)=(-((T(10)-(1-params(12)))*T(12)));
g1(20,24)=1;
g1(20,26)=(-params(15));
g1(20,41)=(-(params(15)*(1-params(9))));
g1(20,52)=(-(params(15)*params(9)));
g1(21,31)=(-1);
g1(21,19)=(-(params(18)*T(16)));
g1(21,39)=1;
g1(21,58)=(-(T(13)*T(16)));
g1(21,49)=(-(T(16)*T(17)));
g1(22,32)=(-1);
g1(22,17)=(-(T(18)*(-T(9))));
g1(22,37)=(-(T(8)*T(18)));
g1(22,19)=(-(params(16)/(1+T(13))));
g1(22,39)=(1+params(16)*T(13))/(1+T(13));
g1(22,58)=(-(T(13)/(1+T(13))));
g1(22,20)=(-T(14));
g1(22,40)=1+T(18);
g1(22,59)=(-(T(13)/(1+T(13))));
g1(22,41)=(-(params(20)*T(18)));
g1(23,3)=(-((1-params(25))*params(23)));
g1(23,23)=(-((-((1-params(25))*params(24)))-(1-params(25))*params(23)));
g1(23,4)=(1-params(25))*params(23);
g1(23,24)=(-((1-params(25))*params(24)+(1-params(25))*params(23)));
g1(23,5)=(-params(25));
g1(23,25)=1;
g1(23,30)=(-1);
g1(23,39)=(-(params(22)*(1-params(25))));
g1(24,6)=(-params(26));
g1(24,26)=1;
g1(24,67)=(-params(37));
g1(25,7)=(-params(28));
g1(25,27)=1;
g1(25,68)=(-params(38));
g1(26,8)=(-params(29));
g1(26,28)=1;
g1(26,67)=(-(params(37)*params(2)));
g1(26,69)=(-params(39));
g1(27,9)=(-params(31));
g1(27,29)=1;
g1(27,70)=(-params(40));
g1(28,10)=(-params(32));
g1(28,30)=1;
g1(28,66)=(-params(41));
g1(29,2)=params(8);
g1(29,22)=(-1);
g1(29,11)=(-params(33));
g1(29,31)=1;
g1(30,22)=1;
g1(30,71)=(-params(42));
g1(31,1)=params(7);
g1(31,21)=(-1);
g1(31,12)=(-params(34));
g1(31,32)=1;
g1(32,21)=1;
g1(32,72)=(-params(43));
g1(33,29)=(-(params(11)*T(4)*T(11)));
g1(33,14)=(-(1-T(11)));
g1(33,34)=1;
g1(33,38)=(-T(11));

end
