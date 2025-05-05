function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = SW_Model.static_g1_tt(T, y, x, params);
end
g1 = zeros(33, 33);
g1(1,6)=1;
g1(1,23)=(-params(9));
g1(1,27)=(-(1-params(9)));
g1(2,22)=1;
g1(2,23)=(-T(7));
g1(3,23)=1;
g1(3,24)=1;
g1(3,26)=(-1);
g1(3,27)=(-1);
g1(4,13)=(-1);
g1(4,22)=(-1);
g1(4,24)=1;
g1(5,9)=(-1);
g1(5,25)=(-(T(9)*1/T(11)));
g1(6,7)=(-(1/T(13)));
g1(6,23)=(-((T(3)-(1-params(12)))/T(3)));
g1(6,25)=1-(1-params(12))/T(3);
g1(6,28)=1;
g1(7,7)=(-1);
g1(7,15)=T(17);
g1(7,28)=T(13);
g1(8,3)=1;
g1(8,8)=(-1);
g1(8,15)=(-(1-params(36)-T(6)));
g1(8,16)=(-T(6));
g1(8,22)=(-((T(3)-(1-params(12)))*T(5)));
g1(9,3)=1;
g1(9,6)=(-params(15));
g1(9,24)=(-(params(15)*params(9)));
g1(9,26)=(-(params(15)*(1-params(9))));
g1(10,15)=(-T(18));
g1(10,26)=(-params(20));
g1(10,27)=1;
g1(11,9)=(-(T(4)*T(11)));
g1(11,13)=1-(1-T(4));
g1(11,16)=(-T(4));
g1(12,6)=1;
g1(12,20)=(-(1-params(9)));
g1(12,29)=1;
g1(12,31)=(-params(9));
g1(13,30)=1;
g1(13,31)=(-T(7));
g1(14,20)=(-1);
g1(14,21)=(-1);
g1(14,31)=1;
g1(14,32)=1;
g1(15,14)=(-1);
g1(15,30)=(-1);
g1(15,32)=1;
g1(16,9)=(-1);
g1(16,33)=(-(T(9)*1/T(11)));
g1(17,5)=1;
g1(17,7)=(-(1/T(13)));
g1(17,19)=(-1);
g1(17,31)=(-((T(3)-(1-params(12)))/T(3)));
g1(17,33)=1-(1-params(12))/T(3);
g1(18,5)=T(13);
g1(18,7)=(-1);
g1(18,17)=T(17);
g1(18,19)=(-T(13));
g1(19,4)=1;
g1(19,8)=(-1);
g1(19,17)=(-(1-params(36)-T(6)));
g1(19,18)=(-T(6));
g1(19,30)=(-((T(3)-(1-params(12)))*T(5)));
g1(20,4)=1;
g1(20,6)=(-params(15));
g1(20,21)=(-(params(15)*(1-params(9))));
g1(20,32)=(-(params(15)*params(9)));
g1(21,11)=(-1);
g1(21,19)=1-1/(1+T(8)*params(18))*(T(8)+params(18));
g1(21,29)=(-(1/(1+T(8)*params(18))*T(14)));
g1(22,12)=(-1);
g1(22,17)=(-(T(16)*T(18)));
g1(22,19)=(-(T(15)+params(16)/(1+T(8))-(1+T(8)*params(16))/(1+T(8))));
g1(22,20)=1-(T(9)+T(15)-T(16));
g1(22,21)=(-(params(20)*T(16)));
g1(23,3)=(1-params(25))*params(24);
g1(23,4)=(-((1-params(25))*params(24)));
g1(23,5)=1-params(25);
g1(23,10)=(-1);
g1(23,19)=(-(params(22)*(1-params(25))));
g1(24,6)=1-params(26);
g1(25,7)=1-params(28);
g1(26,8)=1-params(29);
g1(27,9)=1-params(31);
g1(28,10)=1-params(32);
g1(29,2)=(-(1-params(8)));
g1(29,11)=1-params(33);
g1(30,2)=1;
g1(31,1)=(-(1-params(7)));
g1(31,12)=1-params(34);
g1(32,1)=1;
g1(33,9)=(-(params(11)*T(4)*T(10)));
g1(33,14)=1-(1-T(4));
g1(33,18)=(-T(4));

end
