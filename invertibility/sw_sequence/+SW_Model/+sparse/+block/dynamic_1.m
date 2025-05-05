function [y, T] = dynamic_1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(39)=params(26)*y(6)+params(37)*x(2);
  y(40)=params(28)*y(7)+params(38)*x(3);
  y(41)=params(29)*y(8)+params(39)*x(4)+x(2)*params(37)*params(2);
  y(42)=params(31)*y(9)+params(40)*x(5);
  y(43)=params(32)*y(10)+params(41)*x(1);
  y(35)=params(42)*x(6);
  y(34)=params(43)*x(7);
  y(44)=params(33)*y(11)+y(35)-params(8)*y(2);
  y(45)=params(34)*y(12)+y(34)-params(7)*y(1);
end
