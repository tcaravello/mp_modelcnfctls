function [y, T] = static_6(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(2)=params(42)*x(6);
  y(1)=params(43)*x(7);
end
