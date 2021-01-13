%
%
%  Vector norm for a vector in a product of vector spaces
%
%  Syntax:  >> s = VectorNorm(v, type)
%           >> s = VectorNorm(v)
%
%   Input:  v -- (cell/matrix/polynomial) the vector. If v consists of
%                   components of polyonmials and/or matrices, it should be
%                   a cell.
%        type -- (string, optional) 'Euclidean' or 'infinity'
%                   if not provided, the default is 'Euclidean'
%
%  Output:  s -- the norm
%
%  Example:  >> VectorNorm({'3+2*x-5*y^2',[3 2; 4 1]}, 'infinity')
%            ans  =
%                      5
%
