%
%  Project an m-variate polynomial to an n-variate polynomial
%    using given values on the m-n variables
%
%    syntax:  G = mvPolynProjection(F, x, ix)
%
%     input:  F --- a multivariate polynomial in coef. matrix
%             x --- a vector of m-n values to be used for projection
%            ix --- a vector of indices of the m-n variables
%
%    output:  G --- the projected polynomial
%
%   Example:  To project a polynomial F(x1,x2,x3,x4) to F(x1,1.6,x2,6.8):
%        >> G = mvPolynProject(F,[1.6,6,8],[2,4])
%
