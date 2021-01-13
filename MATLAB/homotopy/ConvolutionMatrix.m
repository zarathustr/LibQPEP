%  Generate an n-column convolution matirx of
%  An n-column convolution matrix of p is the matrix for the 
%  linear transformation L(f) = p*f  for polynomial f of degree n-1.
%  If f is the coefficient vector of f, and C is the convolution matrix,
%  then C*f is the coefficient vector of p*f.
%
%    Syntax:  C = ConvolutionMatrix(p,n)
%
%    Input:   p --- (string or vector) the univariate polynomial
%             n --- (integer) the column dimension of the convolution matrix
%
%   Output:   C --- the convolution matrix
%
%   Example:  >> p = '1+3*x^4-6*x^8';
%             >> C = ConvolutionMatrix(p,3)
%             C =
%             -6     0     0
%              0    -6     0
%              0     0    -6
%              0     0     0
%              3     0     0
%              0     3     0
%              0     0     3
%              0     0     0
%              1     0     0
%              0     1     0
%              0     0     1
