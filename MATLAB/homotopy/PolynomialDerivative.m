% Calculate a (partial) derivative of a polynomial f
%
% Syntax:   (pder is the shortened alias for PolynomialDerivative)
%   >> p = PolynomialDerivative(f)        % if f is univariate
%   >> p = PolynomialDerivative(f,'x',k)  % k-th derivative of f(x)
%   >> p = PolynomialDerivative(f,z)      % partial derivative of f w.r.t. all variables in z
%   >> p = PolynomialDerivative(f,z,k)
%
%  Input:   f --- (string)            polynomial
%           z --- (string or cell)    variable name(s)
%           k --- (integer or vector, optional) derivative orders, 
%                                        assumed to be ones if not provided
%
% Output:   p --- (string)  (partial) derivative of f
%
% Example:  >> f = '3*x^5+2*x+8';
%           >> PolynomialDerivative(f)
%
%           ans =
%
%           2 + 15*x^4
%
%           >> g = '-2*x + 9*x*y^2 - 9*x^2*y^2 + 8*x^3*y^2';
%           >> PolynomialDerivative(g,{'x','y'},[2,1])
%
%           ans =
%
%           -36*y + 96*x*y
