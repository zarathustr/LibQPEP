% Extract the coefficient of a monomial t in the polynomial p
%
% Syntax:   >> c = Coeff(p,t)
%
%   Input:    p --- (string) polynomial
%             t --- (string) monomial
%
%  Output:    c --- (numeric or string) the coefficient of t in p
%
% Example:   >> p = '-2*x + 9*x*y^2 - 9*x^2*y^2 + 8*x^3*y^2';
%            >> c = Coeff(p,'y^2')
%
%            c =
%
%            9*x - 9*x^2 + 8*x^3
%
%            >> d = Coeff(p,'x*y^2')
%
%            d =
%
%            9
