%
% Extract the support of input polynomial(s)
%
% Syntax: (psup is the shortened alias for PolynomialSupport)
%          >>   F = PolynomialSupport(f)
%          >>   F = PolynomialSupport(f,g)
%          >>   F = PolynomialSupport(f,g,h)
%
% Input:   f, g, h --- (strings) any number of polynomial strings
%
% Output:   F --- (cell) monomial support of the input polynomials
%
% Example:  >> f = '3 + 2*x*y + 5*x^2*y^3';
%           >> g = '5+3*x*y-6*x^3*y^2';
%           >> S = PolynomialSupport(f,g)
%           
%           S = 
%           
%               '1'    'x*y'    'x^3*y^2'    'x^2*y^3'
%           
