%
% Clear tiny coefficients, tiny real parts, and tiny imaginary parts
% from a polynomial or a cell array of polynomials
%
% Syntax:  (pclear is the shortened alias of PolynomialClear)
%          >> g = PolynomialClear(f)
%          >> g = PolynomialClear(f,tol)
%
%   Input:   f --- (string/numeric/cell array) polynomial (array)
%          tol --- (numeric)           threshold for being tiny
%
%  Output:   g --- (string/numeric/cell array) cleared polynomial (array)
%
%  Example:
%
%   >> PolynomialClear('(3+2e-13*i) + 2.1e-15*x*y+(1e-14+2*i)*x^3',1e-10)
%
%   ans =
%
%   3 + (0+2i)*x^3
