% Project a polynomial to a vector space of fewnomials
%
% Syntax: (pproj is the shortened alias for PolynomilProject2Space)
%          >> F = PolynomialProject2Space(f, S)
%          >> F = PolynomialProject2Space(f, S, type)
%
% Input:   f --- (string) polynomial
%          S --- (cell/string) The fewnomial space basis info. 
%                   (i) if S is a character string (representing a polynomial
%                           the support of S is the basis
%                  (ii) if S is a cell, it must contain
%                          S{1} -- the variables in a cell
%                          S{2} -- the tuple degree bound
%                          S{3} -- the indices of the monomial basis in
%                               lexicographical order
%       type --- (string) 'polyn' or 'mat' indicating the output type to be
%                   a polynomial string or coefficient matrix (default)
%
% Output:  F --- (string/matrix) the polynomial projected from f to S
%
% Example: >> f = '3+5*x^2-6*x^3+7*x^7+9*x^9-x^5*y';
%          >> S = '1+x+x^3+x^4+x^9';
%          >> F = PolynomialProject2Space(f,S,'polyn');
%           F =
%           3 - 0*x - 6*x^3 - 0*x^4 + 9*x^9
