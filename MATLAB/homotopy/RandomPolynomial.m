%
% Generate a random polynomial
%
% Syntax: >> p = RandomPolynomial(variables,degree,terms)
%     or  >> p = prand(variables,degree,terms)
%
% Input: variables --- (cell or string) variables
%           degree --- (integer or vector) tuple degree bound
%                               for the random polynomial
%            terms --- (integer/string, optional) number of terms or 'all',
%                       if missing, produce random number of terms.
%
%  Output:       p --- (string) random polynomial
%
% Example:  >>  p = RandomPolynomial({'x','y'},[5,2],4)
%           p = 
%                -2*x^2 + 9*x^3*y^2 - 9*x^4*y^2 + 8*x^5*y^2
