%
% Extract the total degree from a polynomial string
%
% Syntax:  >> d = TotalDegree(f)
%          >> d = TotalDegree(f,var)
%
% Input:   f --- (string) polynomial
%        var --- (cell)   optional, variable names 
%
% Output:  d --- (integer) total degree of f 
%
% Example:  >> p = '-2*x^2 + 9*x^3*y^2 - 9*x^4*y^2 + 8*x^5*y^2';
%           >> TotalDegree(p)
%             ans =
%                     7
%
