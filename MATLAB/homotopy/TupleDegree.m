%
% Extract the tuple degree from a polynomial string
%
% Syntax:  >> d = TupleDegree(f,var)
%
% Input:   f --- (string) polynomial
%        var --- (cell)   variable names in order
%
% Output:  d --- (vector) tuple degree of f w.r.t. var
%
% Example:  >> p = '-2*x^2 + 9*x^3*y^2 - 9*x^4*y^2 + 8*x^5*y^2';
%           >> TupleDegree(p,{'x','y'})
%             ans =
%                     5
%                     3
