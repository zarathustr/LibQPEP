%
% Convert a polynomial coefficient matrix to a polynomial string
%
% Syntax:  >> p = CoefMat2PolynString(F,var)
%          >> p = CoefMat2PolynString(F,var,digits)
%
% Input:      F --- (matrix) Coefficient matrix of a polynomial
%           var --- (cell)   variable names in order
%        digits --- (integer) optional, number of digits, default = 15
%
% Output:     p --- (string) polynomial as a string
%
% Example:  >> F = [2 0; 1 3; 6, 8]
%              F =
%                2     0
%                1     3
%                6     8
%           >> CoefMat2PolynString(F,{'x','y'})
%              ans =
%                  6*x^2*y + 8*y^3
