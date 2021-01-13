%
% convert a polynomial string to a coefficient matrix
%
% Syntax: >> F = PolynString2CoefMat(f,var)
%
% Input:    f --- (string) polynomial string
%         var --- (cell) variable names in order
%
% Output:   F --- (matrix) polynomial in coefficient matrix
%
% Example:  >> p = '6*x^2*y + 8*y^3';
%           >> PolynString2CoefMat(p,{'x','y'})
%           ans =
%                  2     0
%                  1     3
%                  6     8
