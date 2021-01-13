%
% Evaluate a polynomial at given values of (some or all) values. The result
% can be a number or a polynomial of remaining variables
%
% Syntax: (peval is the shortened alias for PolynomilaEvaluate)
%          >> s = PolynomialEvaluate(f, var, z)
%
% Input:   f --- (string/cell) polynomial or a cell of polynomials
%        var --- (cell)   variable names
%          z --- (vector) values of the variables in var
%
% Output:  s --- (numeric or string or cell) value of the polynomials
%
% Example 1: >> g = '-2*x + 9*x*y^2 - 9*z^2*y^2 + 8*x^3*y^2*z^4';
%          >> PolynomialEvaluate(g,{'x','z'},[2,1])
%          ans =
%          -4 + 73*y^2
% Example 2: >> F = {'2*x + 7*x*y^2 + y^2*z + 6*x*z^2 - 5*x^3*z^2'
%                     '-6*x^2 - 8*x^3*z - 5*z*y - 4*x*z*y^2'}
%            >> peval({p,q},{'x','y'},[2,3])
%            ans = 
%                '4 + 135*z + 12*z^2'    '-24 - 64*z'
