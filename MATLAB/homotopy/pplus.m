%
%  adding polynomials and polynomial matrices  h = f1 + f2 + ... + fk
%
% Syntax:   (pplus is the shortened alias for PolynomialPlus)
%            >> p = PolynomialPlus(f,g)
%            >> p = PolynomialPlus(f,g,h)
%            >> p = PolynomialPlus(f1,f2,f3,f4)  %
%   Input:    f, g, h, ...  --- (string/numeric/cell/matrix) 
%                 polynomials, numeric matrices and/or polynomial cells,
%                 according to the following addtion rules:
%                 (i) Two matrices having the same size add entrywise
%                (ii) a 1x1 matrix adds to every entry of the other matrix
%
%  Output:    p --- (string) sum of the input polynomials
%
%  Example:
%  >> PolynomialPlus('3+x^2','(5+3i)*x*z',zeros(2,2))
%  ans = 
%
%     '3 + x^2 + (5+3i)*x*z'    '3 + x^2 + (5+3i)*x*z'
%     '3 + x^2 + (5+3i)*x*z'    '3 + x^2 + (5+3i)*x*z'
%
