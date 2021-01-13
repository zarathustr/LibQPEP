%
% Multiply polynomials, polynomial matrices and numeric matrices
%             h = f1 * f2 * ... * fk
%
% Syntax:   (ptimes is the shortened alias for PolynomialTimes)
%            >> p = PolynomialTimes(f,g)
%            >> p = PolynomialTimes(f, g, h)  %
%   Input:    f, g, h, ...  --- (string/numeric/cell/matrix) 
%                 polynomials, numeric matrices and/or polynomial cells,
%                 so that the multiplication is well-defined.
%
%  Output:    p --- (string) product of the input polynomials
%
%  Example:
%
%     >> PolynomialTimes(3,[1 5; -2 4],'x+y',{'y-1';'x+z'})
%
%     ans =
%
%    '-3*x + 15*x^2 - 3*y + 18*x*y + 3*y^2 + 15*x*z + 15*y*z'
%    '6*x + 12*x^2 + 6*y + 6*x*y - 6*y^2 + 12*x*z + 12*y*z'
%
