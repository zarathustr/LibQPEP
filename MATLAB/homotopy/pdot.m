% Calculate the dot product of coefficient vectors of two polynomials
% 
%  Syntax:  (pdot is the shortened alias for PolynomialDotProduct)
%    >> d = PolynomialDotProduct(f,g)        
%    >> d = pdot(f,g)
% 
%   Input:   f --- (string)            polynomial
%            g --- (string)            polynomial
% 
%  Output:   d --- (numeric)  dot product of coeffcient vectors of f and g
% 
%  Example:  >> d = pdot('3+x*y-y^3*x^2+x^5','5-4*x^2*y^3+9*x*y')
%            d =
%                  28
