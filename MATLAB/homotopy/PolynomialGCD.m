% Compute the numerical greatest common divisor of two polynomials.
%   The polynomials can be either univariate or multivariate
%
% Syntax: (pgcd is the shortened alias of PolynomialGCD)
%      >> [u, v, w, residual, condition] = PolynomialGCD(f, g)
%      >> [u, v, w, residual, condition] = PolynomialGCD(f, g, tol)
%      >> [u, v, w, residual, condition, p, q] = PolynomialGCD(f, g)
%      >> [u, v, w, residual, condition, p, q] = PolynomialGCD(f, g, tol)
%
%     Input:   f --- (string)            polynomial one
%              g --- (string)            polynomial two
%            tol --- (vector, optional)  error tolerace in coefficients
%                                          with default 1.0e-10
%
%    Output:   u --- (string)  GCD of f and g
%              v --- (string)  cofactor of f such that f = u*v
%              w --- (string)  cofactor of g such that g = u*v
%       residual --- (numeric) backward error ||(f,g)-(u*v,u*w)||
%      condition --- (numeric) sensitivity measure of the numerical GCD
%            p,q --- (strings) polynomials for p*f+q*g=u in univariate case 
%
% Example:
% >> f = '-45*x*y - 15*x^3*y - 20*x*y^2 + 27*x*y^3 + 9*x^3*y^3 + 12*x*y^4';
% >> g = '45*x^2*y^2 + 15*x^2*y^3 - 27*x^2*y^4 - 9*x^2*y^5';
% >>  [u,v,w,res,cond] = PolynomialGCD(f,g);
% >>  {u; v; w; res; cond}
%
% ans = 
%       '-70*x*y + 42*x*y^3'    
%       '0.642857142857143 + 0.214285714285714*x^2 + 0.285714285714286*y'
%        '-0.642857142857143*x*y - 0.214285714285714*x*y^2'
%        [3.552713678800501e-014]
%        [     1.089132505368521]
%
