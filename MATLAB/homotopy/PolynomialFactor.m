% Numerical factorization of a polynomial, univariate or multivariate
%    (numerical squarefree factorization or 
%     numerical squarefree irreducible factorization)
%
% Syntax: (pfac is the shortened alias for PolynomilaFactor)
%          >> s = PolynomialFactor(f)
%          >> s = PolynomialFactor(f, tol)
%          >> s = PolynomialFactor(f, tol, style)
%          >> s = PolynomialFactor(f, tol, style, sf)
%          >> s = PolynomialFactor(f, tol, style, sf, gn)
%
% Input:   f --- (string) polynomial
%        tol --- (numeric) coefficient error tolerance
%      style --- (character) 'row' or anything else
%                 if style = 'row', the factorizations will be output
%                    in a row
%                 otherwise, the factorization will be output as a cell
%         sf --- (-1,0,or 1) 
%                   sf=1:           squarefree factorization only
%                   sf=0 (default): squarefree+irreducible factorization
%                   sf=-1:          irreducible factorization only 
%         gn --- (boolean) if gn=0, Gauss-Newton iteration will not be used
%
% Output:  p --- (cell array or string) the default output p is a mx2 cell
%                 where p{k,1} is the k-th factor with multiplicity p{k,2}
%                if the input item style = 'row', the output p is a string
%                in a row showing the factorization
%        res --- (numeric) residual, namely backward error
%       fcnd --- (numeric) condition number
%  Example 1: Univariate factorization of (x-1)^20 (x-2)^15 (x-3)^10 (x-4)^5
%  >> p = ptimes(ppower('x-1',20),ppower('x-2',15),ppower('x-3',10),ppower('x-4',5))
%  >> pfac(p)
%
%  ans = 
%
%    'x-4'    [ 5]
%    'x-3'    [10]
%    'x-2'    [15]
%    'x-1'    [20]
%
%  Example 2: A multivariate factorization
%
%  >> p = '1125 + 1500*x*y*z + 500*x^2*y^2*z^2 + 675*x*y^3*z^2 + 900*x^2*y^4*z^3 + 300*x^3*y^5*z^4 + 135*x^2*y^6*z^4 + 180*x^3*y^7*z^5 + 60*x^4*y^8*z^6 + 9*x^3*y^9*z^6 + 12*x^4*y^10*z^7 + 4*x^5*y^11*z^8'
%
%  >> G = pfac(pp,1e-10,'row')
%
%  G =
%  (1125) * (1 + 0.666666666666667*x*y*z)^2 * (1 + 0.2*x*y^3*z^2)^3
