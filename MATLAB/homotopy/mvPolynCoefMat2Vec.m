%
%  Convert a coefficient matrix to a (sparse) coefficient vector
%
%    Syntax:  >> v = mvPolynCoefMat2Vec(F)
%             >> v = mvPolynCoefMat2Vec(F,d)
%
%    Input:      F -- multivariate polynomials in coeff. matrix
%                d -- (optional) degree bound for the polynomial vector space
%                     if not provided, the tuple degree of F is used
%
%   Output:      v -- (sparse vector) the coefficient vector
%
%   Example:
%      >> F = PolynString2CoefMat('3 + 5.5*x^2 - 7*x*y^3',{'x','y'})
%
%      F =
%
%               0    2.0000    1.0000
%               0         0    3.0000
%          3.0000    5.5000   -7.0000
%
%      >> f = mvPolynCoefMat2Vec(F)
%
%      ans =
%
%         (1,1)       3.0000
%         (3,1)       5.5000
%        (11,1)      -7.0000
%
