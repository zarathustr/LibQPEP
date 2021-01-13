% Clear tiny/zero coefficients of a multivariate polynomial f
%
%    Syntax:   
%      >> g = mvPolynClear(f)         % clear zero coefficients of f
%      >> g = mvPolynClear(f,epsilon) % clear coefficients whose magnitudes
%                                           are less than or equal to epsilon
%
%     Input:    f -- multivariate polynomial (as a coeff. matrix)
%         epsilon -- magnitude threshold of coefficients to be cleared
%
%    Output:  multivariate polynomial in coefficient matrix
