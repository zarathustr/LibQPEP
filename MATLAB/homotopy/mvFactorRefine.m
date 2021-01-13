%  For a given initial factorization (stored in F) of multivariate polynomial p,
%  the module  mvFactorRefine  attempts to improve the accuracy by applying
%  the Gauss-Newton iteration module GaussNewton. 
%
%         
%   Syntax:  >>[G,res,cond] = mvFactorRefine(F, p);
%
%  Input:  F -- a kx2 cell containing the factors and multiplicities of the 
%                  initial factorization of p.  For each j = 1, 2, ..., k,
%                  F{j,1} is the j-th factor with multiplicity F{j,2}.
%          p -- multivariate polynomial in coeff. matrix
%
%  Output:
%          G -- a cell containing the refined factors in the same format as F
%        res -- the residual, i.e., backward error of the factorization
%       cond -- the condition number of the factorization
