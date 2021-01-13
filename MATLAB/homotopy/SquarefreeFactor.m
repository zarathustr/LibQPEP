%  Multivariate Squarefree factorization that is capable of identifying 
%  squarefree factors and multiplicities accurately even if the polynomial 
%  is perturbed
%
%   Syntax:  >>[F,res,cond] = SquarefreeFactor(p,tol);
%
%   Input   p -- a multivariate polynomial (as a coeff. matrix) to be factored
%         tol -- (optional) backward error tolerance
%
%  Output   F -- cell of size kx2. For each j = 1, 2, ..., k, F{j,1}, F{j,2}
%                    are a squarefree factor and its multiplicity respectively
%         res -- residual
%        cond -- condition number of the factorization
