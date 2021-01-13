%
%  Univariate factorization that is capable of identifying multiple roots
%  and multiplicities accurately even if the polynomial is perturbed
%
%   Input    p -- coefficient vector of the polynomial to be factored
%          tol -- (optional) backward error tolerance
%    showroots -- (optional, 0 or 1) showing roots or not
%
%  Output    F -- kx3 matrix containing factors with
%                    each row [a,b,m] representing a factor (a*x+b)^m
%
%   Syntax:  >> LHS = RHS 
%            where LHS and RHS can be any one of the following:
%
%               LHS options:   |  RHS options:
%               ---------------|------------------
%                 F            |  uvFactor(p)
%                [F,res]       |  uvFactor(p,tol)
%                [F,res,fcnd]  |  uvFactor(p,tol,1)
%
