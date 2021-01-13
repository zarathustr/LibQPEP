%
%  Compute the factorization of a multivariate polynomial 
%  by Ruppert matrix, projection, Newton polygon, generalized eigenvalue
%  and numerical GCD
%
%  Syntax:   >> fac = NIF(F,tol)
% 
%    Input:    F -- (matrix) multivariate polynomial in coeff. matrices
%                            !!!must be squarefree!!!
%            tol -- (numeric) coefficient error tolerence
%
%    Output: fac -- (cell) irreducible factors of G in coeff. matrices
%  
%  Example:
%
%        >> F
%
%        F =
%
%             0     1     1     2
%             0     1     0     1
%             0     0     2     2
%            15    10     3     2
%
%        >> fac = NIF(F,1e-10);
%        >> fac{:}
%
%        ans =
%        
%                0             1.0000          
%                0             1.0000          
%                0                  0          
%           0.8320 - 0.0035i   0.5547 - 0.0024i
%        
%        
%        ans =
%        
%                0             1.0000          
%                0                  0          
%                0             2.0000          
%           0.9729 + 0.1223i   0.1946 + 0.0245i
%         
