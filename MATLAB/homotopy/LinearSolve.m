%
%  Solve for the general least squares solution z of L(z) = b 
%     Solve for the general numerical solution of  
%        (i)  matrix-vector equation   A*z = b  for z = z0+Kernel(A) 
%        (ii) linear operator equation    L(z) = b  for z = z0+Kernel(L) 
%     within an error tolerance tol, where z0 is the minimum norm solution 
%     and the kernel is the numerical kernle within the error tolerance tol 
%    
%       Syntax for solving A*z = b (within error tolerance tol)   
%       >> [z, K, cnd, res] = LinearSolve(A, b) 
%       >> [z, K, cnd, res] = LinearSolve(A, b, tol)   
%    
%       syntax for solving L(z) = b (within error tolerance tol) 
%       >> [z, K, cnd, res] = LinearSolve({L, Domain, para}, b) 
%       >> [z, K, cnd, res] = LinearSolve({L, Domain, para}, b, tol) 
%       >> [z, K, cnd, res] = LinearSolve({L, Domain, para, r}, b) 
%       >> [z, K, cnd, res] = LinearSolve({L, Domain, para, r}, b, tol) 
%    
%       Input (case i: Matrix-vector equation A*z = b):  
%               A, b  -- coefficient matrix and right-side vector in A*z=b 
%                tol  -- error tolerance 
%    
%       Input (case ii: Linear operator equation L(z) = b):   
%                L  -- Matlab function handle for calculating the linear 
%                       transformation with syntax >> L(X1,...,Xm,Q1,...,Qk)  
%                       where X1,...,Xm and Q1,...,Qk are variables and 
%                       parameters of L respectively. Each one of X1,...,Xm is  
%                          (i) a matrix (including vector), or 
%                         (ii) a polynomial 
%                              to be transformed by L. 
%                       The remaining parameters Q1,...,Qk must be 
%                           provided as cell entries for para. 
%                       If the codomain is a product space, the output of the 
%                           function L must be a cell array 
%           Domain  -- domain of the linear transformation L represented by 
%                         (i) an mxn matrix of 0 and 1's representing 
%                             0 entries and variable entries respectively, or 
%                        (ii) a polynomial with variable coefficients being  
%                               nonzero, or 
%                       (iii) a cell array of matrices in (i) and/or (ii) 
%                       assuming the domain is a product of matrix spaces and 
%                       polynomial spaces 
%             para  -- (cell) parameters needed for running L 
%                       **** L must run with >> L(Domain{:},para{:}) *** 
%                r  -- (integer, optional) if provided, the rank-r TSVD 
%                       solution and kernel is computed (In other words,  
%                       substitute L and b with rank-r projections) 
%                b  -- (matrix, polynomial or cell) The right-hand side of 
%                           the equation L(z) = b, must be in the codomain  
%                           of L 
%               tol -- (numeric) error tolerance (< 1), 
%                         or r = tol > 1 indicating the rank-r projection of 
%                             the linear transformation is used 
%    
%      Output:   z -- the minimum norm solution so that ||L(z)-b|| is minimized 
%                K -- an orthonormal basis for the kernel {y | L(y) = 0} 
%              cnd -- condition number of the representation matrix. 
%              res -- (vector) the residuals of ||L(Z)-b||,||L(K{1}||,...||L(K{n}|| 
