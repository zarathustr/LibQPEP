% 
%   Generic Gauss-Newton iteration routine for the least squares solution of 
%               F(z) = 0 
%   for z in a vector space domain 
%   
%     There are two versions:  
%          (i) generic version: solution consists of components of vectors,  
%                matrices and polynomials. 
%         (ii) column vector version: Iterates on column vectors 
%   
%     Syntax of version (i)  i.e. generic version: 
%         [z,res,fcond] = GaussNewton({@F,domain,para}, @G, z0)  
%         [z,res,fcond] = GaussNewton({@F,domain,para}, @G, z0, trac)  
%         [z,res,fcond] = GaussNewton({@F,domain,para}, @G, z0, trac, tol) 
%         [z,res,fcond] = GaussNewton({@F,domain,para}, {@G,r}, z0)  
%         [z,res,fcond] = GaussNewton({@F,domain,para}, {@G,r}, z0, trac)  
%         [z,res,fcond] = GaussNewton({@F,domain,para}, {@G,r}, z0, trac, tol) 
%   
%     Input:   F -- (function) Matlab function name for calculating y = F(z) 
%                    with syntax >> y = F(z1,...zk,p1,...,pm) where z1,...,zk 
%                    are components of z and p1,...,pm are fixed parameters. 
%                    (F must run on F(domain{:},para{:}). 
%         domain -- domain of the homomorphism represented by 
%                        (i) an mxn matrix of 0 and 1's representing 
%                            0 entries and variable entries respectively, or 
%                       (ii) a polynomial with variable coefficients being  
%                              nonzero, or 
%                      (iii) a cell array of matrices in (i) and/or (ii) 
%                      assuming the domain is a product of matrix spaces and 
%                      polynomial spaces 
%            para -- (cell) parameters needed for running F 
%               G -- (function) Matlab function name for calculating the  
%                     Jacobian u = F_z(z0)y F as a homomorphism with syntax 
%                        >> u = G(y1,...,yk,z1,...zk,p1,...,pm)  
%                     where z1,...,zk,p1,...,pk are the same for F but 
%                     considered fixed, the variables are y1,...,yk in the 
%                     same domain of F.  
%               r -- (integer, optional) if provided, the rank-r projection 
%                       of the Jacobian is used for the iteration 
%              z0 -- (cell/matrix/polynomial) the initial iterate for the 
%                    Gauss-Newton iteration in the domain 
%            trac -- (0,1,or 2, optional) flag for tracking the iteration 
%                     if tracking = 0 or missing, no tracking 
%                     if tracking = 1, track the first component z1 of z 
%                     if tracking = 2, track all components of z 
%              tol -- (numeric, optional) error tolerance of the iterate 
%                     if missing, iteration stops when residule stops 
%                     decreasing 
%   
%   Output:         z  -- the least squares solution in the domain 
%                 res  -- (optional) the residual ||F(z)|| 
%               fcond  -- (optional) the condition number 
%   
%   Syntax of version (ii), i.e. column vector version:   
%       [z,res,fcond] = GaussNewton(@Func, @Jacb, z0, {p1,...,pn},trac) 
%   
%   Input:      Func  -- Matlab function name for F(z) 
%               Jacb  -- Matlab function name for the Jacobian J(z) of F(z) 
%                         both Func and Jacb must be written to accept input 
%                         z0 along with other parameters p1, p2, ..., pn 
%                 z0  -- the initial iterate vector 
%     {p1,p2,...,pn}  -- (cell) vector parameters for Func and Jacb 
%               trac  -- tracking flag, 0, or 1 
%   
%   Output:         z  -- the least squares solution vector 
%                 res  -- (optional) the residual ||F(z)|| 
%               fcond  -- (optional) the condition number 
