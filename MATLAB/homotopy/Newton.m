%   Generic rank-r Newton's iteration routine solving    F(z) = 0  where 
%   F is a smooth mapping between vector spaces, assuming the Jacobian 
%   is of rank r at the solution that can be non-isolated 
%   (c.f. http://homepages.neiu.edu/~zzeng/Papers/Rank-r_Newton.pdf) 
%   
%     Syntax  
%         [z,res,fcond] = Newton({@F,domain,para}, @G, z0)  
%         [z,res,fcond] = Newton({@F,domain,para}, @G, z0, trac)  
%         [z,res,fcond] = Newton({@F,domain,para}, @G, z0, trac, tol) 
%         [z,res,fcond] = Newton({@F,domain,para}, {@G,r}, z0)  
%         [z,res,fcond] = Newton({@F,domain,para}, {@G,r}, z0, trac)  
%         [z,res,fcond] = Newton({@F,domain,para}, {@G,r}, z0, trac, tol) 
%   
%     Input:   F -- (function) Matlab function name for calculating y = F(z) 
%                    with syntax >> y = F(z1,...zk,p1,...,pm) where z1,...,zk 
%                    are components of z and p1,...,pm are fixed parameters. 
%                    (F must run on F(domain{:},para{:}). 
%         domain -- domain of the homomorphism represented by 
%                        (i) an mxn matrix of 0 'and 1's representing 
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
%                    Newton's iteration in the domain 
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
