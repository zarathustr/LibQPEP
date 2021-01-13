% mvGCD computes the approximate greatest common divisor (AGCD) of a given
% multivariate polynomial pair (F,G). The code is optimized under the 
% assumption that  f, g and the GCD triplet are sparse fewnomials. 
%
% Input   F, G -- the polynomial pairs as coefficient matrices
%          tol -- (optional) the residual tolerance, if not sure, use zero 
%                    to invoke the default value   max(||f||,||g||)*1e-10
%
% Output     u,v,w -- the AGCD triplet such that 
%                       u*v approximates f, conv(u,w) approximates g
%              res -- the (weighted) residual
%             cond -- the AGCD condition number
%
% Syntax: 
%
%  >> [U,V,W,res,cond] = mvPolynGCD(f,g)   
%  >> [U,V,W,res,cond] = mvPolynGCD(f,g,tol)  
