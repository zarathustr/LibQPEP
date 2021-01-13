%  compute the univariate numerical GCD of (p,q) given degree k
%
% Syntax: 
%
%  >> [u,v,w,res,cond] = uvGCDfixedDegree(f,g,k)        
%
% Input   f, g -- the polynomial pairs as coefficient vectors
%            k -- the tuple degree of the GCD
%
% Output     u,v,w -- the numerical GCD triplet such that 
%                         conv(u,v) approximates f, conv(u,w) approximates g
%              res -- the (weighted) residual
%             cond -- the numerical GCD condition number
%
% Example:
%
%  >> f = [1  3  6 10 10  9  7  4];
%  >> g = [1  6 10 13 15  15 10  6  3  1];
%  >> [u,v,w,res,cond] = uvGCDfixedDegree(f,g,4);
