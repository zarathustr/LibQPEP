%
% uvGCD computes the numerical greatest common divisor of a given
% polynomial pair (f,g) within a given threshold tol. 
%
%
% Input   f, g -- the polynomial pairs as coefficient vectors
%          tol -- the residual tolerance
%       wtflag -- the weight flag
%                   wtflag = 0: weight norm is not used
%                   wtflag = 1: using weight 2-norm
%
% Output     u,v,w -- the numerical GCD triplet such that 
%                         conv(u,v) approximates f, conv(u,w) approximates g
%              res -- the (weighted) residual
%             cond -- the numerical GCD condition number
%
% Syntax: 
%
%  >> [u,v,w,res,cond] = uvgcd(f,g,tol)    % default
%  >> [u,v,w,res,cond] = uvgcd(f,g,tol,0)  % disable weighted norm
%
% Example:
%
%  >> f = [1  3  6 10 10  9  7  4];
%  >> g = [1  6 10 13 15  15 10  6  3  1];
%  >> [u,v,w,res,cond] = uvgcd(f,g,1.0d-10);
%
