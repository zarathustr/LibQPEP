% mvGCDfixedDegree computes the nearest greatest common divisor of a 
% multivariate polynomial pair (F,G) with a fixed GCD degree. The code is 
% optimized under the assumption that  F, G and the GCD triplet are sparse 
% fewnomials
%
% Syntax:   >> [U,V,W,res,cond] = mvPolynGCD(f,g,gcddeg) 
%
% Input   F, G -- the polynomial pairs as coefficient matrices 
%       gcddeg -- the tuple degree of the GCD, if known
%
% Output     u,v,w -- the NGCD triplet such that 
%                     u*v approximates f, u*w approximates g
%              res -- the (weighted) residual
%             cond -- the AGCD condition number
%
% Example:
%      >> F
%
%      F =
%
%           0     2     1     3
%           0     1     1     2
%           0     0     3     3
%           3    -6     1    -2
%
%      >> G
%
%      G =
%
%           0     1     1     2
%           0     0     1     1
%           0     2     3     5
%           3     9     1     3
%
%      >> [u,v,w,res,cond] = mvGCDfixedDegree(F,G,[1,1,3])
%
%      u =
%
%               0    1.0000
%               0    1.0000
%               0    3.0000
%        -11.6190   -3.8730
%
%
%      v =
%
%               0    2.0000
%               0    1.0000
%               0         0
%         -0.2582    0.5164
%
%
%      w =
%
%               0    1.0000
%               0         0
%               0    2.0000
%         -0.2582   -0.7746
%
%
%      res =
%
%        3.0461e-034
%
%
%      cond =
%
%          0.8969
%          
