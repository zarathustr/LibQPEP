%  <Purpose>
%     Computing the numerical rank of a matrix.
%
%  <Syntax>
%         [r,Basis,C] = NumericalRank(A,tol,HL)
%
%  <Input Parameters>
%  1.   A -- the target matrix;
%  2. tol -- (optional) the rank decision threshold;
%                       default: tol = sqrt(n)*norm(A,1)*eps;
%  3   HL -- (optional)
%            Set to 'high rank' if A is a high rank matrix. (default)
%            Set to 'low rank' if A is a low rank matrix.
%
% <Output Parameters>
%  1. r     -- the numerical rank of A;
%  2. Basis -- (optional)
%
%              For high rank cases...
%                  a matrix whose columns form an orthonormal basis of
%                  the numerical kernel;
%
%              For low rank cases...
%                  a matrix whose columns form an orthonormal basis of
%                  the numerical range;
%
%  3.     C -- (optional) Matlab cell array contains information required
%                         by updating/downdating;
%
%      For high rank cases...
%          C{1,1} = rank : the numerical rank of A;
%          C{2,1} = Basis : matrix whose columns form an orthonormal kernel basis;
%          C{3,1} = Q : the Q in the QR decomposition of the kernel stacked matrix;
%          C{4,1} = R : the R in the QR decomposition of the kernel stacked matrix;
%          C{5,1} = tau : the scaling factor in the kernel stacked matrix;
%          C{6,1} = tol : the rank decision threshold;
%
%      For low rank cases...
%          C{1,1} = rank : the numerical rank of A;
%          C{2,1} = U : U in the USV+E decomposition of A;
%          C{3,1} = V : V in the USV+E decomposition of A;
%          C{4,1} = S : S in the USV+E decomposition of A;
%          C{5,1} = tol : the rank decision threshold;
%
% <Reference>
%   [1] T.Y. Li and Z. Zeng, "A Rank-Revealing Method with Updating, Downdating
%       and Applications", SIAM J. Matrix Anal. and Appl., 26 (2005), pp. 918--946.
%
%   [2] T.L. Lee, T.Y. Li and Z. Zeng, "A Rank-Revealing Method with Updating, 
%       Downdating and Applications, Part II", SIAM J. Matrix Anal. and Appl. (2009).
