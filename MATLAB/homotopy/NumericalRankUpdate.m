%  <Purpose>
%     Computing the numerical rank of a updated matrix.
%
%  <Syntax>
%     [r,Basis,C] = NumericalRankUpdate(A,pth,vec,C,RC)
%
%  <Input Parameters>
%   1.   A -- the target matrix;
%   2. pth -- the index of row/column to be inserted;
%   3. vec -- the row/column vector to be inserted;
%   4.   C -- cell array contains information required by updating/downdating;
%   5.  RC -- Set to 'row', then the pth row will be inserted.
%             Set to 'column', then the pth column will be inserted.
%
%  <Output Parameters>
%   1. r     -- the numerical rank of the updated matrix;
%   2. Basis --
%          For high rank cases...
%              a matrix whose columns form an orthonormal basis of
%              the numerical kernel;
%
%          For low rank cases...
%              a matrix whose columns form an orthonormal basis of
%              the numerical range;
%
%   3.     C -- Matlab cell array
%
%      For high rank cases...
%          C{1,1} = rank : the numerical rank of the updated matrix;
%          C{2,1} = Basis : matrix whose columns form an orthonormal kernel basis;
%          C{3,1} = Q : the Q in the QR decomposition of the kernel stacked matrix;
%          C{4,1} = R : the R in the QR decomposition of the kernel stacked matrix;
%          C{5,1} = tau : scaling factor in the kernel stacked matrix;
%          C{6,1} = tol : the rank decision threshold;
%
%      For low rank cases...
%          C{1,1} = rank : the numerical rank of the updated matrix;
%          C{2,1} = U : the U in the USV+E decomposition of the updated matrix;
%          C{3,1} = V : the V in the USV+E decomposition of the updated matrix;
%          C{4,1} = S : the S in the USV+E decomposition of the updated matrix;
%          C{5,1} = tol : the rank decision threshold;
%
%  <Reference>
%    [1] T.Y. Li and Z. Zeng, "A Rank-Revealing Method with Updating, Downdating
%        and Applications", SIAM J. Matrix Anal. and Appl., 26 (2005), pp. 918--946.
%
%    [2] T.L. Lee, T.Y. Li and Z. Zeng, "A Rank-Revealing Method with Updating, 
%        Downdating and Applications, Part II", SIAM J. Matrix Anal. and Appl. (2009).
