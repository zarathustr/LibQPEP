%
%  Convert a column vector to the vector in the domain D
%   (the inverse of the function RepresentationVector)
%  
%   Syntax:  v = RepresentationInverse(z, B)
%
%   Input:  z -- (vector, sparse) column vector in sparse format
%                    representing the vector v
%           B -- (matrix/string/cell) Basis info of the vector space:
%                (i) If B is an mxn matrix, the vector space consists of
%                    mxn matrices spanned by nonzero entries of B
%                (ii) If B is a character string, the vector space consists
%                    of polynomials spanned by the support of B
%                (iii) If B is a 1xm cell, the vector space is a
%                   a Cartesian product of m vector spaces, each component
%                   is either a matrix space or a polynomial space
%                 Specs of B:
%                 * length(B) = number of spaces in the Cartesian product
%                 * each cell entry B{j} contains one of the following
%                   (i) B{j} is an mxn matrix: the j-th component of u
%                      is an mxn matrix spanned by nonzero entries of B
%                   (ii) B{j} is a character string: the j-th component of
%                      u is a polynomial spanned by the support of B 
%                   (iii)B{j} is a 1x2 cell: the j-th component of u is
%                      a polynomial with
%                             B{j}{1} = cell of variable names such as
%                                          {'x','y','z'}
%                             B{j}{2} = tuple degree such as [3 2 5]
%                                 for degrees in each variable
%                         implying the complete monomial basis of the 
%                         variables under the tuple degree for u{j}
%                   (iv) a 1x3 cell with variable names, tuple degree,
%                         term indices for a polynomial component u{j}
%  Output:  v -- (matrix/string/cell) a (general) vector that is either
%                (i) an mxn matrix, or
%                (ii) a character string representing a polynomial, or
%                (iii) a cell array, each v{i} is either (i) or (ii) above
%
% Example:
%
% >> full(z).'   % ans =
%     0     2     3     1     0     4     5     0     0    -7     0     6
%
% >> u = RepresentationInverse(z, {ones(3,2),'1+x+x*y+x*y^2+x^3+x^5'});
% >> u{:}  % ans =
%      0     1
%      2     0
%      3     4 
% ans =
% 5 - 7*x^5 + 6*x*y^2
%
