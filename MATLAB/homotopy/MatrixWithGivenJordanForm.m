%  Construct a matrix with given Jordan Canonical Form
%
%  Syntax:  There are several choices for either LHS or RHS
%     >>      LHS      =      RHS
%  ------------------- | --------------------------------------
%        A             |   MatrixWithGivenJordanForm(s,j)
%        [A,X]         |   MatrixWithGivenJordanForm(s,j,k)
%        [A,X,Y]       |      
%
%  Input:
%       s -- eigenvalues in each Jordan block
%       j -- size of each Jordan block
%       k -- (optional) number of (simple) eigenvalues of a kxk random matrix
%
% Output:   A -- a matrix whose eigenvalues and Jordan structure are defined by
%                   input s, j and k
%         X,Y -- matrices such that Y*A*X is the Jordan Canonical Form of A
%
%   Example: To generate a matrix with 
%                  eigenvalue  | Jordan block sizes
%                      2       |   4, 3, 2
%                      3       |   5, 1
%   use the following calling sequence:
%      >> [A,X,Y] = MatrixWithGivenJordanForm([2 2 2 3 3],[4 3 2 5 1])
%   and verify it by 
%      >> Y*A*X
