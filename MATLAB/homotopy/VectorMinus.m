%
% Calculate u-v for general vectors u and v in a product space
% 
% Syntax:  >> w = VectorMinus(u,v)
%
%  Input: u, v --- (matrix/string/cell) vectors in a product space given as
%                  (i) mxn matrices of the sime dimension
%                  (ii) polynomials in character strings
%                  (iii) cell array of (i) or (ii)
%
% Output: w --- (matrix/string/cell) the result of u-v stored in the same
%                 way as u and v
%
% Example:
% >> u = {[1 2 3; 4 5 6], '7+8*x*y+9*z^3'};
% >> v = {[1 1 1; 1 1 1],'1+x+x*y+z^3'};
% >> w  = VectorMinus(u,v);
% >> w{:}
% ans =
%
%     0     1     2
%     3     4     5
% ans =
% 6 - x + 7*x*y + 8*z^3
%
