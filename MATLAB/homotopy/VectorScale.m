%
%   Calculate the scalar multiplication of a scalar t and vector u 
%
%  Syntax:  >> v = VectorScale(s,u)
%
%   Input:  t --- (numeric) the scalar
%           u --- (matrix/string/cell) a general vector that is either
%                (i) an mxn matrix, or
%                (ii) a polynomial in character string, or
%                (iii) a 1xk cell array, each u{j} is either (i) or (ii)
%
%  Example:
% >> v = VectorScale(3,{[1 2 3; 4 5 6], '7+8*x*y+9*z^3'});
% >> v{:}
% ans =
%      3     6     9
%     12    15    18
% ans =
% 21 + 24*x*y + 27*z^3
%
