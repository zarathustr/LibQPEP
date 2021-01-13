% generate the Bezout matrix B_n (f,g) of polynomials f and g 
%   that is defined so that
%      
%     f(x)*g(y)-f(y)*g(x)
%     ------------------- = [1,x,...,x^(n-1)]*B_n(f,g)*[1,y,...,y^(n-1)]'
%            x - y
%      
% Syntax:   >> B = BezoutMatrix(f,g)
%           >> B = BezoutMatrix(f,g,n)
%      
% Input:   f, g --- (strings or coefficient vectors) polynomials
%             n --- (integer, optional) size of the Bezout matrix, 
%                      can not be smaller than either degrees of f or g
%                      default value is the larger degree of f and g
%      
% Output:     B --- (matrix) The Bezout matrix
%      
% Example:  >> BezoutMatrix('3*x^3-x','5*x^2+1')
%      
%           ans =
%      
%               -1     0     3
%                0     8     0
%                3     0    15
%      
%          >> BezoutMatrix('3*x^3-x','5*x^2+1',4)
%      
%          ans =
%      
%              -1     0     3     0
%               0     8     0     0
%               3     0    15     0
%               0     0     0     0
%      
