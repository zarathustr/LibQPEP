%  Numerical solution of polynomial systems by hom4ps Homotopy method 
%   
%  Syntax:  >> [S, var] = psolve( P ) 
%      
%  Input:  P --- (cell array) the polynomial system to be solved
%  For example:  
%     >> P = {'-x^5+y^5-3*y-1','5*y^4-3','-20*x+y-z'}
% 
%  Output: S --- (matrix) numerical solutions (as columns of S)
%        var --- (cell array) the array of variables
%   For example, output
%      S = 
%         -0.8264 + 0.6004i  -0.6092 - 1.0165i
%         -0.8801 + 0.0000i  -0.0000 - 0.8801i
%         15.6482 -12.0086i  12.1831 +19.4496i
%
%       var = {'x','y','z'}
%           means there are two numerical solutions 
%       (x,y,z) = 
%       (-0.8264+0.6004i, -0.8801+0.0000i, 15.6482-12.0086i)
%       (-0.6092-1.0165i, -0.8801+0.0000i, 12.1831+19.4496i)
