%  generate the Companion matrix of the polynomials f 
%      
%  Syntax:   >> C = CompanionMatrix(f)
%      
%  Input:      f --- (string or coefficient vector) polynomial
%      
%  Output:     C --- (matrix) The companion matrix
%      
%  Example:  >> CompanionMatrix('x^4 + 2*x^3 + 3*x^2 + 4*x+5')
%      
%            ans =
%      
%                -2    -3    -4    -5
%                 1     0     0     0
%                 0     1     0     0
%                 0     0     1     0
