%
% Convert a multivariate coefficient vector to coefficient matrix
%
%    Syntax:   
%      >> F = mvPolynCoefVec2Mat(f,d)  
%
%     Input:    f -- coefficient vector of a multivariate polynomial 
%               d -- tuple degree of the multivariate polynomial
%
%    Output:    F -- the coefficient matrix of the multivariate polynomial
%
%    Example:
%
%       >> f
%      
%      f =
%
%         (1,1)       3.0000
%         (3,1)       5.5000
%        (11,1)      -7.0000
%
%      >> mvPolynCoefVec2Mat(f,[2,3])
%
%      ans =
%
%               0    2.0000    1.0000
%               0         0    3.0000
%          3.0000    5.5000   -7.0000
