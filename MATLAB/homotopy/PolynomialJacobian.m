%      
% Calculate the Jacobian of a polynomial system
%      
% Syntax:  (pjac is the shortened alias for PolynomialJacobian)
%   >> J = PolynomialJacobian(F,var)  
%   >> J = PolynomialJacobian(F,var,z)
%      
%  Input:   F --- (cell) polynomial system as a cell array of strings
%         var --- (string or cell)    variable name(s)
%           z --- (vector, optional) variable values at which
%                      the Jacobian is calculated
%      
% Output:   J --- Jacobian of F 
%                 If z is not provided, J will be a cell array
%                 If z is provided, J will be a matrix
%      
% Example:  >>  F = {'x + y^2-3','x^2+z^4*y+5','x*z-6'}
%           >>  J = PolynomiaJacobian(F,{'x','y','z'})
%      
%           J = 
%      
%               '1'      '2*y'    '0'      
%               '2*x'    'z^4'    '4*z^3*y'
%               'z'      '0'      'x'      
%      
%           >> J = PolynomiaJacobian(F,{'x','y','z'}, [1,2,3])
%      
%           J =
%      
%                1     4     0
%                2    81   216
%                3     0     1
