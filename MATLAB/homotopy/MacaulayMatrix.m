%   Construct the Macaulay matrix of a polynomial ideal at an 
%   isolated zero
%    
%   Syntax:  >> M = MacaulayMatrix(z, P, var, depth);
% 
% Input:  z --- (1xn vector) the isolated zero
%         P --- (1xm cell) the polynomials generating the ideal
%       var --- (1xn cell) variable names of the ideal 
%     depth --- (integer) depth (differential order) of the Macaulay matrix
% 
% Output  M --- the Macaulay matrix in sparse format
% 
%   Example: For an approximate zero x = 0.57735, y = 0.57735 to the system
%       x^3+y-0.7698 = 0,     x+y^3-0.7698 = 0
%   the call sequence
% 
%     >> z = [0.57735,0.57735];  P = {'x^3+y-0.7698','x+y^3-0.7698'}; 
%     >> var = {'x','y'};  depth = 3; 
%     >> M = MacaulayMatrix(z, P, var, depth)
%   
%   generates the Macaulay matrix of depth 3.
% 
