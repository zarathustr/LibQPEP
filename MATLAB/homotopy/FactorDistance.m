% Calculate the distance between two factorizations 
%   Given two factorizations stored in 2-column cells F and G,
%   where the 1st column contains the polynomial factors and
%   the 2nd column entries are corresponding multiplicities,
%   FactorDistance calculate their distance that is independent
%   of scaling, choice of representative and order of factors.
% 
%  Input:      F --- (m x 2 cell array) 
%                      F{k,1}: the k-th factor of F as a string
%                      F{k,2}: the multiplicity of the k-th factor
%              G --- (n x 2 cell array)
%                     G{k,1}: the k-th factor of G as a string
%                     G{k,2}: the multiplicity of the k-th factor
%       permopt --- (string, optional)  'all':  use all permulations
%                         otherwise:  the current order only
%      
% Example:   
%      
%     f = 
%     
%         '-4888'            [1]
%         'x*y-2*x+3*y-1'    [3]
%         '3*y^2-5*x^3+4'    [6]
%     g = 
%     
%         '16'                          [1]
%         '352 + 264*y^2 - 440*x^3'     [6]
%         '-6 - 12*x + 18*y + 6*x*y'    [3]
%     
%     >> FactorDistance(f,g)
%     
%     ans =
%     
%       5.2663e-016
