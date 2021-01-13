% 
%  Sort the terms of a multivariate polynomial f according to  
%  the lexicographical order, and combine the like terms 
% 
%    Syntax:  >> G = mvPolynSort(F) 
%             >> G = mvPolynSort(F, tord) 
%             >> [G,idg] = mvPolynSort(F) 
%             >> [G,idg] = mvPolynSort(F, tord) 
% 
%    Input:   F -- (matrix) multivariate polynomial in coeff. matrix 
%          tord -- (string, optional) 'lex' or 'grlex', term order type 
%                   The default is 'lex' if not provided 
%
%   Output:   G --- (matrix) sorted coeff. matrix of the input polynomial  
%           idg --- (vector) the indices of G terms 
% 
%   Examples: 
% 
%   >> F  
%       
%   F = 
%     
%     0     0     4     3     3     0 
%     0     0     1     0     0     1 
%     0     0     1     0     0     2 
%     8     7    -7     2     4    -5 
%       
%   >> [g,idg] = mvPolynSort(F) 
% 
%   G = 
%  
%     0     3     4     0 
%     0     0     1     1 
%     0     0     1     2 
%    15     6    -7    -5 
% 
%   idg = 
% 
%     1     4    20    26 
% 
%   >> [G,idg] = mvPolynSort(F,'grlex') 
% 
%   G = 
% 
%     0     0     3     4 
%     0     1     0     1 
%     0     2     0     1 
%    15    -5     6    -7 
% 
% 
%   idg = 
% 
%     1    12    20    80 
% 
