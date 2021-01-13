% 
%  Find the index of a given monomial in the graded lexicographical order
% 
%  Syntax: >> idx = GrlexMonomialIndex(m) 
%          >> idx = GrlexMonomialIndex(m,d) 
% 
%  Input  m -- (vector) the exponents of the monomial 
%                        (e.g. [3;1;2] represents x^3*y*z^2) 
%         d -- (optional) total degree of the monomial 
%              Case 1: If d is not provided, the output index is the  
%                        location of the monomial in all possible terms 
%              Case 2: If d is provided, the output index is the location 
%                        of the monomial of degree d 
% 
% Output  idx -- the index of the monomial 
% 
% Example:   
%  >> GrlexMonomialIndex([2;0;1;1]) 
%    ans = 
%     
%        62 
%     
%    >> GrlexMonomialIndex([2;0;1;1],4) 
%     
%    ans = 
%     
%        27 
% 
