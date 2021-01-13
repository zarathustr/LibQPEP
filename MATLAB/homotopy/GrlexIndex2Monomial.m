% 
%   Generate a i-th monomial in graded lexicographical order 
% 
%   Syntax:  m = GrlexIndex2Monomial(i, n) 
%            m = GrlexIndex2Monomial(i, n, d) 
% 
%   Input:   i -- index of the monomial 
%            n -- the number of variables 
%            d -- (optional)  
%                 Case 1. (default) If not provided, the monomial is the  
%                          i-th one with all degrees allowed 
%                 Case 2:  If d is provided, the monomial is the i-th one 
%                          among monomials of degree d 
% 
%  Output:   m -- (column vector) exponents of the monomial 
% 
%  Example:  
% >> GrlexIndex2Monomial(18,3) 
%    ans = 
%          2 
%          0 
%          1 
% 
%  >> GrlexIndex2Monomial(5,4,3) 
%     ans = 
%          0 
%          1 
%          0 
%          2 
% 
