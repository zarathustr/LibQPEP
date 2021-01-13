%      
%  Generate the monomial basis 
% 
%  Syntax: B = mvPolynMonomBasis(d)               under tuple degree 
%          B = mvPolynMonomBasis(d, n)            n-variate, degree <= d 
%          B = mvPolynMonomBasis(d, n, tord)      as above, w/ grlex 
%          B = mvPolynMonomBasis(d, n, tord, dg)  as above, w/ fixed deg. 
%
%  Input:   d --- (integer or vector) total degree bound (if d is 1x1) or
%                    tuple degree bound (if d is a vector)
%           n --- (integer, optional) number of variables if d is total degree
%        tord --- (string, optional) term order,  'lex' or 'grlex'
%                  Case 1:  tord = 'lex' or missing, terms are in
%                              lexicographical order
%                  Case 2:  tord = 'grlex', terms are in graded
%                              lexicographical order
%           dg -- (string, optional) fixed degree activator
%                    If provided with dg = 'fixed', only monomials of fixed
%                    degree d are generated
%
%  Output:  B --- the coefficient matrix whose columns 
%                    represent the basis
%
%  Example:
% >> mvPolynMonomBasis([2,3])
%
% ans =
%
%     0     1     2     0     1     2     0     1     2     0     1     2
%     0     0     0     1     1     1     2     2     2     3     3     3
%     1     1     1     1     1     1     1     1     1     1     1     1
%
% >>mvPolynMonomBasis(3,2) 
% 
% ans = 
% 
%     0     1     2     3     0     1     2     0     1     0 
%     0     0     0     0     1     1     1     2     2     3 
%     1     1     1     1     1     1     1     1     1     1 
% 
% >> mvPolynMonomBasis(3,2,'grlex') 
% 
% ans = 
% 
%     0     0     1     0     1     2     0     1     2     3 
%     0     1     0     2     1     0     3     2     1     0 
%     1     1     1     1     1     1     1     1     1     1 
% 
