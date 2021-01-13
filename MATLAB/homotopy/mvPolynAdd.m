% 
%Compute the addition f+g of two multivariate polynomials f and g,  
%     or the linear combination  s*f + t*g  if scalars s and t are also given 
%   
%     Syntax:  >> h = mvPolynAdd(f,g) 
%              >> h = mvPolynAdd(f,g,s,t) 
%              >> h = mvPolynAdd(f,g,s,t,tord) 
%  
%     Input:   f, g -- multivariate polynomials in coeff. matrices 
%                      (assuming like terms have been combined 
%              s, t -- (optional) complex numbers 
%              tord -- (optional) if tord = 'grlex', the output terms will 
%                         be sorted in graded lexico order 
%  
%    Output:   multivariate polynomial  f+g,  or s*f + t*g 
% 
