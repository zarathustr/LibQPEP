% Fast identification if a polynomial is squarefree
%      
%  Syntax: >> [sf,e,d] = SquarefreeTest(f)
%          >> [sf,e,d] = SquarefreeTest(f,tol)
%      
%  Input:    f -- (string) polynomial as a string
%          tol -- (numeric, optional) relative error tolerance on 
%                   coefficients, 1e-10 is the default value if not given
%      
%  Output   sf -- true of false on squarefree
%            e -- (vector) the multiplicities of squarefree factors
%            d -- (matrix) the j-th column is the tuple degree of 
%                     j-th squarefree factor
%      
% Example: For the polynomial (x*y+1)*(x^4+y^2+2)^3
%  >> f = ptimes('x*y+1',ppower('x^4+y+2',3));
%  >> [sf,e,d] = IsSquarefree(f)
%  sf =
%      0
%  e =
%      1     3
%  d =
%      1     4
%      1     2
