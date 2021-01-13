%      
%  A squarefree test for polynomial f 
%
%  Syntax: >> [sf,e,d] = SquarefreeTest(f)
%          >> [sf,e,d] = SquarefreeTest(f,tol)
%
%  Input:    f -- (matrix) polynomial in a coefficient matrix
%          tol -- (numeric, optional) relative error tolerance on 
%                   coefficients, 1e-10 is the default value if not given
%   
%  Output   sf -- true of false on squarefree
%            e -- (vector) the multiplicities of squarefree factors
%            d -- (matrix) the j-th column is the tuple degree of 
%                     j-th squarefree factor 
%
