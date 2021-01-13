%  generate the k-th Sylverster matrix of f and g
%  The k-th Sylvester matrix is the matrix of the linear transformation
%     L(v,w) = f*w+g*v
%  for deg(v) = deg(f)-k,  deg(w) = deg(g)-k.  
%  The default value for k is k=1 if not provided, which is the most
%   well-known Sylvester matrix.
%
%  Syntax:   >> S = SylvesterMatrix(f,g)
%            >> S = SylvesterMatrix(f,g,k)
%
%  Input:   f, g --- (strings or coefficient vectors) polynomials
%              k --- (integer, optional) see above
%
%  Output:     S --- (matrix) The Sylvester matrix
%
%  Example:  >> SylvesterMatrix('1+3*x-6*x^3','2+5*x^2')
%             ans =
%             -6     0     5     0     0
%              0    -6     0     5     0
%              3     0     2     0     5
%              1     3     0     2     0
%              0     1     0     0     2
%
%           >> SylvesterMatrix('1+3*x-6*x^3','2+5*x^2',2)
%
%         ans =
%
%             -6     5     0
%              0     0     5
%              3     2     0
%              1     0     2
