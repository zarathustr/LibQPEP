% Display a basis for dual space D computed by Multiplicity
%
%  Syntax:  ShowDual(D)
%
%  Input:   D --- (cell array) the basis of dual space computed by the
%                     function Multiplicity
%
% Output:   Screen display of the dual basis
%
%  Example:  >> [m,D,H] = Multiplicity({'x^2+sin(y^2)-2*x+1', 'x-cos(y)'}, ...
%                                      {'x','y'}, [1,0] )
%            >> ShowDual(D)
%                     d00 
% 
%                     d01 
% 
%                     0.447214*d10 -0.894427*d02 
%
%                     -0.816497*d03 +0.408248*d11 
%            
