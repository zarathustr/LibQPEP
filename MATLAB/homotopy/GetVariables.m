%
% Extract variable names from any number of polynomial strings
%
% Syntax:  (pvar is the shortened alias for GetVariables)
%          >> z = GetVariables(f)
%          >> z = GetVariables(f,g,h)
% 
%  to extract variables from a cell array f of polynomial strings:
%          >> z = GetVariables(f{:})
%
%  Input:   f --- (string) polynomial
% Output:   z --- (cell)   variable names of f
% 
% Example:  >> f = '-2*x + 9*x*y^2 - 9*x^2*y^2 + 8*x^3*y^2';
%           >> g = '5 - x^3*z';
%           >> z = GetVariables(f,g)
%           z = 
%               'x'    'y'    'z'
