%
%
%  Compute the numerical kernel of a general homomorphism L(Z) = 0
%
%   Syntax: >> Z = NumericalKernel(A,tol)  %           >> Z = NumericalKernel({@H,Domain,para},tol)  %
%    Input: H -- (cell/matrix) the homomorphism cell or a matrix
%                For a general homomorphism, H must be a cell 
%                         {@F,Domain,para}   
%                consists of the input for HomomorphismMatrix. Refer to
%                the document of HomomorphismMatrix for its input
%
%         tol -- (numeric, optional) the error tolerance. If not provided,
%                the default value of 1e-10*||H|| is used.
%
%    Example: For given polynomial f and g, find the cofactor pair (v,w) 
%       such that  f*v-g*w = 0.
%       Step 1:  Write a Matlab function for the homomorphism (v,w)->f*v-g*w
%          function p = SylvesterMap(v,w,f,g)
%          %          p = pminus(ptimes(f,v),ptimes(g,w));
%
%       Step 2: Execute the following
%          >> f = '2+5*x+6*x^2+4*x^3+x^4';  %          >> g = '6+7*x++2*x^2+2*x^3+x^4'; %          >> Domain = {'1+x+x^2','1+x+x^2'}; %          >> para = {f,g}  %          >>  Z = NumericalKernel({@SylvesterMap,Domain, para},1e-10)
%          Z = 
%            Column 1
%           '-0.801783725737273 + 0.267261241912424*x - 0.267261241912424*x^2'
%            Column 2
%           '-0.267261241912425 - 0.267261241912424*x - 0.267261241912425*x^2'
%            
