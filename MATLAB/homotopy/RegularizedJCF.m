%  Computing an numerical Jordan decomposition
%                      X*J*inv(X)
%  of a given matrix A within a distance threshold theta, formulated under
%  a 'three-strikes principle'.  For details, see
%     Z. Zeng and T.Y. Li,  A numerical algorithm for computing the Jordan
%        Canonical Form, preprint, 2007
%
%  In a nutshell, let A be a matrix whose entries are given approximately 
%  with error magnitude being small.  The exact JCF of A will be degraded.
%  However, it is possible to recover the underlying Jordan structure and
%  approximate X and J by NumericalJoranForm
%
%  Syntax:  There are several choices for either LHS or RHS
%        >>   LHS      =      RHS
%  ------------------- | --------------------------------------
%   [J,X]              |   RegularizedJCF(A)
%   [J,X,e,s,t]        |   RegularizedJCF(A,theta)
%                      |   RegularizedJCF(A,theta,tau,gap)
%
% Input:    A -- matrix whose AJCF is to be computed
%       theta -- (optional) distance threshold
%         tau -- (optional) deflation threshold, simple eigenvalues
%                  whose geometric condition numbers above tau will be deflated
%         gap -- (optional) singular value gap used in rank revealing
%
% Output:   J -- The numerical Jordan Canonical Form
%           X -- The principle vector matrix
%           e -- the list of distinct eigenvalues
%           s -- The Jordan block sizes of the eigenvalues
%           t -- the residuals (backward errors) and the staircase condition
%                  numbers of the eigenvalues
