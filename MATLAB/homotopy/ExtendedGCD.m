%
% ExtendedGCD computates the extended GCD in approximate sense
% For a given polynomial pair (f,g), it computes a polynomial pair (p,q)
% such that 
%         p*f + q*g = GCD(f,g)  in approximate sense
%
%  Calling syntax:  >> [p, q] = ExtendedGCD(f, g, u);
%
%  Input   f, g -- coefficient vectors of f and g, respectively
%  Input      u -- coefficient vector of GCD(f,g)
%                    >>[u,v,w] = uvGCD(f,g,tol)
%
%  Output  p, q -- polynomials such that p*f+q*g = GCD(f,g) approximately
%
