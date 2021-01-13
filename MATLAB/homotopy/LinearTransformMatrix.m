%function [M, R] = LinearTransformMatrix(F, Domain, para)
%
%  Generic routine for generating the linear transformation matrix
%
%    Syntax:  >> M = LinearTransformMatrix(@F, Domain, para)
%             >> [M,R] = LinearTransformMatrix(@F, Domain, para)
%
%    Input:  F  -- Matlab function name for calculate the linear
%                  transformation with syntax  F(Z,Q1,...,Qk) where Z is 
%                       (i) a matrix (including vector), or
%                      (ii) a polynomial, or
%                     (iii) a cell array of matrices and/or polynomials
%                    to be transformed by the linear transformation.
%                    The remaining parameters Q1,...,Qk must be 
%                    provided as cell entries for para.
%        Domain  -- domain of the linear transformation represented by
%                      (i) an mxn matrix of 0' and 1's representing
%                          0 entries and variable entries respectively, or
%                     (ii) a polynomial with variable coefficients being 
%                            nonzero, or
%                    (iii) a cell array of matrices in (i) and/or (ii)
%          para  -- (cell) parameters needed for running F
%
%   Output:  M  -- the (sparse) matrix for the linear transformation
%            R  -- (cell) codomain of the linear transformation, each cell
%                         entry contains 
%                      (i) an mxn matrix of 0' and 1's representing
%                          zero entries and variable entries respectively,
%                     (ii) or, a 3-dimensional cell with variable names, tuple
%                     degree and term indices for a polynomial component
%
