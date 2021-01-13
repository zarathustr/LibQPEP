%  MULTIPLICITY --> Computing the multiplicity and multiplicity structure
%  of a system of nonlinear equations at an isolated zero.
%
%  <Synopsis>
%          m = Multiplicity(f, variables, zero)
%          m = Multiplicity(f, variables, zero, options)
%  [m, D, H] = Multiplicity(f, variables, zero)
%  [m, D, H] = Multiplicity(f, variables, zero, options)
%
%  <Input Parameters>
%    1. f         --> a cell array containing the system of nonlinear
%                     equations as strings, e.g.,
%                         >>   f = { 'x^2 + sin(y) -1',  'x-cos(y)' };
%
%    2. variables --> a cell array containing the unknown variables as
%                     strings, e.g.,
%                         >>   variables = { 'x',  'y' };
%
%    3. zero      --> a vector containing numerical zero (root) of f, e.g.,
%                         >>   zero = [1, 0];
%
%    4. options   --> an optional parameter which includes the configuration
%                     settings:
%                     Display: Set to 2 to have all output (the dual space
%                              and Hilbert function) printed to the screen,
%                              and set to 1 to have depth and Hilbert
%                              function printed to the screen. Otherwise
%                              set the default value 0;
%                         Tol: The threshold for numerical rank-revealing.
%                              Singular values above Tol are counted as
%                              nonzero. The default value is 1e-8;
%                     EqnType: The equation type for MULTIPLICITY,
%                              polynomial system ('Poly') or nonlinear
%                              functions ('Nonlinear'). The default value is
%                              'Nonlinear'. By setting the value to 'Poly',
%                              MULTIPLICITY will transfer the polynomial
%                              system to the matrix representation
%                              internally and speed up the computation.
%                      MaxInt: Maximum multiplicity allowed in the recursive
%                              computation. If a zero is not isolated, the
%                              multiplicity will be infinity.  The code can
%                              be used for identifying a nonisolated zero by
%                              setting MaxInt to a known upper bound (e.g.
%                              the Bezout number). The default value is 1000.
%
%                     All the configuration settings may be changed by using
%                     optset function, i.e.,
%                         >>   options = optset('para', value);
%                     para could be any configuration setting, value is set
%                     to para. See OPTSET for details. Any configuration
%                     setting that is not changed will be set to its default
%                     value.
%
%  <Output Parameters>
%    1. m        --> the multiplicity of f at the zero;
%
%    2. D        --> a basis for the dual space as a cell array with
%                    each component being a matrix in the Matlab form of
%                                   D{i} = [c_1, j_1;
%                                           c_2, j_2;
%                                             ...
%                                           c_n, j_n ];
%                    representing a differential functional
%                            D_i = c_1*d_{j_1} + ... + c_n*d_{j_n}
%                    where d_{j_i}'s are differential monomial functionals
%                    (e.g. For a system of equations with variables {x,y,z}
%                     at the zero x=a, y=b, z=c, the functional d_{[i,j,k]}
%                     applied to any function g is the value of the partial
%                     derivative
%
%                                                         i+j+k
%                                             1          d
%                         d_{[i,j,k]}(g) = -------- * ----------- g(a,b,c)
%                                          i! j! k!     i   j   k
%                                                     dx  dy  dz
%
%                     The dual space is the vector space consists of such
%                     differential functionals that vanish on the nonlinear
%                     system while satisfying the so-called closedness
%                     condition);
%
%    3. H        --> values of the Hilbert function in a vector.
%
%  <Examples>
%   Consider the nonlinear system
%
%               sin(x)*cos(y)-x       = 0
%               sin(y)*sin(x)^2 - y^2 = 0
%
%   at the zero (x, y) = (0, 0), the multiplicity can be computed by the
%   following statements:
%
%     >> f = {'sin(x)*cos(y) - x', 'sin(y)*sin(x)^2 - y^2'};
%     >> variables = {'x', 'y'};
%     >> zero = [0, 0];
%     >> m = Multiplicity(f, variables, zero)
%
%   To create an options structure with Tol = 1e-10, MaxInt = 100:
%     >> options = optset('Tol', 1e-10, 'MaxInt', 100);
%     >>  m = Multiplicity(f, variables, zero, options)
%
%  <Algorithm>
%   This code implements a modified closedness subspace method for
%   multiplicity identification with a newly developed equation-by-equation
%   strategy for improving efficiency.
%
%  <References>
%  [1] An algorithm and software for computing multiplicity structures at
%      zeros of nonlinear systems, W. Hao, A. J. Sommesse and Z. Zeng
%  [2] The closedness subspace method for computing the multiplicity
%      structure of a polynomial system, Z. Zeng.
%
%  See also optset, cell
