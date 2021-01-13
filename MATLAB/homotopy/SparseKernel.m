%  Computer the numerical kernel within a given threshold 
%  of a sparse matrix
%
%  Syntax:
%           K = SparseKernel(A,thres)
%           [K,R] = SparseKernel(A,thres)
%           [K,R,s] = SparseKernel(A,thres)
%
%  Input:       A --- the matrix, sparse or full
%           thres --- the rank threshold. Any singular value below thres 
%                        is considered zero
%
%  Output:      K --- the unitary matrix whose columns form an orthonormal
%                        basis for Kernel(A)
%               R --- the unitary matrix whose columns form an orthonormal
%                        basis for Corange(A)
%               s --- the vector of the estimated singular values
