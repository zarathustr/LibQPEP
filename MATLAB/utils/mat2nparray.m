% 
% LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems (QPEPs),
%          It also gives highly accurate uncertainty description of the solutions.
%
%
% Article: 
%      Wu, J., Zheng, Y., Gao, Z., Jiang, Y., Hu, X., Zhu, Y., Jiao, J., Liu, M. (2020)
%           Quadratic Pose Estimation Problems: Unified Solutions, 
%           Solvability/Observability Analysis and Uncertainty Description 
%           in A Globally Optimal Framework.
%
%
% Authors:      Jin Wu and Ming Liu
% Affiliation:  Hong Kong University of Science and Technology (HKUST)
% Emails:       jin_wu_uestc@hotmail.com; eelium@ust.hk
% Websites:     https://zarathustr.github.io
%               https://ram-lab.com
%         


function result = mat2nparray( matarray )
              %mat2nparray Convert a Matlab array into an nparray
              %   Convert an n-dimensional Matlab array into an equivalent nparray  
              data_size=size(matarray);
              if length(data_size)==1
                  % 1-D vectors are trivial
                  result=py.numpy.array(matarray);
              elseif length(data_size)==2
                  % A transpose operation is required either in Matlab, or in Python due
                  % to the difference between row major and column major ordering
                  transpose=matarray';
                  % Pass the array to Python as a vector, and then reshape to the correct
                  % size
                  result=py.numpy.reshape(transpose(:)', int32(data_size));
              else
                  % For an n-dimensional array, transpose the first two dimensions to
                  % sort the storage ordering issue
                  transpose=permute(matarray,[length(data_size):-1:1]);
                  % Pass it to python, and then reshape to the python style of matrix
                  % sizing
                  result=py.numpy.reshape(transpose(:)', int32(fliplr(size(transpose))));
              end
          end