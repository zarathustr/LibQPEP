% estimateWorldCameraPose Estimate camera pose from 3-D to 2-D point correspondences
%   [worldOrientation, worldLocation] = estimateWorldCameraPose(imagePoints, worldPoints, cameraParams)
%   returns the orientation and location of a calibrated camera in the world coordinate
%   system in which worldPoints are defined. 
% 
%   The function solves the Perspective-n-Point (PnP) problem using the P3P algorithm. 
%   The function eliminates spurious correspondences using the M-estimator SAmple Consensus 
%   (MSAC) algorithm
% 
%   Inputs            Description
%   ------            -----------
%   imagePoints       M-by-2 array of [x,y] coordinates of undistorted image points,
%                     with M >= 4.  
%  
%   worldPoints       M-by-3 array of [x,y,z] coordinates of world points.
%  
%   cameraParams      a cameraParameters or cameraIntrinsics object
%  
%   Outputs           Description
%   -------            -----------
%   worldOrientation  orientation of the camera in the world coordinates specified
%                     as a 3-by-3 matrix
%   worldLocation     location of the camera in the world coordinates specified as a
%                     1-by-3 vector
%   
%   [..., inlierIdx] = estimateWorldCameraPose(...) additionally returns the indices of the 
%   inliers used to compute the camera pose. inlierIdx is an M-by-1 logical vector, 
%   with the values of true corresponding to the inliers.
% 
%   [..., status] = estimateWorldCameraPose(...) additionally returns a status code. If the 
%   status output is not specified, the function will issue an error if the number of 
%   input points or the number of inliers is less than 4. The status can have the following 
%   values:
%  
%       0: No error.
%       1: imagePoints and worldPoints do not contain enough points.
%          At least 4 points are required.
%       2: Not enough inliers found. At least 4 inliers are required.
% 
%   [...] = estimateWorldCameraPose(..., Name, Value) specifies additional name-value pair 
%   arguments described below:
% 
%   'MaxNumTrials'          Positive integer scalar.
%                           Specifies the number of random trials for finding
%                           the outliers. The actual number of trials depends
%                           on imagePoints, worldPoints, and the values of the 
%                           MaxReprojectionError and Confidence parameters. Increasing 
%                           this value will improve the robustness of the output 
%                           at the expense of additional computation.
%  
%                           Default: 1000
%  
%   'Confidence'            Scalar value greater than 0 and less than 100.
%                           Specifies the desired confidence (in percentage)
%                           for finding the maximum number of inliers. Increasing 
%                           this value will improve the robustness of the output 
%                           at the expense of additional computation.
%  
%                           Default: 99
%  
%   'MaxReprojectionError'  Positive numeric scalar specifying the reprojection
%                           error threshold in pixels for finding outliers. Increasing 
%                           this value will make the algorithm converge faster, 
%                           but may reduce the accuracy of the result.
%  
%                           Default: 1
%  
%    Notes
%    -----
%    The function does not account for lens distortion. You can either
%    undistort the images using the undistortImage function before 
%    detecting the image points, or you can undistort the image points 
%    themselves using the undistortPoints function.
% 
%    Class Support
%    -------------
%    imagePoints and worldPoints must be of the same class, which can be 
%    double or single. cameraParams must be a cameraParameters or
%    cameraIntrinsics object. orientation and location are of class double.
% 
%    Example: Determine camera pose from world-to-image correspondences
%    ------------------------------------------------------------------
%    data = load('worldToImageCorrespondences.mat');
%    [worldOrientation, worldLocation] = estimateWorldCameraPose(...
%        data.imagePoints, data.worldPoints, data.cameraParams);
%    pcshow(data.worldPoints, 'VerticalAxis', 'Y', 'VerticalAxisDir', 'down', ...
%        'MarkerSize', 30);
%    hold on
%    plotCamera('Size', 10, 'Orientation', worldOrientation, 'Location',...
%        worldLocation);
%    
%    See also relativeCameraPose, viewSet, triangulateMultiview, bundleAdjustment, 
%             plotCamera, pcshow, extrinsics, triangulate,
%             cameraPoseToExtrinsics

%  Copyright 2015 The MathWorks, Inc.

% References:
% ----------
% [1] X.-S. Gao, X.-R. Hou, J. Tang, and H.-F. Cheng, "Complete Solution 
%     Classification for the Perspective-Three-Point Problem", IEEE Trans. 
%     Pattern Analysis and Machine Intelligence, vol. 25, no. 8, pp. 930-943, 
%     2003.

%#codegen

function [orientation, location, inlierIdx, status, err] = ...
    estimateWorldCameraPoseNew(imagePoints, worldPoints, cameraParams, varargin)

[params, outputClass, imagePts, worldPts] = parseInputs(imagePoints, worldPoints, cameraParams, ...
    varargin{:});

if isa(cameraParams, 'cameraIntrinsics')
    cameraParams = cameraParams.CameraParameters;
end

% List of status codes
statusCode = struct(...
    'NoError',           int32(0),...
    'NotEnoughPts',      int32(1),...
    'NotEnoughInliers',  int32(2));

% Additional RANSAC parameters
params.sampleSize = 4;
params.recomputeModelFromInliers = false;
params.defaultModel.R = nan(3);
params.defaultModel.t = nan(1, 3);

% RANSAC function handles
funcs.fitFunc = @solveCameraPose;
funcs.evalFunc = @evalCameraPose;
funcs.checkFunc = @check;

numPoints = size(worldPts, 1);
if numPoints < params.sampleSize
   status = statusCode.NotEnoughPts;
   [orientation, location, inlierIdx] = badPose(numPoints);
else
    % Compute the pose using RANSAC
    points = pack(imagePts, worldPts);
    [isFound, pose, inlierIdx] = vision.internal.ransac.msac(...
        points, params, funcs, cameraParams.IntrinsicMatrix, outputClass);
    
%     isFound = true;
%     pose = solveCameraPose_world_image(worldPts, imagePts
    if isFound
        % Convert from extrinsics to orientation and location
        orientation =  pose.R';
        location    = -pose.t * pose.R';
        err = pose.err / numPoints;
        status = statusCode.NoError;
    else
        % Could not compute the pose
        status = statusCode.NotEnoughInliers;
        [orientation, location, inlierIdx] = badPose(numPoints);
        err = nan;
    end
end

if nargout < 4
    checkRuntimeStatus(statusCode, status);
end

%==========================================================================
% Check runtime status and report error if there is one
%==========================================================================
function checkRuntimeStatus(statusCode, status)
coder.internal.errorIf(status==statusCode.NotEnoughPts, ...
    'vision:points:notEnoughMatchedPts', 'imagePoints', 'worldPoints', 4);

coder.internal.errorIf(status==statusCode.NotEnoughInliers, ...
    'vision:points:notEnoughInlierMatches', 'imagePoints', ...
    'worldPoints');

%--------------------------------------------------------------------------
function [orientation, location, inlierIdx] = badPose(numPoints)
orientation = nan(3);
location = nan(1, 3);
inlierIdx = false(numPoints, 1);



function pose = solveCameraPose(points, varargin)
[worldPoints, imagePoints] = unpack(points);

intrinsicMatrix = varargin{1};
% N = size(worldPoints, 1);
% [R,T,Xc,best_solution] = efficient_pnp_gauss([worldPoints, ones(N, 1)], [imagePoints, ones(N, 1)], intrinsicMatrix', true);

% Get up to 4 solutions for the pose using 3 points
[Rs, Ts] = vision.internal.calibration.solveP3P(...
    imagePoints, worldPoints, intrinsicMatrix);

% Choose the best solution using the 4th point
p = [worldPoints(4, :), 1]; % homogeneous coordinates
q = imagePoints(4, :);
% pose_pnp = chooseBestSolution(p, q, R', T', intrinsicMatrix);

pose.R = nan(3);
pose.t = nan(1, 3);
if ~isempty(Rs)
    pose = chooseBestSolution(p, q, Rs, Ts, intrinsicMatrix);
end


function pose = solveCameraPose_world_image(worldPoints, imagePoints, K)

intrinsicMatrix = K;
% N = size(worldPoints, 1);
% [R,T,Xc,best_solution] = efficient_pnp_gauss([worldPoints, ones(N, 1)], [imagePoints, ones(N, 1)], intrinsicMatrix', true);

% Get up to 4 solutions for the pose using 3 points
[Rs, Ts] = vision.internal.calibration.solveP3P(...
    imagePoints, worldPoints, intrinsicMatrix);

% Choose the best solution using the 4th point
p = [worldPoints(4, :), 1]; % homogeneous coordinates
q = imagePoints(4, :);
% pose_pnp = chooseBestSolution(p, q, R', T', intrinsicMatrix);

pose.R = nan(3);
pose.t = nan(1, 3);
if ~isempty(Rs)
    pose = chooseBestSolution(p, q, Rs, Ts, intrinsicMatrix);
end
pose.err = 0;
% pose = pose_pnp;

%--------------------------------------------------------------------------
function pose = solveCameraPose_(points, varargin)
[worldPoints, imagePoints] = unpack(points);

intrinsicMatrix = varargin{1};
NN = size(imagePoints, 1);
bbb = [imagePoints, ones(NN, 1)] / (intrinsicMatrix);
rrr = worldPoints(:, 1 : 3);
% for i = 1 : NN
%     bbb(i, :) = bbb(i, :) ./ norm(bbb(i, :));
%     rrr(i, :) = rrr(i, :) ./ norm(rrr(i, :));
% end
% [RR, tt, ss, xs] = apnp_algebraic(bbb, rrr, 0, 0, intrinsicMatrix', eye(4));

% Get up to 4 solutions for the pose using 3 points
[Rs, Ts] = vision.internal.calibration.solveP3P(...
    imagePoints, worldPoints, intrinsicMatrix);

% Choose the best solution using the 4th point
p = [worldPoints(4, :), 1]; % homogeneous coordinates
q = imagePoints(4, :);

pose.R = nan(3);
pose.t = nan(1, 3);
pose.err = 1;
if ~isempty(Rs)
    pose = chooseBestSolution(p, q, Rs, Ts, intrinsicMatrix);
end


function [R_, t_, s_, xs, min_val] = apnp_algebraic(bbb, rrr, Xw, U, K, Rt)
    N = size(bbb, 1);
    b = zeros(3, N);
    r = zeros(3, N);
    RR = Rt(1 : 3, 1 : 3);
    tt = Rt(1 : 3, 4);
    for i = 1 : N
        b(:, i) = [bbb(i, 1 : 3)]';
        r(:, i) = rrr(i, 1 : 3)';
%         
%         b(:, i) - RR * r(:, i) - tt
    end
    b_bar = zeros(3, 1);
    r_bar = zeros(3, 1);
    for i = 1 : N
        b_bar = b_bar + 1 / N * b(:, i); 
        r_bar = r_bar + 1 / N * r(:, i); 
    end
    
    P = zeros(4, 4);
    for i = 1 : N
        for j = 1 : 3
            str = sprintf('M = M%d_matrix(b(:, i) - b_bar);', j);
            eval(str);
            P = P + 1 / N * r(j, i) * M;
        end
    end
    [V, ~] = eig(P);

    q11 = V(:, 1); q_1 = q11 ./ norm(q11);
    q22 = V(:, 2); q_2 = q22 ./ norm(q22);
    q33 = V(:, 3); q_3 = q33 ./ norm(q33);
    q44 = V(:, 4); q_4 = q44 ./ norm(q44);
    x_count = 4 + 1;
    xs_ = [
        q_1';
        q_2';
        q_3';
        q_4';
        ];
    xs = zeros(4, 9);
    
    for i = 1 : x_count - 1
        q0 = xs_(i, 1);
        q1 = xs_(i, 2);
        q2 = xs_(i, 3);
        q3 = xs_(i, 4);
        R = quat2dcm([q0, q1, q2, q3]);
        aa = 0;
        bb = 0;
        for j = 1 : N
            aa = aa + 1 / N * (b(:, j)' * (b(:, j) - b_bar));
            bb = bb + 1 / N * (b(:, j)' * R * (r(:, j) - r_bar));
        end
        ss = bb / aa;
        s = ss;
        t = s * b_bar - R * r_bar;
        xs(i, :) = [q0, q1, q2, q3, t', s, ss];
    end
    
    Ls = zeros(x_count - 1, 1);
    for i = 1 : x_count - 1
        RR = quat2dcm(xs(i, 1 : 4));
        tt = xs(i, 5 : 7)';
        Ls(i) = J_func(RR, tt, bbb, rrr, K);
    end
    [minimum, idx] = sort(Ls);
    R_ = quat2dcm(xs(idx(1), 1 : 4));
    t_ = xs(idx(1), 5 : 7)';
    s_ = xs(idx(1), 9);
    min_val = minimum;


function J = J_func(R, t, b, r, K)
J = 0;
len = size(b, 2);
for j = 1 : len
    xb = b(j, :);
    bb = xb';
    xr = r(j, :)';
    xr = (R * xr + t);
    res = bb - xr;
    J = J + 1 / len * trace(res' * res);
end


%--------------------------------------------------------------------------
% Find the solution that results in smallest squared reprojection error for
% the 4th point.
% worldPoint must be in homogeneous coordinates (1-by-4)
function pose = chooseBestSolution(worldPoint, imagePoint, Rs, Ts, intrinsicMatrix)

pose.R = zeros(3);
pose.t = zeros(1, 3);

numSolutions = size(Ts, 1);
errors = zeros(numSolutions, 1, 'like', worldPoint);

for i = 1:numSolutions    
    cameraMatrix = [Rs(:,:,i); Ts(i,:)] * intrinsicMatrix;
    projectedPoint = worldPoint * cameraMatrix;
    projectedPoint = projectedPoint(1:2) ./ projectedPoint(3);
    d = imagePoint - projectedPoint;
    errors(i) = d * d';
end

[~, idx] = min(errors);
idx = idx(1); % in case we have two identical errors
pose.t = Ts(idx, :);
pose.R = Rs(:,:,idx);
pose.err = errors(idx);
    
%--------------------------------------------------------------------------
% Compute reprojection errors
function dis = evalCameraPose(pose, points, varargin)
[worldPoints, imagePoints] = unpack(points);

intrinsicMatrix = varargin{1};
cameraMatrix = [pose.R; pose.t] * intrinsicMatrix;

% Project world points into the image
numPoints = size(worldPoints, 1);
worldPointsHomog = [worldPoints, ones(numPoints, 1, 'like', worldPoints)];
projectedPointsHomog = worldPointsHomog * cameraMatrix;
projectedPoints  = bsxfun(@rdivide, projectedPointsHomog(:, 1:2), ...
    projectedPointsHomog(:, 3));

% Compute reprojection errors
diffs = imagePoints - projectedPoints;
dis = sum(diffs.^2, 2);

%--------------------------------------------------------------------------
% Pack points into a single entity for RANSAC
function points = pack(imagePoints, worldPoints)
points = [imagePoints, worldPoints];

%--------------------------------------------------------------------------
% Unpack the points
function [worldPoints, imagePoints] = unpack(points)
imagePoints = points(:, 1:2);
worldPoints = points(:, 3:end);

%--------------------------------------------------------------------------
function r = check(pose, varargin)
r = ~isempty(pose) && ~isempty(pose.R) && ~isempty(pose.t);

%--------------------------------------------------------------------------
function [ransacParams, outputClass, imagePts, worldPts] = ...
    parseInputs(imagePoints, worldPoints, cameraParams, varargin)

validatePoints(imagePoints, worldPoints);
imagePts = double(imagePoints);
worldPts = double(worldPoints);
outputClass = class(imagePts);
validateattributes(cameraParams, {'cameraParameters','cameraIntrinsics'}, ...
    {'scalar'}, mfilename, 'cameraParams');

defaults = struct('MaxNumTrials',1000, 'Confidence',99, 'MaxDistance', 1);

if isempty(coder.target)
    ransacParams = parseRANSACParamsMatlab(defaults, varargin{:});
else
    ransacParams = parseRANSACParamsCodegen(defaults, varargin{:});
end

%--------------------------------------------------------------------------
function ransacParams = parseRANSACParamsMatlab(defaults, varargin)
parser = inputParser;
parser.FunctionName = mfilename;
parser.addParameter('MaxNumTrials', defaults.MaxNumTrials, @checkMaxNumTrials);
parser.addParameter('Confidence', defaults.Confidence, @checkConfidence);
parser.addParameter('MaxReprojectionError', defaults.MaxDistance, @checkMaxDistance);

parser.parse(varargin{:});
ransacParams.confidence = parser.Results.Confidence;
ransacParams.maxDistance = parser.Results.MaxReprojectionError^2;
ransacParams.maxNumTrials = parser.Results.MaxNumTrials;
ransacParams.verbose = false;

%--------------------------------------------------------------------------
function ransacParams = parseRANSACParamsCodegen(defaults, varargin)
% Instantiate an input parser
parms = struct( ...
    'MaxNumTrials',       uint32(0), ...
    'Confidence',         uint32(0), ...
    'MaxReprojectionError',        uint32(0));

popt = struct( ...
    'CaseSensitivity', false, ...
    'StructExpand',    true, ...
    'PartialMatching', false);

% Specify the optional parameters
optarg = eml_parse_parameter_inputs(parms, popt, varargin{:});
ransacParams.maxNumTrials = eml_get_parameter_value(optarg.MaxNumTrials,...
    defaults.MaxNumTrials, varargin{:});
ransacParams.confidence   = eml_get_parameter_value(optarg.Confidence,...
    defaults.Confidence, varargin{:});
ransacParams.maxDistance  = eml_get_parameter_value(...
    optarg.MaxReprojectionError, defaults.MaxDistance, varargin{:});
ransacParams.verbose  = false;

checkMaxNumTrials(ransacParams.maxNumTrials);
checkConfidence  (ransacParams.confidence);
checkMaxDistance (ransacParams.maxDistance);

ransacParams.maxDistance = ransacParams.maxDistance^2;

%--------------------------------------------------------------------------
function validatePoints(imagePoints, worldPoints)
validateattributes(imagePoints, {'double', 'single'}, ...
    {'real', 'nonsparse', 'nonempty', '2d', 'ncols', 2}, ...
    mfilename, 'imagePoints');

validateattributes(worldPoints, {'double', 'single'}, ...
    {'real', 'nonsparse', 'nonempty', '2d', 'ncols', 3}, ...
    mfilename, 'worldPoints');

coder.internal.errorIf(~isa(imagePoints, class(worldPoints)), ...
    'vision:points:ptsClassMismatch', 'imagePoints', 'worldPoints');
coder.internal.errorIf(size(imagePoints, 1) ~= size(worldPoints, 1), ...
    'vision:points:numPtsMismatch', 'imagePoints', 'worldPoints');

%--------------------------------------------------------------------------
function tf = checkMaxNumTrials(value)
validateattributes(value, {'numeric'}, ...
    {'scalar', 'nonsparse', 'real', 'integer', 'positive'}, mfilename, ...
    'MaxNumTrials');
tf = true;

%--------------------------------------------------------------------------
function tf = checkConfidence(value)
validateattributes(value, {'numeric'}, ...
    {'scalar', 'nonsparse', 'real', 'positive', '<', 100}, mfilename, ...
    'Confidence');
tf = true;

%--------------------------------------------------------------------------
function tf = checkMaxDistance(value)
validateattributes(value,{'single','double'}, ...
    {'real', 'nonsparse', 'scalar','nonnegative','finite'}, mfilename, ...
    'MaxDistance');
tf = true;


