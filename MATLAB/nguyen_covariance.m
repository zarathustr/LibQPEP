% 
% LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems (QPEPs),
%          It also gives highly accurate uncertainty description of the solutions.
%
%
% Articles: 
%      Wu, J., Zheng, Y., Gao, Z., Jiang, Y., Hu, X., Zhu, Y., Jiao, J., Liu, M. (2022)
%           Quadratic Pose Estimation Problems: Globally Optimal Solutions, 
%           Solvability/Observability Analysis and Uncertainty Description.
%           IEEE Transactions on Robotics.
%           https://doi.org/10.1109/TRO.2022.3155880
%
%      Nguyen, H., Pham, Q. (2018) On the Covariance of AX = XB.                             
%           IEEE Transactions on Robotics.
%
%
% Authors:      Jin Wu and Ming Liu
% Affiliation:  Hong Kong University of Science and Technology (HKUST)
% Emails:       jin_wu_uestc@hotmail.com; eelium@ust.hk
% Websites:     https://zarathustr.github.io
%               https://ram-lab.com
%          
%
% nguyen_covariance.m: Comparison with Nguyen et al. 2018.
%
% Note: this file requires MATLAB version >= R2015a and Python 2.7,
%       Please modify pybinpath and systlibpath to adapt to your Python
%       executable and library searching path


clear all
close all
clc

warning('off');
format long g

% pybinpath = '/usr/local/bin/python2.7';
% systlibpath = '/usr/local/lib';

% pybinpath = '/usr/bin/python2.7';
% systlibpath = '/usr/lib';

pybinpath = '/usr/bin/python2.7';
systlibpath = '/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python';

str = sprintf('! export DYLD_LIBRARY_PATH=%s:DYLD_LIBRARY_PATH', systlibpath);
eval(str);

if(verLessThan('matlab', '8.5.0'))
   error('The MATLAB version is older than R2015a, please update.');
end

addpath('func_files');
addpath('solvers');
addpath('utils');
addpath('calib');
addpath('sdpt3');
run(fullfile('sdpt3', 'install_sdpt3.m'))
addpath('sedumi');
if(~verLessThan('matlab', '8.0.0'))
    run(fullfile('sedumi', 'install_sedumi.m'))
end
p = genpath('YALMIP');
addpath(p);

[env_, ~, isloaded] = pyversion;
if(~isloaded)
    if(verLessThan('matlab', '9.8.0'))
        pyversion(pybinpath);
    else
        pyenv('Version', pybinpath, 'ExecutionMode', 'InProcess');
    end
    insert(py.sys.path,int32(0), './cope')
end

cope = py.importlib.import_module('cope');
axxb = cope.axxbcovariance;

filename = fullfile('calib', 'pattern_tfs');
pattern_tfs =  py.pickle.load(py.open( filename, 'rb' ) );
filename = fullfile('calib', 'robot_tfs');
robot_tfs =  py.pickle.load(py.open( filename, 'rb' ) );

cam_pose = zeros(4, 4, size(pattern_tfs, 2));
robot_pose = zeros(4, 4, size(pattern_tfs, 2));
for i = 1 : size(pattern_tfs, 2)
    cam_pose(:, :, i) = double_(py.numpy.double(pattern_tfs(i)), 4, 4);
    robot_pose(:, :, i) = double_(py.numpy.double(robot_tfs(i)), 4, 4);
end


sigmaRa = 1e-10 * eye(3);
sigmata = 1e-10 * eye(3);

sigmaRb = [  4.15625435e-05,  -2.88693145e-05,  -6.06526440e-06;
            -2.88693145e-05,   3.20952008e-04,  -1.44817304e-06;
            -6.06526440e-06,  -1.44817304e-06,   1.43937081e-05];
        
sigmatb = [  1.95293655e-04,   2.12627214e-05,  -1.06674886e-05;
             2.12627214e-05,   4.44314426e-05,   3.86787591e-06;
            -1.06674886e-05,   3.86787591e-06,   2.13069579e-05];
        
datasize = py.len(pattern_tfs);
ksamples = 3;
iters = 100;
Rxlist = py.list;
sigmaRx_list = py.list;
txlist = py.list;
sigmatx_list = py.list;

t_funcs = {@t1_hand_eye_func_new, @t2_hand_eye_func_new, @t3_hand_eye_func_new};
j = 1;
n = 0;
while(true)
    AAA = zeros(4, 4, ksamples);
    BBB = zeros(4, 4, ksamples);
    alpha = py.list;
    beta = py.list;
    ta = py.list;
    tb = py.list;
    for i = 1 : ksamples
        rand_number_1 = 0;
        rand_number_2 = 0;
        
        while(rand_number_1 == 0 || rand_number_2 == 0)
            rand_number_1 = int32(py.numpy.random.uniform(0,datasize));
            rand_number_2 = int32(py.numpy.random.uniform(0,datasize));
            while rand_number_1 == rand_number_2
                rand_number_2 = int32(py.numpy.random.uniform(0,datasize));
            end
        end
        
        A = py.numpy.dot(robot_tfs(rand_number_1), py.numpy.linalg.inv(robot_tfs(rand_number_2)));
        B = py.numpy.dot(pattern_tfs(rand_number_1), py.numpy.linalg.inv(pattern_tfs(rand_number_2)));
        
        
        A = double(py.array.array('d', py.numpy.nditer(A)));
        A = reshape(A, [4 4])';
        
        B = double(py.array.array('d', py.numpy.nditer(B)));
        B = reshape(B, [4 4])';
        
        alpha.append(py.SE3lib.RotToVec(mat2nparray(A(1 : 3, 1 : 3))))
        beta.append(py.SE3lib.RotToVec(mat2nparray(B(1 : 3, 1 : 3))))
        ta.append(mat2nparray(A(1 : 3, 4)))
        tb.append(mat2nparray(B(1 : 3, 4)))

        AAA(:, :, i) = A;
        BBB(:, :, i) = B;
    end
    
    res = axxb.FCParkSolution(alpha, beta, ta, tb);
    Rxinit = pymat(res(1), 3, 3);
    txinit = pymat(res(2), 3, 1);
    rot_res = axxb.IterativeSolutionRot(beta, alpha, ...
                                      mat2nparray(sigmaRa), mat2nparray(sigmaRb), ...
                                      mat2nparray(Rxinit));
    Rxhat = py.numpy.double(rot_res{1});
    sigmaRx = py.numpy.double(rot_res{2});
    rot_converged = rot_res{3};
    betahat = py.numpy.double(rot_res{4});
    alphahat = py.numpy.double(rot_res{5});
    if(strcmp(char(rot_res{6}), 'None'))
        continue;
    end
    
    sigmaRbeta = py.numpy.double(rot_res{6});
    sigmabeta = py.numpy.double(rot_res{7});
    sigmaRahat = py.numpy.double(rot_res{8});
    sigmaRRa = py.numpy.double(rot_res{9});
    res = axxb.IterativeSolutionTrans(betahat, alphahat, ta, tb, ...
                                  Rxhat, sigmaRahat, ...
                                  mat2nparray(sigmaRb), mat2nparray(sigmata), ...
                                  mat2nparray(sigmatb), sigmaRx, sigmaRbeta, ...
                                  mat2nparray(txinit), 10);
    txhat = res{1};
    sigmatx = res{2};
    trans_converged = res{3};
    
    if(rot_converged && trans_converged)
        Rxlist.append(Rxhat)
        sigmaRx_list.append(sigmaRx)
        txlist.append(txhat.reshape(uint8(3)))
        sigmatx_list.append(sigmatx)
    else
        continue;
    end
    
    n = n + 1;
    if(n > iters)
        break;
    end
    
        [W, Q] = hand_eye_WQ_new(A, B);
        Ws(:, :, j) = W;
        Qs(:, :, j) = Q;
        [D, G, c, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG] = hand_eye_DGc_new(AAA, BBB);
        R = andreff(AAA, BBB);
        q = dcm2quat(real(R)).';
        if(q(1) < 0)
            q = - q;
        end
        [~, ~, X_] = QPEP_lm_single(q, 100, 5e-3, ...
                          @eq_hand_eye_func_new, @Jacob_hand_eye_func_new, t_funcs, ...
                          coef_f_q_sym, coefs_tq, pinvG);
        q = dcm2quat(X_(1 : 3, 1 : 3)).';
        t = X_(1 : 3, 4);
        if(q(1) < 0)
            q = - q;
        end
        
        Rs(:, :, j) = X_(1 : 3, 1 : 3);
        qs(:, j) = q;
        ts(:, j) = t;
        y = y_func_hand_eye_new(q);
        v = v_func_hand_eye_new(q);
        ys(:, j) = y;
        vs(:, j) = v;
        Ds(:, :, j) = D;
        ds(:, j) = vec(D);
        Gs(:, :, j) = G;
        gs(:, j) = vec(G);
        cs(:, j) = c;
        Dvs(:, j) = D * v;
        Gys(:, j) = G * y;
        ps(:, j) = D * v + G * y + c;
        DvGys(:, j) = D * v + G * y;
        j = j + 1;
    
 
end

logRx_list = py.list;
for i = 1 : size(Rxlist, 2)
    Rx = Rxlist(i);
    vv = py.SE3lib.RotToVec(Rx{1});
    if(isfinite(double_(vv, 3, 1)))
        logRx_list.append(vv);
    end
end
avg_log = py.numpy.average(logRx_list, uint8(0))
avg_Rx = py.SE3lib.VecToRot(avg_log)
inv_avg_Rx = py.numpy.linalg.inv(py.numpy.double(avg_Rx))
avg_tx = py.numpy.average(txlist, uint8(0))


xiRx_list = py.list;
RRs = zeros(3, 3, size(Rxlist, 2));
qRs = zeros(4, size(Rxlist, 2));
for i = 1 : size(Rxlist, 2) - 1
    Rx = Rxlist(i);
    RRs(:, :, i) = double_(py.numpy.double(Rx), 3, 3);
    qRs(:, i) = dcm2quat(Rs(:, :, i));
    if(qRs(1, i) < 0)
        qRs(:, i) = - qRs(:, i);
    end
    xiRx_list.append(py.SE3lib.RotToVec(py.numpy.dot(Rx{1}, inv_avg_Rx)))
end
real_sigmaRx_using_avg = cov(double_(py.numpy.double(xiRx_list), size(Rxlist, 2) - 1, 3))
avg_est_sigmaRx = double_(py.numpy.average(sigmaRx_list, uint8(0)), 3, 3)

figure(1);
plot_R_cov(real_sigmaRx_using_avg, avg_est_sigmaRx, RRs)


num = j - 1;
Sigma_q_stat = covx(qs.', qs.');
Sigma_y_stat = covx(ys.', ys.');
Sigma_v_stat = covx(vs.', vs.');

Sigma_d_stat = covx(ds.', ds.');
Sigma_g_stat = covx(gs.', gs.');
Sigma_c_stat = covx(cs.', cs.');

Sigma_d_q_stat = covx(ds.', qs.');
Sigma_q_d_stat = covx(qs.', ds.');
Sigma_g_q_stat = covx(gs.', qs.');
Sigma_q_g_stat = covx(qs.', gs.');
Sigma_c_q_stat = covx(cs.', qs.');
Sigma_q_c_stat = covx(qs.', cs.');

Sigma_d_g_stat = covx(ds.', gs.');
Sigma_g_d_stat = covx(gs.', ds.');
Sigma_d_c_stat = covx(ds.', cs.');
Sigma_c_d_stat = covx(cs.', ds.');
Sigma_c_g_stat = covx(cs.', gs.');
Sigma_g_c_stat = covx(gs.', cs.');

Sigma_p_stat = covx(ps.', ps.');

Sigma_Dv_stat = covx(Dvs.', Dvs.');
Sigma_Gy_stat = covx(Gys.', Gys.');
Sigma_DvGy_stat = covx(DvGys.', DvGys.');

mean_D = zeros(3, 28);
mean_d = mean(ds, 2);
mean_G = zeros(3, 9);
mean_g = mean(gs, 2);
mean_Dv = mean(Dvs, 2);
mean_Gy = mean(Gys, 2);
mean_W = zeros(4, 64);
mean_Q = zeros(4, 4);
mean_xi = zeros(3, 1);

for i = 1 : num
    mean_D = mean_D + 1 / num * Ds(:, :, i);
    mean_G = mean_G + 1 / num * Gs(:, :, i);
    mean_W = mean_W + 1 / num * Ws(:, :, i);
    mean_Q = mean_Q + 1 / num * Qs(:, :, i);
    mean_xi = mean_xi + 1 / num * vex(logR(Rs(:, :, i)));
end
mean_R = expm(skew(mean_xi));

Sigma_R_stat = 0;
for i = 1 : num
    res = Rs(:, :, i) - mean_R;
    Sigma_R_stat = Sigma_R_stat + 1 / num * res * res.';
end


q = qs(:, num);
q = q ./ norm(q);
P1 = P1_matrix(q);
P2 = P2_matrix(q);
P3 = P3_matrix(q);
D = mean_D;
G = mean_G;
W = mean_W;
Q = mean_Q;
y = y_func_hand_eye_new(q);
v = v_func_hand_eye_new(q);
V_cal = kron(v.', eye(3));
Y_cal = kron(y.', eye(3));

syms q0_ q1_ q2_ q3_
q_ = [q0_; q1_; q2_; q3_];
y_ = y_func_hand_eye_new(q_);
v_ = v_func_hand_eye_new(q_);

partial_v_q_sym = jacobian(v_, q_);
partial_v_q_sym_func = matlabFunction(partial_v_q_sym, 'Vars', {q_});
partial_v_q_sym_val = partial_v_q_sym_func(q);

partial_y_q_sym = jacobian(y_, q_);
partial_y_q_sym_func = matlabFunction(partial_y_q_sym, 'Vars', {q_});
partial_y_q_sym_val = partial_y_q_sym_func(q);

global F cov_left
F = D * partial_v_q_sym_val + G * partial_y_q_sym_val;

cov_left = V_cal * Sigma_d_stat * V_cal.' + Y_cal * Sigma_g_stat * Y_cal.' + Sigma_c_stat + ...
           V_cal * Sigma_d_g_stat * Y_cal.' + V_cal * Sigma_d_c_stat + Y_cal * Sigma_g_c_stat + ...
           Sigma_c_d_stat * V_cal.' + Sigma_c_g_stat * Y_cal.' + Y_cal * Sigma_g_d_stat * V_cal.';

          
scalings = 1;
solver = 'sedumi';
verbose = false;
for i = 1 : length(scalings)
    num_ = ksamples;
    scaling = scalings(i);
    AA = zeros(num_ * 9, 16);
    bb = zeros(num_ * 9, 1);
    for j = 1 : num_
        [D_, G_, c_] = hand_eye_DGc_new(AAA(:, :, j), BBB(:, :, j));
        F_ = D_ * partial_v_q_sym_val + G_ * partial_y_q_sym_val;
        right2 = Sigma_c_stat - (V_cal * Sigma_d_g_stat * Y_cal.' + ...
              Y_cal * Sigma_g_d_stat * V_cal.') - ...
              (V_cal * Sigma_d_stat * V_cal.' + Y_cal * Sigma_g_stat * Y_cal.');
        bb(9 * (j - 1) + 1 : 9 * j) = vec(right2);
        AA(9 * (j - 1) + 1 : 9 * j, :) = 3 * kron(F_, F_);
    end
    
    [XX, ff] = optimize_quat_cov2(q, F, scaling * cov_left, solver, verbose);
    scale = abs(norm(cov_left, 'fro') / norm(F * XX * F.', 'fro'));
    
    figure;
    plot_q_cov(Sigma_q_stat, XX * scale, qs);
end
