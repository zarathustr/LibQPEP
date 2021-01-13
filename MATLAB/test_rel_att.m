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
%
% test_rel_att.m: The QPEP illustration of range-based
%                 relative attitude estimtion



clear all
close all
clc

if(verLessThan('matlab', '8.0.0'))
   error('The MATLAB version is too old to be supported.'); 
end

addpath('func_files');
addpath('solvers');
addpath('utils');
addpath('homotopy');

p00 = randn(3, 1);
d00 = sqrt(p00.' * p00);
C0 = orthonormalize(randn(3, 3));

len = 10;
p1 = zeros(3, len);
p2 = zeros(3, len);
epsilon = zeros(len, 1);
noise = 1e-2; % Noise level
for i = 1 : len
    p1(:, i) = randn(3, 1);
    p2(:, i) = randn(3, 1);
    epsilon(i) = (p1(:, i) - p00).' * C0 * p2(:, i) + p1(:, i).' * p00 + noise * randn(1, 1);
end

syms q0 q1 q2 q3
q = [q0; q1; q2; q3];
syms t1 t2 t3;
t = [t1; t2; t3];
R = q2R(q);
rr = R.' * t;
syms r1 r2 r3 
r = [r1; r2; r3];
eqs = sym(zeros(len + 2, 1));
for i = 1 : len
    eqs(i) = p1(:, i).' * R * p2(:, i) + p1(:, i).' * t - r.' * p2(:, i) - epsilon(i);
end
eqs(len + 1) = q.' * q - 1;
eqs(len + 2) = r.' * r - t.' * t;
eqs = expand(eqs);
x = [q; t; r];
H = expand(jacobian(eqs, x).' * eqs);
assumeAlso(q.' * q == 1);
assumeAlso(r.' * r == t.' * t);
eq = vpa(expand(simplify(H)), 32);
ss = vpasolve(eq(5 : 10), [t; r]);
eq_ = eq(1 : 4);
eq_ = subs(eq_, t1, ss.t1);
eq_ = subs(eq_, t2, ss.t2);
eq_ = subs(eq_, t3, ss.t3);
t_func = matlabFunction([ss.t1; ss.t2; ss.t3], 'Vars', {q});
eq_ = subs(eq_, r1, ss.r1);
eq_ = subs(eq_, r2, ss.r2);
eq_ = subs(eq_, r3, ss.r3);
eq_ = vpa(expand(eval(eq_)), 32);
syms lambda
eqs = [
    eq_;
    q.' * q - 1;
    ]

str = '';
for i = 1 : length(eqs)
    str = strcat(str, sprintf(' PP{%d} = char(vpa(%%s, 32));', i));
end
    
str_ = sprintf(str, char(eqs(1)), ...
                    char(eqs(2)), ...
                    char(eqs(3)), ...
                    char(eqs(4)), ...
                    char(eqs(5)));
eval(str_);
[S, vars] = psolve(PP);
S = S.';
SS = S;
for i = 1 : length(vars)
    if(strcmp(vars{i}, 'q0'))
        SS(:, 1) = S(:, i);
    elseif(strcmp(vars{i}, 'q1'))
        SS(:, 2) = S(:, i);
    elseif(strcmp(vars{i}, 'q2'))
        SS(:, 3) = S(:, i);
    elseif(strcmp(vars{i}, 'q3'))
        SS(:, 4) = S(:, i);
    elseif(strcmp(vars{i}, 'lambda'))
        SS(:, 5) = S(:, i);
    end
end
S = real(SS);
xs_ = S;
sols = SS.';



num = size(sols, 2);
sol = zeros(4, num);
ts = zeros(3, num);
Ls = zeros(len, 1);
for i = 1 : num
    sol(:, i) = real(sols(1 : 4, i));
    sol(:, i) = sol(:, i) ./ norm(sol(:, i));
    C = q2R(sol(:, i));
    t = t_func(sol(:, i));
    ts(:, i) = t;
    loss = 0;
    for j = 1 : len
        loss = loss + (epsilon(j) - (p1(:, j) - t).' * C * p2(:, j) - p1(:, j).' * t)^2;
    end
    Ls(i) = loss;
end
[~, idx] = sort(Ls);

q_ = positive_quat(sol(:, idx(1)).')
q_true = positive_quat(dcm2quat(C0))

t_ = ts(:, idx(1)).'
t_true = p00.'








