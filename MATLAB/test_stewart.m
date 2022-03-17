% 
% LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems (QPEPs),
%          It also gives highly accurate uncertainty description of the solutions.
%
%
% Article: 
%      Wu, J., Zheng, Y., Gao, Z., Jiang, Y., Hu, X., Zhu, Y., Jiao, J., Liu, M. (2020)
%           Quadratic Pose Estimation Problems: Globally Optimal Solutions, 
%           Solvability/Observability Analysis and Uncertainty Description.
%           IEEE Transactions on Robotics.
%           https://doi.org/10.1109/TRO.2022.3155880
%
%
% Authors:      Jin Wu and Ming Liu
% Affiliation:  Hong Kong University of Science and Technology (HKUST)
% Emails:       jin_wu_uestc@hotmail.com; eelium@ust.hk
% Websites:     https://zarathustr.github.io
%               https://ram-lab.com
%
%
% test_stewart.m: The QPEP illustration of forwart kinematics of hexapod
%                 Stewart platform


clear all
close all
clc

if(verLessThan('matlab', '8.0.0'))
   error('The MATLAB version is too old to be supported.'); 
end

format long g

addpath('func_files');
addpath('solvers');
addpath('utils');
addpath('homotopy');

R0 = angle2dcm(-8 * pi / 180, 12 * pi / 180, -15 * pi / 180, 'XYZ');
q0 = dcm2quat(R0).';
if(q0(1) < 0)
    q0 = - q0;
end
t0 = 1e-2 * randn(3, 1);
X0 = inv([R0, t0;
     zeros(1, 3), 1]);
 
base = [
     0.1448888739433600, 1,  0.0388228567653781;
    -0.0388228567653781, 1,  0.1448888739433600;
    -0.1060660171779820, 1,  0.1060660171779820;
    -0.1060660171779820, 1, -0.1060660171779820;
    -0.0388228567653781, 1, -0.1448888739433600;
     0.1448888739433600, 1, -0.0388228567653781;
    ].';
plat = [
     0.0707106781186548, 1, 0.07071067811865480;
     0.0258819045102521, 1, 0.09659258262890680;
    -0.0965925826289068, 1, 0.02588190451025210;
    -0.0965925826289068, 1, -0.0258819045102521;
     0.0258819045102521, 1, -0.0965925826289068;
     0.0707106781186548, 1, -0.0707106781186548;
    ].';



conv = [
    1, 0, 0;
    0, 0, 1;
    0, -1, 0;
    ];
height = 0.15;
base = conv * base;
plat = conv * plat;
plat(3, :) = plat(3, :) + height;
plat00 = plat;
base00 = base;



leg0 = zeros(6, 1);
plat0 = zeros(3, 6);
for i = 1 : 6
    plat0(:, i) = R0 * plat(:, i)  + t0;
    res = base(:, i) - plat0(:, i);
    leg0(i) = sqrt(res.' * res);
end
colors = linspecer(8);

figure(1);
subplot(1, 2, 1);
plot3(base(1, :), base(2, :), base(3, :), 'LineStyle', 'None', 'Marker', '.', 'MarkerSize', 10); hold on
plot3(base(1, 1 : 6), base(2, 1 : 6), base(3, 1 : 6), 'LineStyle', '-', 'LineWidth', 2, 'Marker', 'None'); hold on
plot3([base(1, 1), base(1, 6)], [base(2, 1), base(2, 6)], [base(3, 1), base(3, 6)], 'LineStyle', '-', 'LineWidth', 2, 'Marker', 'None'); hold on
plot3(plat0(1, :), plat0(2, :), plat0(3, :), 'LineStyle', 'None', 'Marker', '.', 'MarkerSize', 10); hold on
plot3(plat0(1, 1 : 6), plat0(2, 1 : 6), plat0(3, 1 : 6), 'LineStyle', '-', 'LineWidth', 2, 'Marker', 'None'); hold on
plot3([plat0(1, 1), plat0(1, 6)], [plat0(2, 1), plat0(2, 6)], [plat0(3, 1), plat0(3, 6)], 'LineStyle', '-', 'LineWidth', 2, 'Marker', 'None'); hold on
for i = 1 : 6
    plot3([base(1, i), plat0(1, i)], [base(2, i), plat0(2, i)], [base(3, i), plat0(3, i)], 'LineStyle', '-', 'LineWidth', 4, 'Marker', 'None'); hold on
end
hold on
fill3(base(1, :), base(2, :), base(3, :), colors(8, :)); hold on
fill3(plat0(1, :), plat0(2, :), plat0(3, :), colors(5, :)); hold off
grid on
grid minor
title('Reference Result', 'Interpreter', 'LaTeX', 'FontSize', 14);


base_ = base.';
plat_ = plat.';
counter = 1;
base = [];
plat = [];
for i = 1 : 6
    base = [base; base_(i, :)];
    plat = [plat; plat_(i, :)];
    if(i < 6)
        tmp = (base_(i, :) + base_(i + 1, :)) / 2;
        base = [base; tmp];
        
        tmp = (plat_(i, :) + plat_(i + 1, :)) / 2;
        plat = [plat; tmp];
    end
        
end
base = base.';
plat = plat.';
len = size(base, 2);


leg0 = zeros(len, 1);
plat0 = zeros(3, len);
for i = 1 : len
    plat0(:, i) = R0 * plat(:, i)  + t0;
    res = base(:, i) - plat0(:, i);
    leg0(i) = sqrt(res.' * res);
end


syms q0 q1 q2 q3
q = [q0; q1; q2; q3];
syms t1 t2 t3;
t = [t1; t2; t3];
R = q2R(q);
syms r1 r2 r3 r4
r = [r1; r2; r3; r4];
rr = r(1 : 3);
eqs = sym(zeros(len + 3, 1));
for i = 1 : len
    eqs(i) = base(:, i).' * base(:, i) - 2 * plat(:, i).' * R.' * base(:, i) + ...
             plat(:, i).' * plat(:, i) - 2 * base(:, i).' * t + 2 * plat(:, i).' * rr + r4 - leg0(i)^2;
end
eqs(len + 1) = q.' * q - 1;
eqs(len + 2) = rr.' * rr - r4;
eqs(len + 3) = t.' * t - r4;
eqs = expand(eqs);
x = [q; t; r];
H = expand(jacobian(eqs, x).' * eqs);
assumeAlso(q.' * q == 1);
assumeAlso(rr.' * rr - r4 == 0);
assumeAlso(t.' * t - r4 == 0);
eq = vpa(expand(simplify(H)), 32);

eq_ = eq(1 : 4);
G = jacobian(eq(5 : 11), [t; r]);
ss = - pinv(G) * (eq(5 : 11) - G * [t; r]);
ss = neglect_tiny_terms(ss, 32);
ss = ss.';
eq_ = eq(1 : 4);
eq_ = subs(eq_, t1, ss(1));
eq_ = subs(eq_, t2, ss(2));
eq_ = subs(eq_, t3, ss(3));
eq_ = subs(eq_, r1, ss(4));
eq_ = subs(eq_, r2, ss(5));
eq_ = subs(eq_, r3, ss(6));
eq_ = subs(eq_, r4, ss(7));
t_func = matlabFunction([ss(1); ss(2); ss(3)], 'Vars', {q});
r_func = matlabFunction([ss(4); ss(5); ss(6)], 'Vars', {q});
r4_func = matlabFunction(ss(7), 'Vars', {q});
eq_ = vpa(expand(eval(eq_)), 32);
syms lambda
eq_ = [
    neglect_tiny_terms(eq_, 32).';
    q.' * q - 1;
    ]

eqs = [
    eq_(1 : 4) + lambda * q;
    eq_(5);
    ];
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
Ls = 1e50 * ones(num, 1);
q_true = dcm2quat(R0).';
if(q_true(1) < 0)
    q_true = - q_true;
end
for i = 1 : num
    sol(:, i) = real(sols(1 : 4, i));
    sol(:, i) = sol(:, i) ./ norm(sol(:, i));
    if(sol(1, i) < 0)
        sol(:, i) = - sol(:, i);
    end
    C = q2R(sol(:, i));
    t = t_func(sol(:, i));
    r4 = r4_func(sol(:, i));
    ts(:, i) = t;
    res = abs(q_true - sol(:, i));
    loss = res.' * res;
    Ls(i) = loss;
end
[~, idx] = sort(Ls);

q_ = sol(:, idx(1)).'
q_true_ = q_true.'


R_ = quat2dcm(sol(:, idx(1)).');
t_ = t0;
plat0 = zeros(3, 6);
for i = 1 : 6
    plat0(:, i) = R_ * plat00(:, i)  + t_;
end
base = base00;



subplot(1, 2, 2);
plot3(base(1, :), base(2, :), base(3, :), 'LineStyle', 'None', 'Marker', '.', 'MarkerSize', 10); hold on
plot3(base(1, 1 : 6), base(2, 1 : 6), base(3, 1 : 6), 'LineStyle', '-', 'LineWidth', 2, 'Marker', 'None'); hold on
plot3([base(1, 1), base(1, 6)], [base(2, 1), base(2, 6)], [base(3, 1), base(3, 6)], 'LineStyle', '-', 'LineWidth', 2, 'Marker', 'None'); hold on
plot3(plat0(1, :), plat0(2, :), plat0(3, :), 'LineStyle', 'None', 'Marker', '.', 'MarkerSize', 10); hold on
plot3(plat0(1, 1 : 6), plat0(2, 1 : 6), plat0(3, 1 : 6), 'LineStyle', '-', 'LineWidth', 2, 'Marker', 'None'); hold on
plot3([plat0(1, 1), plat0(1, 6)], [plat0(2, 1), plat0(2, 6)], [plat0(3, 1), plat0(3, 6)], 'LineStyle', '-', 'LineWidth', 2, 'Marker', 'None'); hold on
for i = 1 : 6
    plot3([base(1, i), plat0(1, i)], [base(2, i), plat0(2, i)], [base(3, i), plat0(3, i)], 'LineStyle', '-', 'LineWidth', 4, 'Marker', 'None'); hold on
end
hold on
fill3(base(1, :), base(2, :), base(3, :), colors(8, :)); hold on
fill3(plat0(1, :), plat0(2, :), plat0(3, :), colors(5, :)); hold off
grid on
grid minor
title('QPEP Result', 'Interpreter', 'LaTeX', 'FontSize', 14);

if(~ispc())
    set(gcf, 'Position', [634 780 1159 320])
end








