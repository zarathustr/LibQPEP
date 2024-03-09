clear all
close all
clc

syms q [4, 1]
syms dq [4, 1]

u = kron(q, kron(q, q));
qq = q + dq;
uu = kron(qq, kron(qq, qq));

du = expand(uu - u);

J_u_q = jacobian(u, q);
J_u_q_func = matlabFunction(J_u_q, 'Vars', {q});
expand(du - kron(dq, kron(dq, dq)) - J_u_q_func(dq) * q - J_u_q_func(q) * dq)