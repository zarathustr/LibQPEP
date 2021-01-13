function J = J_func_hand_eye(X, A, B)
J = 0;
len = size(A, 3);
for i = 1 : len
    res = A(:, :, i) * X - X * B(:, :, i);
    J = J + 1 / len * trace(res.' * res);
end
end