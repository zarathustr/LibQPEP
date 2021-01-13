function [R, t, X] = andreff(AA,BB)
% Solves the problem AX=XB
% using the formulation of
%
% On-line Hand-Eye Calibration.
% N. Andreff, R. Horaud, B. Espiau 
%
% Mili Shah
% July 2014

[~,n] = size(AA); n = n/4;

A = zeros(12*n,12);
b = zeros(12*n,1);
for i = 1:n
    Ra = AA(1:3,4*i-3:4*i-1);
    Rb = BB(1:3,4*i-3:4*i-1);
    ta = AA(1:3,4*i);
    tb = BB(1:3,4*i);
    A(12*i-11:12*i-3,1:9) = eye(9) - kron(Rb,Ra);
    A(12*i-2:12*i,:) = [kron(tb',eye(3)) eye(3)-Ra];
    b(12*i-2:12*i) = ta;
end
x = pinv(A) * b;

X = reshape(x(1:9),3,3)';
X = sign(det(X))/abs(det(X))^(1/3)*X;

[u, ~, v] = svds(X, 3); X = u*v'; if det(X)<0, X = u*diag([1 1 -1])*v'; end
X = [X' x(10:12);[0 0 0 1]];
R = X(1 : 3, 1 : 3);
t = X(1 : 3, 4);

end
