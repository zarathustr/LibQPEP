function x = dcm2Lie(Rs)
len = size(Rs, 3);
x = zeros(len, 3);
for i = 1 : len
    x(i, :) = vex(logR(Rs(:, :, i)));
end
end