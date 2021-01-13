function g = logR(R)
phi = acos((trace(R) - 1) / 2);
g = phi / (2 * sin(phi)) * (R - R.');
end