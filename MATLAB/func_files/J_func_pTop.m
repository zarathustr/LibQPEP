function J = J_func_pTop(q, t, r, b, nv)
len = size(r, 1);
R = q2R(q);
J = 0;
for i = 1 : len
   rr = r(i, :).';
   bb = b(i, :).';
   J = J + 1 / len * (nv(i, :) * (R * rr + t - bb))^2; 
end
end