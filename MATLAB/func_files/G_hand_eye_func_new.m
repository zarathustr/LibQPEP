function obj = G_hand_eye_func_new(in1)
coef_J76 = in1(:,76);
coef_J77 = in1(:,77);
coef_J78 = in1(:,78);
coef_J80 = in1(:,80);
coef_J81 = in1(:,81);
coef_J83 = in1(:,83);
t2 = -coef_J77;
t3 = -coef_J78;
t4 = -coef_J81;
obj = reshape([coef_J76.*-2.0,t2,t3,t2,coef_J80.*-2.0,t4,t3,t4,coef_J83.*-2.0],[3,3]);
end
