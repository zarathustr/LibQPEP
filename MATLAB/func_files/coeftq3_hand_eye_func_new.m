function obj = coeftq3_hand_eye_func_new(in1)
coef_J13 = in1(:,13);
coef_J23 = in1(:,23);
coef_J30 = in1(:,30);
coef_J35 = in1(:,35);
coef_J45 = in1(:,45);
coef_J52 = in1(:,52);
coef_J57 = in1(:,57);
coef_J64 = in1(:,64);
coef_J69 = in1(:,69);
coef_J74 = in1(:,74);
coef_J84 = in1(:,84);
obj = [coef_J13,coef_J23,coef_J30,coef_J35,coef_J45,coef_J52,coef_J57,coef_J64,coef_J69,coef_J74,coef_J84];
end
