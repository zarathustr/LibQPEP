function obj = coeftq1_hand_eye_func_new(in1)
coef_J11 = in1(:,11);
coef_J21 = in1(:,21);
coef_J28 = in1(:,28);
coef_J33 = in1(:,33);
coef_J43 = in1(:,43);
coef_J50 = in1(:,50);
coef_J55 = in1(:,55);
coef_J62 = in1(:,62);
coef_J67 = in1(:,67);
coef_J72 = in1(:,72);
coef_J79 = in1(:,79);
obj = [coef_J11,coef_J21,coef_J28,coef_J33,coef_J43,coef_J50,coef_J55,coef_J62,coef_J67,coef_J72,coef_J79];
end
