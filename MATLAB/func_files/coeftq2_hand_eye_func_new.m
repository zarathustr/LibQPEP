function obj = coeftq2_hand_eye_func_new(in1)
coef_J12 = in1(:,12);
coef_J22 = in1(:,22);
coef_J29 = in1(:,29);
coef_J34 = in1(:,34);
coef_J44 = in1(:,44);
coef_J51 = in1(:,51);
coef_J56 = in1(:,56);
coef_J63 = in1(:,63);
coef_J68 = in1(:,68);
coef_J73 = in1(:,73);
coef_J82 = in1(:,82);
obj = [coef_J12,coef_J22,coef_J29,coef_J34,coef_J44,coef_J51,coef_J56,coef_J63,coef_J68,coef_J73,coef_J82];
end
