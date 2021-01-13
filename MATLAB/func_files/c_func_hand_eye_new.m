function obj = c_func_hand_eye_new(in1,in2,in3,in4)
coef_f1_q_sym1 = in2(:,1);
coef_f2_q_sym1 = in3(:,1);
coef_f3_q_sym1 = in4(:,1);
coef_f1_q_sym11 = in2(:,11);
coef_f2_q_sym11 = in3(:,11);
coef_f3_q_sym11 = in4(:,11);
obj = [-coef_f1_q_sym1-coef_f1_q_sym11;-coef_f2_q_sym1-coef_f2_q_sym11;-coef_f3_q_sym1-coef_f3_q_sym11];
