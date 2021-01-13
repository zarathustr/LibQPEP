function obj = y_func_hand_eye_new(in1)
q0 = in1(1,:);
q1 = in1(2,:);
q2 = in1(3,:);
q3 = in1(4,:);
obj = [q0.*q1;q0.*q2;q0.*q3;q1.^2;q1.*q2;q1.*q3;q2.^2;q2.*q3;q3.^2];
