function obj = v_func_hand_eye_new(in1)
q0 = in1(1,:);
q1 = in1(2,:);
q2 = in1(3,:);
q3 = in1(4,:);
t2 = q0.^3;
t3 = q1.^2;
t4 = q1.^3;
t5 = q2.^2;
t6 = q2.^3;
t7 = q3.^2;
t8 = q3.^3;
obj = [q1.*t2;q2.*t2;q3.*t2;q0.*t4;q0.*q2.*t3;q0.*q3.*t3;q0.*q1.*t5;q0.*q1.*q2.*q3;q0.*q1.*t7;q0.*t6;q0.*q3.*t5;q0.*q2.*t7;q0.*t8;t3.^2;q2.*t4;q3.*t4;t3.*t5;q2.*q3.*t3;t3.*t7;q1.*t6;q1.*q3.*t5;q1.*q2.*t7;q1.*t8;t5.^2;q3.*t6;t5.*t7;q2.*t8;t7.^2];
end
