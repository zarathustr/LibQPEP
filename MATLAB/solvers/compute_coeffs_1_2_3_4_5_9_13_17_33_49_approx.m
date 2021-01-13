% 
% LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems (QPEPs),
%          It also gives highly accurate uncertainty description of the solutions.
%
%
% Article: 
%      Wu, J., Zheng, Y., Gao, Z., Jiang, Y., Hu, X., Zhu, Y., Jiao, J., Liu, M. (2020)
%           Quadratic Pose Estimation Problems: Unified Solutions, 
%           Solvability/Observability Analysis and Uncertainty Description 
%           in A Globally Optimal Framework.
%
%
% Authors:      Jin Wu and Ming Liu
% Affiliation:  Hong Kong University of Science and Technology (HKUST)
% Emails:       jin_wu_uestc@hotmail.com; eelium@ust.hk
% Websites:     https://zarathustr.github.io
%               https://ram-lab.com


function [coeffs] = compute_coeffs_1_2_3_4_5_9_13_17_33_49_approx(data)
if(isnumeric(data))
    coeffs = zeros(1, 62);
else
    coeffs = sym(zeros(1, 62));
end
coeffs(1) = 3 * data(1);
coeffs(2) = data(19);
coeffs(3) = 6 * data(4);
coeffs(4) = 3 * data(22);
coeffs(5) = 3 * data(10);
coeffs(6) = 3 * data(28);
coeffs(7) = data(37);
coeffs(8) = 6 * data(7);
coeffs(9) = 3 * data(25);
coeffs(10) = 6 * data(13);
coeffs(11) = 6 * data(31);
coeffs(12) = 3 * data(40);
coeffs(13) = 3 * data(16);
coeffs(14) = 3 * data(34);
coeffs(15) = 3 * data(43);
coeffs(16) = data(46);
coeffs(17) = -data(49);
coeffs(18) = -data(52);
coeffs(19) = -data(55);
coeffs(20) = -data(58);
coeffs(21) = 3 * data(2);
coeffs(22) = data(20);
coeffs(23) = 6 * data(5);
coeffs(24) = 3 * data(23);
coeffs(25) = 3 * data(11);
coeffs(26) = 3 * data(29);
coeffs(27) = data(38);
coeffs(28) = 6 * data(8);
coeffs(29) = 3 * data(26);
coeffs(30) = 6 * data(14);
coeffs(31) = 6 * data(32);
coeffs(32) = 3 * data(41);
coeffs(33) = 3 * data(17);
coeffs(34) = 3 * data(35);
coeffs(35) = 3 * data(44);
coeffs(36) = data(47);
coeffs(37) = -data(50);
coeffs(38) = -data(53);
coeffs(39) = -data(56);
coeffs(40) = -data(59);
coeffs(41) = 3 * data(3);
coeffs(42) = data(21);
coeffs(43) = 6 * data(6);
coeffs(44) = 3 * data(24);
coeffs(45) = 3 * data(12);
coeffs(46) = 3 * data(30);
coeffs(47) = data(39);
coeffs(48) = 6 * data(9);
coeffs(49) = 3 * data(27);
coeffs(50) = 6 * data(15);
coeffs(51) = 6 * data(33);
coeffs(52) = 3 * data(42);
coeffs(53) = 3 * data(18);
coeffs(54) = 3 * data(36);
coeffs(55) = 3 * data(45);
coeffs(56) = data(48);
coeffs(57) = -data(51);
coeffs(58) = -data(54);
coeffs(59) = -data(57);
coeffs(60) = -data(60);
coeffs(61) = 1;
coeffs(62) = -1;
end

