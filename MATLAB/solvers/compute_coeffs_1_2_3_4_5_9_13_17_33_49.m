% 
% LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems (QPEPs),  
%          It also gives highly accurate uncertainty description of the solutions.
%
%
% Article: 
%      Wu,   J.,   Zheng,   Y.,   Gao,   Z.,   Jiang,   Y.,   Hu,   X.,   Zhu,   Y.,   Jiao,   J.,   Liu,   M. (2020)
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


function [coeffs] = compute_coeffs_1_2_3_4_5_9_13_17_33_49(data)
coeffs = zeros(1,  82);
coeffs(1) = 3 * data(1);
coeffs(2) = data(25);
coeffs(3) = 6 * data(5);
coeffs(4) = 3 * data(29);
coeffs(5) = 3 * data(13);
coeffs(6) = 3 * data(37);
coeffs(7) = data(49);
coeffs(8) = 6 * data(9);
coeffs(9) = 3 * data(33);
coeffs(10) = 6 * data(17);
coeffs(11) = 6 * data(41);
coeffs(12) = 3 * data(53);
coeffs(13) = 3 * data(21);
coeffs(14) = 3 * data(45);
coeffs(15) = 3 * data(57);
coeffs(16) = data(61);
coeffs(17) = 1;
coeffs(18) = -data(65);
coeffs(19) = -data(69);
coeffs(20) = -data(73);
coeffs(21) = -data(77);
coeffs(22) = 3 * data(2);
coeffs(23) = data(26);
coeffs(24) = 6 * data(6);
coeffs(25) = 3 * data(30);
coeffs(26) = 3 * data(14);
coeffs(27) = 3 * data(38);
coeffs(28) = data(50);
coeffs(29) = 6 * data(10);
coeffs(30) = 3 * data(34);
coeffs(31) = 6 * data(18);
coeffs(32) = 6 * data(42);
coeffs(33) = 3 * data(54);
coeffs(34) = 3 * data(22);
coeffs(35) = 3 * data(46);
coeffs(36) = 3 * data(58);
coeffs(37) = data(62);
coeffs(38) = -data(66);
coeffs(39) = -data(70);
coeffs(40) = -data(74);
coeffs(41) = -data(78);
coeffs(42) = 3 * data(3);
coeffs(43) = data(27);
coeffs(44) = 6 * data(7);
coeffs(45) = 3 * data(31);
coeffs(46) = 3 * data(15);
coeffs(47) = 3 * data(39);
coeffs(48) = data(51);
coeffs(49) = 6 * data(11);
coeffs(50) = 3 * data(35);
coeffs(51) = 6 * data(19);
coeffs(52) = 6 * data(43);
coeffs(53) = 3 * data(55);
coeffs(54) = 3 * data(23);
coeffs(55) = 3 * data(47);
coeffs(56) = 3 * data(59);
coeffs(57) = data(63);
coeffs(58) = -data(67);
coeffs(59) = -data(71);
coeffs(60) = -data(75);
coeffs(61) = -data(79);
coeffs(62) = 3 * data(4);
coeffs(63) = data(28);
coeffs(64) = 6 * data(8);
coeffs(65) = 3 * data(32);
coeffs(66) = 3 * data(16);
coeffs(67) = 3 * data(40);
coeffs(68) = data(52);
coeffs(69) = 6 * data(12);
coeffs(70) = 3 * data(36);
coeffs(71) = 6 * data(20);
coeffs(72) = 6 * data(44);
coeffs(73) = 3 * data(56);
coeffs(74) = 3 * data(24);
coeffs(75) = 3 * data(48);
coeffs(76) = 3 * data(60);
coeffs(77) = data(64);
coeffs(78) = -data(68);
coeffs(79) = -data(72);
coeffs(80) = -data(76);
coeffs(81) = -data(80);
coeffs(82) = -1;
end