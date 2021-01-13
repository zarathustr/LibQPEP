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


function [coeffs] = compute_coeffs_solver_WQ_approx(data)
if(isnumeric(data))
    coeffs = zeros(1, 74);
else
    coeffs = sym(zeros(1, 74));
end
coeffs(1) = data(1);
coeffs(2) = data(4) + data(13) + data(49);
coeffs(3) = data(16) + data(52) + data(61);
coeffs(4) = data(64);
coeffs(5) = data(7) + data(25) + data(97);
coeffs(6) = data(19) + data(28) + data(55) + data(73) + data(100) + data(109);
coeffs(7) = data(67) + data(76) + data(112);
coeffs(8) = data(31) + data(103) + data(121);
coeffs(9) = data(79) + data(115) + data(124);
coeffs(10) = data(127);
coeffs(11) = data(10) + data(37) + data(145);
coeffs(12) = data(22) + data(40) + data(58) + data(85) + data(148) + data(157);
coeffs(13) = data(70) + data(88) + data(160);
coeffs(14) = data(34) + data(43) + data(106) + data(133) + data(151) + data(169);
coeffs(15) = data(82) + data(91) + data(118) + data(136) + data(163) + data(172);
coeffs(16) = data(130) + data(139) + data(175);
coeffs(17) = data(46) + data(154) + data(181);
coeffs(18) = data(94) + data(166) + data(184);
coeffs(19) = data(142) + data(178) + data(187);
coeffs(20) = data(190);
coeffs(21) = -data(193);
coeffs(22) = -data(196);
coeffs(23) = -data(199);
coeffs(24) = -data(202);
coeffs(25) = data(2);
coeffs(26) = data(5) + data(14) + data(50);
coeffs(27) = data(17) + data(53) + data(62);
coeffs(28) = data(65);
coeffs(29) = data(8) + data(26) + data(98);
coeffs(30) = data(20) + data(29) + data(56) + data(74) + data(101) + data(110);
coeffs(31) = data(68) + data(77) + data(113);
coeffs(32) = data(32) + data(104) + data(122);
coeffs(33) = data(80) + data(116) + data(125);
coeffs(34) = data(128);
coeffs(35) = data(11) + data(38) + data(146);
coeffs(36) = data(23) + data(41) + data(59) + data(86) + data(149) + data(158);
coeffs(37) = data(71) + data(89) + data(161);
coeffs(38) = data(35) + data(44) + data(107) + data(134) + data(152) + data(170);
coeffs(39) = data(83) + data(92) + data(119) + data(137) + data(164) + data(173);
coeffs(40) = data(131) + data(140) + data(176);
coeffs(41) = data(47) + data(155) + data(182);
coeffs(42) = data(95) + data(167) + data(185);
coeffs(43) = data(143) + data(179) + data(188);
coeffs(44) = data(191);
coeffs(45) = -data(194);
coeffs(46) = -data(197);
coeffs(47) = -data(200);
coeffs(48) = -data(203);
coeffs(49) = data(3);
coeffs(50) = data(6) + data(15) + data(51);
coeffs(51) = data(18) + data(54) + data(63);
coeffs(52) = data(66);
coeffs(53) = data(9) + data(27) + data(99);
coeffs(54) = data(21) + data(30) + data(57) + data(75) + data(102) + data(111);
coeffs(55) = data(69) + data(78) + data(114);
coeffs(56) = data(33) + data(105) + data(123);
coeffs(57) = data(81) + data(117) + data(126);
coeffs(58) = data(129);
coeffs(59) = data(12) + data(39) + data(147);
coeffs(60) = data(24) + data(42) + data(60) + data(87) + data(150) + data(159);
coeffs(61) = data(72) + data(90) + data(162);
coeffs(62) = data(36) + data(45) + data(108) + data(135) + data(153) + data(171);
coeffs(63) = data(84) + data(93) + data(120) + data(138) + data(165) + data(174);
coeffs(64) = data(132) + data(141) + data(177);
coeffs(65) = data(48) + data(156) + data(183);
coeffs(66) = data(96) + data(168) + data(186);
coeffs(67) = data(144) + data(180) + data(189);
coeffs(68) = data(192);
coeffs(69) = -data(195);
coeffs(70) = -data(198);
coeffs(71) = -data(201);
coeffs(72) = -data(204);
coeffs(73) = 1;
coeffs(74) = -1;
