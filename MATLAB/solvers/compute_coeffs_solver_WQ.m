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

function [coeffs] = compute_coeffs_solver_WQ(data)
if(isnumeric(data))
    coeffs = zeros(1, 98);
else
    coeffs = sym(zeros(1, 98));
end
coeffs(1) = data(1);
coeffs(2) = data(5) + data(17) + data(65);
coeffs(3) = data(21) + data(69) + data(81);
coeffs(4) = data(85);
coeffs(5) = data(9) + data(33) + data(129);
coeffs(6) = data(25) + data(37) + data(73) + data(97) + data(133) + data(145);
coeffs(7) = data(89) + data(101) + data(149);
coeffs(8) = data(41) + data(137) + data(161);
coeffs(9) = data(105) + data(153) + data(165);
coeffs(10) = data(169);
coeffs(11) = data(13) + data(49) + data(193);
coeffs(12) = data(29) + data(53) + data(77) + data(113) + data(197) + data(209);
coeffs(13) = data(93) + data(117) + data(213);
coeffs(14) = data(45) + data(57) + data(141) + data(177) + data(201) + data(225);
coeffs(15) = data(109) + data(121) + data(157) + data(181) + data(217) + data(229);
coeffs(16) = data(173) + data(185) + data(233);
coeffs(17) = data(61) + data(205) + data(241);
coeffs(18) = data(125) + data(221) + data(245);
coeffs(19) = data(189) + data(237) + data(249);
coeffs(20) = data(253);
coeffs(21) = 1;
coeffs(22) = -data(257);
coeffs(23) = -data(261);
coeffs(24) = -data(265);
coeffs(25) = -data(269);
coeffs(26) = data(2);
coeffs(27) = data(6) + data(18) + data(66);
coeffs(28) = data(22) + data(70) + data(82);
coeffs(29) = data(86);
coeffs(30) = data(10) + data(34) + data(130);
coeffs(31) = data(26) + data(38) + data(74) + data(98) + data(134) + data(146);
coeffs(32) = data(90) + data(102) + data(150);
coeffs(33) = data(42) + data(138) + data(162);
coeffs(34) = data(106) + data(154) + data(166);
coeffs(35) = data(170);
coeffs(36) = data(14) + data(50) + data(194);
coeffs(37) = data(30) + data(54) + data(78) + data(114) + data(198) + data(210);
coeffs(38) = data(94) + data(118) + data(214);
coeffs(39) = data(46) + data(58) + data(142) + data(178) + data(202) + data(226);
coeffs(40) = data(110) + data(122) + data(158) + data(182) + data(218) + data(230);
coeffs(41) = data(174) + data(186) + data(234);
coeffs(42) = data(62) + data(206) + data(242);
coeffs(43) = data(126) + data(222) + data(246);
coeffs(44) = data(190) + data(238) + data(250);
coeffs(45) = data(254);
coeffs(46) = -data(258);
coeffs(47) = -data(262);
coeffs(48) = -data(266);
coeffs(49) = -data(270);
coeffs(50) = data(3);
coeffs(51) = data(7) + data(19) + data(67);
coeffs(52) = data(23) + data(71) + data(83);
coeffs(53) = data(87);
coeffs(54) = data(11) + data(35) + data(131);
coeffs(55) = data(27) + data(39) + data(75) + data(99) + data(135) + data(147);
coeffs(56) = data(91) + data(103) + data(151);
coeffs(57) = data(43) + data(139) + data(163);
coeffs(58) = data(107) + data(155) + data(167);
coeffs(59) = data(171);
coeffs(60) = data(15) + data(51) + data(195);
coeffs(61) = data(31) + data(55) + data(79) + data(115) + data(199) + data(211);
coeffs(62) = data(95) + data(119) + data(215);
coeffs(63) = data(47) + data(59) + data(143) + data(179) + data(203) + data(227);
coeffs(64) = data(111) + data(123) + data(159) + data(183) + data(219) + data(231);
coeffs(65) = data(175) + data(187) + data(235);
coeffs(66) = data(63) + data(207) + data(243);
coeffs(67) = data(127) + data(223) + data(247);
coeffs(68) = data(191) + data(239) + data(251);
coeffs(69) = data(255);
coeffs(70) = -data(259);
coeffs(71) = -data(263);
coeffs(72) = -data(267);
coeffs(73) = -data(271);
coeffs(74) = data(4);
coeffs(75) = data(8) + data(20) + data(68);
coeffs(76) = data(24) + data(72) + data(84);
coeffs(77) = data(88);
coeffs(78) = data(12) + data(36) + data(132);
coeffs(79) = data(28) + data(40) + data(76) + data(100) + data(136) + data(148);
coeffs(80) = data(92) + data(104) + data(152);
coeffs(81) = data(44) + data(140) + data(164);
coeffs(82) = data(108) + data(156) + data(168);
coeffs(83) = data(172);
coeffs(84) = data(16) + data(52) + data(196);
coeffs(85) = data(32) + data(56) + data(80) + data(116) + data(200) + data(212);
coeffs(86) = data(96) + data(120) + data(216);
coeffs(87) = data(48) + data(60) + data(144) + data(180) + data(204) + data(228);
coeffs(88) = data(112) + data(124) + data(160) + data(184) + data(220) + data(232);
coeffs(89) = data(176) + data(188) + data(236);
coeffs(90) = data(64) + data(208) + data(244);
coeffs(91) = data(128) + data(224) + data(248);
coeffs(92) = data(192) + data(240) + data(252);
coeffs(93) = data(256);
coeffs(94) = -data(260);
coeffs(95) = -data(264);
coeffs(96) = -data(268);
coeffs(97) = -data(272);
coeffs(98) = -1;
end
