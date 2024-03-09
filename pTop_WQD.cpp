// LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems,
//          It also gives highly accurate uncertainty description of the solutions.
//
// Author: Jin Wu
// Affiliation: Hong Kong University of Science and Technology (HKUST)
// Emails: jin_wu_uestc@hotmail.com; jwucp@connect.ust.hk
// Reference: Wu, J., et al. (2022) Quadratic Pose Estimation Problems: 
//                                  Globally Optimal Solutions, 
//                                  Solvability/Observability Analysis,
//                                  and Uncertainty Description.
//                                  IEEE Transactions on Robotics.
//                                  https://doi.org/10.1109/TRO.2022.3155880
//
//
// pTop_WQD.cpp: Functions for computing matrices of Point-to-Plane registration problems

#include "pTop_WQD.h"
#include "utils.h"
#include <numeric>
#include <thread>

void mixed_pTop_func(Eigen::Matrix<double, 1, 85> &coef_J_pure,
                     const Eigen::Matrix<double, 9, 1> &pack) {
    double r1 = pack(0);
    double r2 = pack(1);
    double r3 = pack(2);
    double b1 = pack(3);
    double b2 = pack(4);
    double b3 = pack(5);
    double n1 = pack(6);
    double n2 = pack(7);
    double n3 = pack(8);

    double s_l_1_r_ = n3 * n3;
    double s_l_2_r_ = n2 * n2;
    double s_l_3_r_ = n1 * n1;
    double s_l_7_r_ = -n1 * r1 + n2 * r2 + n3 * r3;
    double s_l_8_r_ = n1 * r1 - n2 * r2 + n3 * r3;
    double s_l_9_r_ = b1 * n1 + b2 * n2 + b3 * n3;
    double s_l_10_r_ = n1 * r1 + n2 * r2 - n3 * r3;
    double s_l_12_r_ = n1 * r1 + n2 * r2 + n3 * r3;
    double s_l_13_r_ = n1 * r2 * 4.0 - n2 * r1 * 4.0;
    double s_l_14_r_ = n2 * r3 * 2.0 + n3 * r2 * 2.0;
    double s_l_15_r_ = n1 * r2 * 4.0 + n2 * r1 * 4.0;
    double s_l_16_r_ = n1 * r3 * 4.0 - n3 * r1 * 4.0;
    double s_l_17_r_ = n2 * r3 * 2.0 - n3 * r2 * 2.0;
    double s_l_24_r_ = n1 * r3 * 2.0 + n3 * r1 * 2.0;
    double s_l_25_r_ = n1 * r3 * 2.0 - n3 * r1 * 2.0;
    double s_l_48_r_ = n1 * r2 * 2.0 + n2 * r1 * 2.0;
    double s_l_54_r_ = n1 * r2 * 2.0 - n2 * r1 * 2.0;

    coef_J_pure(0, 0) = s_l_12_r_ * s_l_12_r_;
    coef_J_pure(0, 1) = s_l_12_r_ * s_l_17_r_ * 2.0;
    coef_J_pure(0, 2) = -s_l_12_r_ * s_l_16_r_;
    coef_J_pure(0, 3) = s_l_12_r_ * s_l_13_r_;
    coef_J_pure(0, 4) = s_l_7_r_ * s_l_12_r_ * -2.0 + s_l_17_r_ * s_l_17_r_;
    coef_J_pure(0, 5) = s_l_12_r_ * s_l_15_r_ - s_l_16_r_ * s_l_17_r_;
    coef_J_pure(0, 6) = s_l_13_r_ * s_l_17_r_ + s_l_12_r_ * s_l_24_r_ * 2.0;
    coef_J_pure(0, 7) = s_l_8_r_ * s_l_12_r_ * -2.0 + s_l_25_r_ * s_l_25_r_;
    coef_J_pure(0, 8) = s_l_12_r_ * s_l_14_r_ * 2.0 - s_l_13_r_ * s_l_25_r_;
    coef_J_pure(0, 9) = s_l_10_r_ * s_l_12_r_ * -2.0 + s_l_54_r_ * s_l_54_r_;
    coef_J_pure(0, 10) = n1 * s_l_12_r_ * 2.0;
    coef_J_pure(0, 11) = n2 * s_l_12_r_ * 2.0;
    coef_J_pure(0, 12) = n3 * s_l_12_r_ * 2.0;
    coef_J_pure(0, 13) = s_l_9_r_ * s_l_12_r_ * -2.0;
    coef_J_pure(0, 14) = s_l_7_r_ * s_l_17_r_ * -2.0;
    coef_J_pure(0, 15) = s_l_7_r_ * s_l_16_r_ + s_l_15_r_ * s_l_17_r_;
    coef_J_pure(0, 16) = -s_l_7_r_ * s_l_13_r_ + s_l_17_r_ * s_l_24_r_ * 2.0;
    coef_J_pure(0, 17) = s_l_8_r_ * s_l_17_r_ * -2.0 - s_l_15_r_ * s_l_25_r_;
    coef_J_pure(0, 18) = s_l_14_r_ * s_l_17_r_ * 2.0 - s_l_16_r_ * s_l_24_r_ + s_l_13_r_ * s_l_48_r_;
    coef_J_pure(0, 19) = s_l_10_r_ * s_l_17_r_ * -2.0 + s_l_13_r_ * s_l_24_r_;
    coef_J_pure(0, 20) = n1 * s_l_17_r_ * 2.0;
    coef_J_pure(0, 21) = n2 * s_l_17_r_ * 2.0;
    coef_J_pure(0, 22) = n3 * s_l_17_r_ * 2.0;
    coef_J_pure(0, 23) = s_l_9_r_ * s_l_17_r_ * -2.0;
    coef_J_pure(0, 24) = s_l_8_r_ * s_l_16_r_;
    coef_J_pure(0, 25) = -s_l_8_r_ * s_l_13_r_ - s_l_14_r_ * s_l_16_r_;
    coef_J_pure(0, 26) = s_l_10_r_ * s_l_16_r_ + s_l_13_r_ * s_l_14_r_;
    coef_J_pure(0, 27) = n1 * s_l_25_r_ * -2.0;
    coef_J_pure(0, 28) = n2 * s_l_25_r_ * -2.0;
    coef_J_pure(0, 29) = n3 * s_l_25_r_ * -2.0;
    coef_J_pure(0, 30) = s_l_9_r_ * s_l_16_r_;
    coef_J_pure(0, 31) = -s_l_10_r_ * s_l_13_r_;
    coef_J_pure(0, 32) = n1 * s_l_54_r_ * 2.0;
    coef_J_pure(0, 33) = n2 * s_l_54_r_ * 2.0;
    coef_J_pure(0, 34) = n3 * s_l_54_r_ * 2.0;
    coef_J_pure(0, 35) = -s_l_9_r_ * s_l_13_r_;
    coef_J_pure(0, 36) = s_l_7_r_ * s_l_7_r_;
    coef_J_pure(0, 37) = -s_l_7_r_ * s_l_15_r_;
    coef_J_pure(0, 38) = s_l_7_r_ * s_l_24_r_ * -2.0;
    coef_J_pure(0, 39) = s_l_7_r_ * s_l_8_r_ * 2.0 + s_l_48_r_ * s_l_48_r_;
    coef_J_pure(0, 40) = s_l_7_r_ * s_l_14_r_ * -2.0 + s_l_15_r_ * s_l_24_r_;
    coef_J_pure(0, 41) = s_l_7_r_ * s_l_10_r_ * 2.0 + s_l_24_r_ * s_l_24_r_;
    coef_J_pure(0, 42) = n1 * s_l_7_r_ * -2.0;
    coef_J_pure(0, 43) = n2 * s_l_7_r_ * -2.0;
    coef_J_pure(0, 44) = n3 * s_l_7_r_ * -2.0;
    coef_J_pure(0, 45) = s_l_7_r_ * s_l_9_r_ * 2.0;
    coef_J_pure(0, 46) = -s_l_8_r_ * s_l_15_r_;
    coef_J_pure(0, 47) = s_l_14_r_ * s_l_15_r_ - s_l_8_r_ * s_l_24_r_ * 2.0;
    coef_J_pure(0, 48) = -s_l_10_r_ * s_l_15_r_ + s_l_14_r_ * s_l_24_r_ * 2.0;
    coef_J_pure(0, 49) = n1 * s_l_48_r_ * 2.0;
    coef_J_pure(0, 50) = n2 * s_l_48_r_ * 2.0;
    coef_J_pure(0, 51) = n3 * s_l_48_r_ * 2.0;
    coef_J_pure(0, 52) = -s_l_9_r_ * s_l_15_r_;
    coef_J_pure(0, 53) = s_l_10_r_ * s_l_24_r_ * -2.0;
    coef_J_pure(0, 54) = n1 * s_l_24_r_ * 2.0;
    coef_J_pure(0, 55) = n2 * s_l_24_r_ * 2.0;
    coef_J_pure(0, 56) = n3 * s_l_24_r_ * 2.0;
    coef_J_pure(0, 57) = s_l_9_r_ * s_l_24_r_ * -2.0;
    coef_J_pure(0, 58) = s_l_8_r_ * s_l_8_r_;
    coef_J_pure(0, 59) = s_l_8_r_ * s_l_14_r_ * -2.0;
    coef_J_pure(0, 60) = s_l_8_r_ * s_l_10_r_ * 2.0 + s_l_14_r_ * s_l_14_r_;
    coef_J_pure(0, 61) = n1 * s_l_8_r_ * -2.0;
    coef_J_pure(0, 62) = n2 * s_l_8_r_ * -2.0;
    coef_J_pure(0, 63) = n3 * s_l_8_r_ * -2.0;
    coef_J_pure(0, 64) = s_l_8_r_ * s_l_9_r_ * 2.0;
    coef_J_pure(0, 65) = s_l_10_r_ * s_l_14_r_ * -2.0;
    coef_J_pure(0, 66) = n1 * s_l_14_r_ * 2.0;
    coef_J_pure(0, 67) = n2 * s_l_14_r_ * 2.0;
    coef_J_pure(0, 68) = n3 * s_l_14_r_ * 2.0;
    coef_J_pure(0, 69) = s_l_9_r_ * s_l_14_r_ * -2.0;
    coef_J_pure(0, 70) = s_l_10_r_ * s_l_10_r_;
    coef_J_pure(0, 71) = n1 * s_l_10_r_ * -2.0;
    coef_J_pure(0, 72) = n2 * s_l_10_r_ * -2.0;
    coef_J_pure(0, 73) = n3 * s_l_10_r_ * -2.0;
    coef_J_pure(0, 74) = s_l_9_r_ * s_l_10_r_ * 2.0;
    coef_J_pure(0, 75) = s_l_3_r_;
    coef_J_pure(0, 76) = n1 * n2 * 2.0;
    coef_J_pure(0, 77) = n1 * n3 * 2.0;
    coef_J_pure(0, 78) = n1 * s_l_9_r_ * -2.0;
    coef_J_pure(0, 79) = s_l_2_r_;
    coef_J_pure(0, 80) = n2 * n3 * 2.0;
    coef_J_pure(0, 81) = n2 * s_l_9_r_ * -2.0;
    coef_J_pure(0, 82) = s_l_1_r_;
    coef_J_pure(0, 83) = n3 * s_l_9_r_ * -2.0;
    coef_J_pure(0, 84) = s_l_9_r_ * s_l_9_r_;
}


void mixed2_pTop_func(Eigen::Matrix<double, 3, 3> &G,
                      Eigen::Matrix<double, 1, 11> &coeftq1,
                      Eigen::Matrix<double, 1, 11> &coeftq2,
                      Eigen::Matrix<double, 1, 11> &coeftq3,
                      Eigen::Matrix<double, 1, 36> &coef_Jacob1_qt,
                      Eigen::Matrix<double, 1, 36> &coef_Jacob2_qt,
                      Eigen::Matrix<double, 1, 36> &coef_Jacob3_qt,
                      Eigen::Matrix<double, 1, 36> &coef_Jacob4_qt,
                      const Eigen::Matrix<double, 1, 85> &coef_J_pure) {
    double coef_J1 = coef_J_pure(0);
    double coef_J2 = coef_J_pure(1);
    double coef_J3 = coef_J_pure(2);
    double coef_J4 = coef_J_pure(3);
    double coef_J5 = coef_J_pure(4);
    double coef_J6 = coef_J_pure(5);
    double coef_J7 = coef_J_pure(6);
    double coef_J8 = coef_J_pure(7);
    double coef_J9 = coef_J_pure(8);
    double coef_J10 = coef_J_pure(9);
    double coef_J11 = coef_J_pure(10);
    double coef_J12 = coef_J_pure(11);
    double coef_J13 = coef_J_pure(12);
    double coef_J14 = coef_J_pure(13);
    double coef_J15 = coef_J_pure(14);
    double coef_J16 = coef_J_pure(15);
    double coef_J17 = coef_J_pure(16);
    double coef_J18 = coef_J_pure(17);
    double coef_J19 = coef_J_pure(18);
    double coef_J20 = coef_J_pure(19);
    double coef_J21 = coef_J_pure(20);
    double coef_J22 = coef_J_pure(21);
    double coef_J23 = coef_J_pure(22);
    double coef_J24 = coef_J_pure(23);
    double coef_J25 = coef_J_pure(24);
    double coef_J26 = coef_J_pure(25);
    double coef_J27 = coef_J_pure(26);
    double coef_J28 = coef_J_pure(27);
    double coef_J29 = coef_J_pure(28);
    double coef_J30 = coef_J_pure(29);
    double coef_J31 = coef_J_pure(30);
    double coef_J32 = coef_J_pure(31);
    double coef_J33 = coef_J_pure(32);
    double coef_J34 = coef_J_pure(33);
    double coef_J35 = coef_J_pure(34);
    double coef_J36 = coef_J_pure(35);
    double coef_J37 = coef_J_pure(36);
    double coef_J38 = coef_J_pure(37);
    double coef_J39 = coef_J_pure(38);
    double coef_J40 = coef_J_pure(39);
    double coef_J41 = coef_J_pure(40);
    double coef_J42 = coef_J_pure(41);
    double coef_J43 = coef_J_pure(42);
    double coef_J44 = coef_J_pure(43);
    double coef_J45 = coef_J_pure(44);
    double coef_J46 = coef_J_pure(45);
    double coef_J47 = coef_J_pure(46);
    double coef_J48 = coef_J_pure(47);
    double coef_J49 = coef_J_pure(48);
    double coef_J50 = coef_J_pure(49);
    double coef_J51 = coef_J_pure(50);
    double coef_J52 = coef_J_pure(51);
    double coef_J53 = coef_J_pure(52);
    double coef_J54 = coef_J_pure(53);
    double coef_J55 = coef_J_pure(54);
    double coef_J56 = coef_J_pure(55);
    double coef_J57 = coef_J_pure(56);
    double coef_J58 = coef_J_pure(57);
    double coef_J59 = coef_J_pure(58);
    double coef_J60 = coef_J_pure(59);
    double coef_J61 = coef_J_pure(60);
    double coef_J62 = coef_J_pure(61);
    double coef_J63 = coef_J_pure(62);
    double coef_J64 = coef_J_pure(63);
    double coef_J65 = coef_J_pure(64);
    double coef_J66 = coef_J_pure(65);
    double coef_J67 = coef_J_pure(66);
    double coef_J68 = coef_J_pure(67);
    double coef_J69 = coef_J_pure(68);
    double coef_J70 = coef_J_pure(69);
    double coef_J71 = coef_J_pure(70);
    double coef_J72 = coef_J_pure(71);
    double coef_J73 = coef_J_pure(72);
    double coef_J74 = coef_J_pure(73);
    double coef_J75 = coef_J_pure(74);
    double coef_J76 = coef_J_pure(75);
    double coef_J77 = coef_J_pure(76);
    double coef_J78 = coef_J_pure(77);
    double coef_J79 = coef_J_pure(78);
    double coef_J80 = coef_J_pure(79);
    double coef_J81 = coef_J_pure(80);
    double coef_J82 = coef_J_pure(81);
    double coef_J83 = coef_J_pure(82);
    double coef_J84 = coef_J_pure(83);
    double coef_J85 = coef_J_pure(84);

    G(0, 0) = coef_J76 * -2.0;
    G(0, 1) = -coef_J77;
    G(0, 2) = -coef_J78;
    G(1, 0) = -coef_J77;
    G(1, 1) = coef_J80 * -2.0;
    G(1, 2) = -coef_J81;
    G(2, 0) = -coef_J78;
    G(2, 1) = -coef_J81;
    G(2, 2) = coef_J83 * -2.0;

    coeftq1(0, 0) = coef_J11;
    coeftq1(0, 1) = coef_J21;
    coeftq1(0, 2) = coef_J28;
    coeftq1(0, 3) = coef_J33;
    coeftq1(0, 4) = coef_J43;
    coeftq1(0, 5) = coef_J50;
    coeftq1(0, 6) = coef_J55;
    coeftq1(0, 7) = coef_J62;
    coeftq1(0, 8) = coef_J67;
    coeftq1(0, 9) = coef_J72;
    coeftq1(0, 10) = coef_J79;

    coeftq2(0, 0) = coef_J12;
    coeftq2(0, 1) = coef_J22;
    coeftq2(0, 2) = coef_J29;
    coeftq2(0, 3) = coef_J34;
    coeftq2(0, 4) = coef_J44;
    coeftq2(0, 5) = coef_J51;
    coeftq2(0, 6) = coef_J56;
    coeftq2(0, 7) = coef_J63;
    coeftq2(0, 8) = coef_J68;
    coeftq2(0, 9) = coef_J73;
    coeftq2(0, 10) = coef_J82;

    coeftq3(0, 0) = coef_J13;
    coeftq3(0, 1) = coef_J23;
    coeftq3(0, 2) = coef_J30;
    coeftq3(0, 3) = coef_J35;
    coeftq3(0, 4) = coef_J45;
    coeftq3(0, 5) = coef_J52;
    coeftq3(0, 6) = coef_J57;
    coeftq3(0, 7) = coef_J64;
    coeftq3(0, 8) = coef_J69;
    coeftq3(0, 9) = coef_J74;
    coeftq3(0, 10) = coef_J84;

    coef_Jacob1_qt(0, 0) = coef_J1 * 4.0;
    coef_Jacob1_qt(0, 1) = coef_J2 * 3.0;
    coef_Jacob1_qt(0, 2) = coef_J3 * 3.0;
    coef_Jacob1_qt(0, 3) = coef_J4 * 3.0;
    coef_Jacob1_qt(0, 4) = coef_J5 * 2.0;
    coef_Jacob1_qt(0, 5) = coef_J6 * 2.0;
    coef_Jacob1_qt(0, 6) = coef_J7 * 2.0;
    coef_Jacob1_qt(0, 7) = coef_J8 * 2.0;
    coef_Jacob1_qt(0, 8) = coef_J9 * 2.0;
    coef_Jacob1_qt(0, 9) = coef_J10 * 2.0;
    coef_Jacob1_qt(0, 10) = coef_J11 * 2.0;
    coef_Jacob1_qt(0, 11) = coef_J12 * 2.0;
    coef_Jacob1_qt(0, 12) = coef_J13 * 2.0;
    coef_Jacob1_qt(0, 13) = coef_J14 * 2.0;
    coef_Jacob1_qt(0, 14) = coef_J15;
    coef_Jacob1_qt(0, 15) = coef_J16;
    coef_Jacob1_qt(0, 16) = coef_J17;
    coef_Jacob1_qt(0, 17) = coef_J18;
    coef_Jacob1_qt(0, 18) = coef_J19;
    coef_Jacob1_qt(0, 19) = coef_J20;
    coef_Jacob1_qt(0, 20) = coef_J21;
    coef_Jacob1_qt(0, 21) = coef_J22;
    coef_Jacob1_qt(0, 22) = coef_J23;
    coef_Jacob1_qt(0, 23) = coef_J24;
    coef_Jacob1_qt(0, 24) = coef_J25;
    coef_Jacob1_qt(0, 25) = coef_J26;
    coef_Jacob1_qt(0, 26) = coef_J27;
    coef_Jacob1_qt(0, 27) = coef_J28;
    coef_Jacob1_qt(0, 28) = coef_J29;
    coef_Jacob1_qt(0, 29) = coef_J30;
    coef_Jacob1_qt(0, 30) = coef_J31;
    coef_Jacob1_qt(0, 31) = coef_J32;
    coef_Jacob1_qt(0, 32) = coef_J33;
    coef_Jacob1_qt(0, 33) = coef_J34;
    coef_Jacob1_qt(0, 34) = coef_J35;
    coef_Jacob1_qt(0, 35) = coef_J36;

    coef_Jacob2_qt(0, 0) = coef_J2;
    coef_Jacob2_qt(0, 1) = coef_J5 * 2.0;
    coef_Jacob2_qt(0, 2) = coef_J6;
    coef_Jacob2_qt(0, 3) = coef_J7;
    coef_Jacob2_qt(0, 4) = coef_J15 * 3.0;
    coef_Jacob2_qt(0, 5) = coef_J16 * 2.0;
    coef_Jacob2_qt(0, 6) = coef_J17 * 2.0;
    coef_Jacob2_qt(0, 7) = coef_J18;
    coef_Jacob2_qt(0, 8) = coef_J19;
    coef_Jacob2_qt(0, 9) = coef_J20;
    coef_Jacob2_qt(0, 10) = coef_J21;
    coef_Jacob2_qt(0, 11) = coef_J22;
    coef_Jacob2_qt(0, 12) = coef_J23;
    coef_Jacob2_qt(0, 13) = coef_J24;
    coef_Jacob2_qt(0, 14) = coef_J37 * 4.0;
    coef_Jacob2_qt(0, 15) = coef_J38 * 3.0;
    coef_Jacob2_qt(0, 16) = coef_J39 * 3.0;
    coef_Jacob2_qt(0, 17) = coef_J40 * 2.0;
    coef_Jacob2_qt(0, 18) = coef_J41 * 2.0;
    coef_Jacob2_qt(0, 19) = coef_J42 * 2.0;
    coef_Jacob2_qt(0, 20) = coef_J43 * 2.0;
    coef_Jacob2_qt(0, 21) = coef_J44 * 2.0;
    coef_Jacob2_qt(0, 22) = coef_J45 * 2.0;
    coef_Jacob2_qt(0, 23) = coef_J46 * 2.0;
    coef_Jacob2_qt(0, 24) = coef_J47;
    coef_Jacob2_qt(0, 25) = coef_J48;
    coef_Jacob2_qt(0, 26) = coef_J49;
    coef_Jacob2_qt(0, 27) = coef_J50;
    coef_Jacob2_qt(0, 28) = coef_J51;
    coef_Jacob2_qt(0, 29) = coef_J52;
    coef_Jacob2_qt(0, 30) = coef_J53;
    coef_Jacob2_qt(0, 31) = coef_J54;
    coef_Jacob2_qt(0, 32) = coef_J55;
    coef_Jacob2_qt(0, 33) = coef_J56;
    coef_Jacob2_qt(0, 34) = coef_J57;
    coef_Jacob2_qt(0, 35) = coef_J58;

    coef_Jacob3_qt(0, 0) = coef_J3;
    coef_Jacob3_qt(0, 1) = coef_J6;
    coef_Jacob3_qt(0, 2) = coef_J8 * 2.0;
    coef_Jacob3_qt(0, 3) = coef_J9;
    coef_Jacob3_qt(0, 4) = coef_J16;
    coef_Jacob3_qt(0, 5) = coef_J18 * 2.0;
    coef_Jacob3_qt(0, 6) = coef_J19;
    coef_Jacob3_qt(0, 7) = coef_J25 * 3.0;
    coef_Jacob3_qt(0, 8) = coef_J26 * 2.0;
    coef_Jacob3_qt(0, 9) = coef_J27;
    coef_Jacob3_qt(0, 10) = coef_J28;
    coef_Jacob3_qt(0, 11) = coef_J29;
    coef_Jacob3_qt(0, 12) = coef_J30;
    coef_Jacob3_qt(0, 13) = coef_J31;
    coef_Jacob3_qt(0, 14) = coef_J38;
    coef_Jacob3_qt(0, 15) = coef_J40 * 2.0;
    coef_Jacob3_qt(0, 16) = coef_J41;
    coef_Jacob3_qt(0, 17) = coef_J47 * 3.0;
    coef_Jacob3_qt(0, 18) = coef_J48 * 2.0;
    coef_Jacob3_qt(0, 19) = coef_J49;
    coef_Jacob3_qt(0, 20) = coef_J50;
    coef_Jacob3_qt(0, 21) = coef_J51;
    coef_Jacob3_qt(0, 22) = coef_J52;
    coef_Jacob3_qt(0, 23) = coef_J53;
    coef_Jacob3_qt(0, 24) = coef_J59 * 4.0;
    coef_Jacob3_qt(0, 25) = coef_J60 * 3.0;
    coef_Jacob3_qt(0, 26) = coef_J61 * 2.0;
    coef_Jacob3_qt(0, 27) = coef_J62 * 2.0;
    coef_Jacob3_qt(0, 28) = coef_J63 * 2.0;
    coef_Jacob3_qt(0, 29) = coef_J64 * 2.0;
    coef_Jacob3_qt(0, 30) = coef_J65 * 2.0;
    coef_Jacob3_qt(0, 31) = coef_J66;
    coef_Jacob3_qt(0, 32) = coef_J67;
    coef_Jacob3_qt(0, 33) = coef_J68;
    coef_Jacob3_qt(0, 34) = coef_J69;
    coef_Jacob3_qt(0, 35) = coef_J70;

    coef_Jacob4_qt(0, 0) = coef_J4;
    coef_Jacob4_qt(0, 1) = coef_J7;
    coef_Jacob4_qt(0, 2) = coef_J9;
    coef_Jacob4_qt(0, 3) = coef_J10 * 2.0;
    coef_Jacob4_qt(0, 4) = coef_J17;
    coef_Jacob4_qt(0, 5) = coef_J19;
    coef_Jacob4_qt(0, 6) = coef_J20 * 2.0;
    coef_Jacob4_qt(0, 7) = coef_J26;
    coef_Jacob4_qt(0, 8) = coef_J27 * 2.0;
    coef_Jacob4_qt(0, 9) = coef_J32 * 3.0;
    coef_Jacob4_qt(0, 10) = coef_J33;
    coef_Jacob4_qt(0, 11) = coef_J34;
    coef_Jacob4_qt(0, 12) = coef_J35;
    coef_Jacob4_qt(0, 13) = coef_J36;
    coef_Jacob4_qt(0, 14) = coef_J39;
    coef_Jacob4_qt(0, 15) = coef_J41;
    coef_Jacob4_qt(0, 16) = coef_J42 * 2.0;
    coef_Jacob4_qt(0, 17) = coef_J48;
    coef_Jacob4_qt(0, 18) = coef_J49 * 2.0;
    coef_Jacob4_qt(0, 19) = coef_J54 * 3.0;
    coef_Jacob4_qt(0, 20) = coef_J55;
    coef_Jacob4_qt(0, 21) = coef_J56;
    coef_Jacob4_qt(0, 22) = coef_J57;
    coef_Jacob4_qt(0, 23) = coef_J58;
    coef_Jacob4_qt(0, 24) = coef_J60;
    coef_Jacob4_qt(0, 25) = coef_J61 * 2.0;
    coef_Jacob4_qt(0, 26) = coef_J66 * 3.0;
    coef_Jacob4_qt(0, 27) = coef_J67;
    coef_Jacob4_qt(0, 28) = coef_J68;
    coef_Jacob4_qt(0, 29) = coef_J69;
    coef_Jacob4_qt(0, 30) = coef_J70;
    coef_Jacob4_qt(0, 31) = coef_J71 * 4.0;
    coef_Jacob4_qt(0, 32) = coef_J72 * 2.0;
    coef_Jacob4_qt(0, 33) = coef_J73 * 2.0;
    coef_Jacob4_qt(0, 34) = coef_J74 * 2.0;
    coef_Jacob4_qt(0, 35) = coef_J75 * 2.0;

}


void mixed3_pTop_func(Eigen::Matrix<double, 4, 24> &coef_f_q_sym,
                      Eigen::Matrix<double, 4, 64> &W,
                      Eigen::Matrix<double, 4, 4> &Q,
                      const Eigen::Matrix<double, 3, 3> &pinvG,
                      const Eigen::Matrix<double, 3, 11> &coefs_tq,
                      const Eigen::Matrix<double, 4, 36> &coef_Jacob_qt_syms) {
    double pinvG1_1 = pinvG(0, 0);
    double pinvG1_2 = pinvG(0, 1);
    double pinvG1_3 = pinvG(0, 2);
    double pinvG2_1 = pinvG(1, 0);
    double pinvG2_2 = pinvG(1, 1);
    double pinvG2_3 = pinvG(1, 2);
    double pinvG3_1 = pinvG(2, 0);
    double pinvG3_2 = pinvG(2, 1);
    double pinvG3_3 = pinvG(2, 2);
    double coefs_tq1_1 = coefs_tq(0, 0);
    double coefs_tq1_2 = coefs_tq(0, 1);
    double coefs_tq1_3 = coefs_tq(0, 2);
    double coefs_tq1_4 = coefs_tq(0, 3);
    double coefs_tq1_5 = coefs_tq(0, 4);
    double coefs_tq1_6 = coefs_tq(0, 5);
    double coefs_tq1_7 = coefs_tq(0, 6);
    double coefs_tq1_8 = coefs_tq(0, 7);
    double coefs_tq1_9 = coefs_tq(0, 8);
    double coefs_tq1_10 = coefs_tq(0, 9);
    double coefs_tq1_11 = coefs_tq(0, 10);
    double coefs_tq2_1 = coefs_tq(1, 0);
    double coefs_tq2_2 = coefs_tq(1, 1);
    double coefs_tq2_3 = coefs_tq(1, 2);
    double coefs_tq2_4 = coefs_tq(1, 3);
    double coefs_tq2_5 = coefs_tq(1, 4);
    double coefs_tq2_6 = coefs_tq(1, 5);
    double coefs_tq2_7 = coefs_tq(1, 6);
    double coefs_tq2_8 = coefs_tq(1, 7);
    double coefs_tq2_9 = coefs_tq(1, 8);
    double coefs_tq2_10 = coefs_tq(1, 9);
    double coefs_tq2_11 = coefs_tq(1, 10);
    double coefs_tq3_1 = coefs_tq(2, 0);
    double coefs_tq3_2 = coefs_tq(2, 1);
    double coefs_tq3_3 = coefs_tq(2, 2);
    double coefs_tq3_4 = coefs_tq(2, 3);
    double coefs_tq3_5 = coefs_tq(2, 4);
    double coefs_tq3_6 = coefs_tq(2, 5);
    double coefs_tq3_7 = coefs_tq(2, 6);
    double coefs_tq3_8 = coefs_tq(2, 7);
    double coefs_tq3_9 = coefs_tq(2, 8);
    double coefs_tq3_10 = coefs_tq(2, 9);
    double coefs_tq3_11 = coefs_tq(2, 10);
    double coef_Jacob1_qt_syms1 = coef_Jacob_qt_syms(0, 0);
    double coef_Jacob1_qt_syms2 = coef_Jacob_qt_syms(0, 1);
    double coef_Jacob1_qt_syms3 = coef_Jacob_qt_syms(0, 2);
    double coef_Jacob1_qt_syms4 = coef_Jacob_qt_syms(0, 3);
    double coef_Jacob1_qt_syms5 = coef_Jacob_qt_syms(0, 4);
    double coef_Jacob1_qt_syms6 = coef_Jacob_qt_syms(0, 5);
    double coef_Jacob1_qt_syms7 = coef_Jacob_qt_syms(0, 6);
    double coef_Jacob1_qt_syms8 = coef_Jacob_qt_syms(0, 7);
    double coef_Jacob1_qt_syms9 = coef_Jacob_qt_syms(0, 8);
    double coef_Jacob1_qt_syms10 = coef_Jacob_qt_syms(0, 9);
    double coef_Jacob1_qt_syms11 = coef_Jacob_qt_syms(0, 10);
    double coef_Jacob1_qt_syms12 = coef_Jacob_qt_syms(0, 11);
    double coef_Jacob1_qt_syms13 = coef_Jacob_qt_syms(0, 12);
    double coef_Jacob1_qt_syms14 = coef_Jacob_qt_syms(0, 13);
    double coef_Jacob1_qt_syms15 = coef_Jacob_qt_syms(0, 14);
    double coef_Jacob1_qt_syms16 = coef_Jacob_qt_syms(0, 15);
    double coef_Jacob1_qt_syms17 = coef_Jacob_qt_syms(0, 16);
    double coef_Jacob1_qt_syms18 = coef_Jacob_qt_syms(0, 17);
    double coef_Jacob1_qt_syms19 = coef_Jacob_qt_syms(0, 18);
    double coef_Jacob1_qt_syms20 = coef_Jacob_qt_syms(0, 19);
    double coef_Jacob1_qt_syms21 = coef_Jacob_qt_syms(0, 20);
    double coef_Jacob1_qt_syms22 = coef_Jacob_qt_syms(0, 21);
    double coef_Jacob1_qt_syms23 = coef_Jacob_qt_syms(0, 22);
    double coef_Jacob1_qt_syms24 = coef_Jacob_qt_syms(0, 23);
    double coef_Jacob1_qt_syms25 = coef_Jacob_qt_syms(0, 24);
    double coef_Jacob1_qt_syms26 = coef_Jacob_qt_syms(0, 25);
    double coef_Jacob1_qt_syms27 = coef_Jacob_qt_syms(0, 26);
    double coef_Jacob1_qt_syms28 = coef_Jacob_qt_syms(0, 27);
    double coef_Jacob1_qt_syms29 = coef_Jacob_qt_syms(0, 28);
    double coef_Jacob1_qt_syms30 = coef_Jacob_qt_syms(0, 29);
    double coef_Jacob1_qt_syms31 = coef_Jacob_qt_syms(0, 30);
    double coef_Jacob1_qt_syms32 = coef_Jacob_qt_syms(0, 31);
    double coef_Jacob1_qt_syms33 = coef_Jacob_qt_syms(0, 32);
    double coef_Jacob1_qt_syms34 = coef_Jacob_qt_syms(0, 33);
    double coef_Jacob1_qt_syms35 = coef_Jacob_qt_syms(0, 34);
    double coef_Jacob1_qt_syms36 = coef_Jacob_qt_syms(0, 35);
    double coef_Jacob2_qt_syms1 = coef_Jacob_qt_syms(1, 0);
    double coef_Jacob2_qt_syms2 = coef_Jacob_qt_syms(1, 1);
    double coef_Jacob2_qt_syms3 = coef_Jacob_qt_syms(1, 2);
    double coef_Jacob2_qt_syms4 = coef_Jacob_qt_syms(1, 3);
    double coef_Jacob2_qt_syms5 = coef_Jacob_qt_syms(1, 4);
    double coef_Jacob2_qt_syms6 = coef_Jacob_qt_syms(1, 5);
    double coef_Jacob2_qt_syms7 = coef_Jacob_qt_syms(1, 6);
    double coef_Jacob2_qt_syms8 = coef_Jacob_qt_syms(1, 7);
    double coef_Jacob2_qt_syms9 = coef_Jacob_qt_syms(1, 8);
    double coef_Jacob2_qt_syms10 = coef_Jacob_qt_syms(1, 9);
    double coef_Jacob2_qt_syms11 = coef_Jacob_qt_syms(1, 10);
    double coef_Jacob2_qt_syms12 = coef_Jacob_qt_syms(1, 11);
    double coef_Jacob2_qt_syms13 = coef_Jacob_qt_syms(1, 12);
    double coef_Jacob2_qt_syms14 = coef_Jacob_qt_syms(1, 13);
    double coef_Jacob2_qt_syms15 = coef_Jacob_qt_syms(1, 14);
    double coef_Jacob2_qt_syms16 = coef_Jacob_qt_syms(1, 15);
    double coef_Jacob2_qt_syms17 = coef_Jacob_qt_syms(1, 16);
    double coef_Jacob2_qt_syms18 = coef_Jacob_qt_syms(1, 17);
    double coef_Jacob2_qt_syms19 = coef_Jacob_qt_syms(1, 18);
    double coef_Jacob2_qt_syms20 = coef_Jacob_qt_syms(1, 19);
    double coef_Jacob2_qt_syms21 = coef_Jacob_qt_syms(1, 20);
    double coef_Jacob2_qt_syms22 = coef_Jacob_qt_syms(1, 21);
    double coef_Jacob2_qt_syms23 = coef_Jacob_qt_syms(1, 22);
    double coef_Jacob2_qt_syms24 = coef_Jacob_qt_syms(1, 23);
    double coef_Jacob2_qt_syms25 = coef_Jacob_qt_syms(1, 24);
    double coef_Jacob2_qt_syms26 = coef_Jacob_qt_syms(1, 25);
    double coef_Jacob2_qt_syms27 = coef_Jacob_qt_syms(1, 26);
    double coef_Jacob2_qt_syms28 = coef_Jacob_qt_syms(1, 27);
    double coef_Jacob2_qt_syms29 = coef_Jacob_qt_syms(1, 28);
    double coef_Jacob2_qt_syms30 = coef_Jacob_qt_syms(1, 29);
    double coef_Jacob2_qt_syms31 = coef_Jacob_qt_syms(1, 30);
    double coef_Jacob2_qt_syms32 = coef_Jacob_qt_syms(1, 31);
    double coef_Jacob2_qt_syms33 = coef_Jacob_qt_syms(1, 32);
    double coef_Jacob2_qt_syms34 = coef_Jacob_qt_syms(1, 33);
    double coef_Jacob2_qt_syms35 = coef_Jacob_qt_syms(1, 34);
    double coef_Jacob2_qt_syms36 = coef_Jacob_qt_syms(1, 35);
    double coef_Jacob3_qt_syms1 = coef_Jacob_qt_syms(2, 0);
    double coef_Jacob3_qt_syms2 = coef_Jacob_qt_syms(2, 1);
    double coef_Jacob3_qt_syms3 = coef_Jacob_qt_syms(2, 2);
    double coef_Jacob3_qt_syms4 = coef_Jacob_qt_syms(2, 3);
    double coef_Jacob3_qt_syms5 = coef_Jacob_qt_syms(2, 4);
    double coef_Jacob3_qt_syms6 = coef_Jacob_qt_syms(2, 5);
    double coef_Jacob3_qt_syms7 = coef_Jacob_qt_syms(2, 6);
    double coef_Jacob3_qt_syms8 = coef_Jacob_qt_syms(2, 7);
    double coef_Jacob3_qt_syms9 = coef_Jacob_qt_syms(2, 8);
    double coef_Jacob3_qt_syms10 = coef_Jacob_qt_syms(2, 9);
    double coef_Jacob3_qt_syms11 = coef_Jacob_qt_syms(2, 10);
    double coef_Jacob3_qt_syms12 = coef_Jacob_qt_syms(2, 11);
    double coef_Jacob3_qt_syms13 = coef_Jacob_qt_syms(2, 12);
    double coef_Jacob3_qt_syms14 = coef_Jacob_qt_syms(2, 13);
    double coef_Jacob3_qt_syms15 = coef_Jacob_qt_syms(2, 14);
    double coef_Jacob3_qt_syms16 = coef_Jacob_qt_syms(2, 15);
    double coef_Jacob3_qt_syms17 = coef_Jacob_qt_syms(2, 16);
    double coef_Jacob3_qt_syms18 = coef_Jacob_qt_syms(2, 17);
    double coef_Jacob3_qt_syms19 = coef_Jacob_qt_syms(2, 18);
    double coef_Jacob3_qt_syms20 = coef_Jacob_qt_syms(2, 19);
    double coef_Jacob3_qt_syms21 = coef_Jacob_qt_syms(2, 20);
    double coef_Jacob3_qt_syms22 = coef_Jacob_qt_syms(2, 21);
    double coef_Jacob3_qt_syms23 = coef_Jacob_qt_syms(2, 22);
    double coef_Jacob3_qt_syms24 = coef_Jacob_qt_syms(2, 23);
    double coef_Jacob3_qt_syms25 = coef_Jacob_qt_syms(2, 24);
    double coef_Jacob3_qt_syms26 = coef_Jacob_qt_syms(2, 25);
    double coef_Jacob3_qt_syms27 = coef_Jacob_qt_syms(2, 26);
    double coef_Jacob3_qt_syms28 = coef_Jacob_qt_syms(2, 27);
    double coef_Jacob3_qt_syms29 = coef_Jacob_qt_syms(2, 28);
    double coef_Jacob3_qt_syms30 = coef_Jacob_qt_syms(2, 29);
    double coef_Jacob3_qt_syms31 = coef_Jacob_qt_syms(2, 30);
    double coef_Jacob3_qt_syms32 = coef_Jacob_qt_syms(2, 31);
    double coef_Jacob3_qt_syms33 = coef_Jacob_qt_syms(2, 32);
    double coef_Jacob3_qt_syms34 = coef_Jacob_qt_syms(2, 33);
    double coef_Jacob3_qt_syms35 = coef_Jacob_qt_syms(2, 34);
    double coef_Jacob3_qt_syms36 = coef_Jacob_qt_syms(2, 35);
    double coef_Jacob4_qt_syms1 = coef_Jacob_qt_syms(3, 0);
    double coef_Jacob4_qt_syms2 = coef_Jacob_qt_syms(3, 1);
    double coef_Jacob4_qt_syms3 = coef_Jacob_qt_syms(3, 2);
    double coef_Jacob4_qt_syms4 = coef_Jacob_qt_syms(3, 3);
    double coef_Jacob4_qt_syms5 = coef_Jacob_qt_syms(3, 4);
    double coef_Jacob4_qt_syms6 = coef_Jacob_qt_syms(3, 5);
    double coef_Jacob4_qt_syms7 = coef_Jacob_qt_syms(3, 6);
    double coef_Jacob4_qt_syms8 = coef_Jacob_qt_syms(3, 7);
    double coef_Jacob4_qt_syms9 = coef_Jacob_qt_syms(3, 8);
    double coef_Jacob4_qt_syms10 = coef_Jacob_qt_syms(3, 9);
    double coef_Jacob4_qt_syms11 = coef_Jacob_qt_syms(3, 10);
    double coef_Jacob4_qt_syms12 = coef_Jacob_qt_syms(3, 11);
    double coef_Jacob4_qt_syms13 = coef_Jacob_qt_syms(3, 12);
    double coef_Jacob4_qt_syms14 = coef_Jacob_qt_syms(3, 13);
    double coef_Jacob4_qt_syms15 = coef_Jacob_qt_syms(3, 14);
    double coef_Jacob4_qt_syms16 = coef_Jacob_qt_syms(3, 15);
    double coef_Jacob4_qt_syms17 = coef_Jacob_qt_syms(3, 16);
    double coef_Jacob4_qt_syms18 = coef_Jacob_qt_syms(3, 17);
    double coef_Jacob4_qt_syms19 = coef_Jacob_qt_syms(3, 18);
    double coef_Jacob4_qt_syms20 = coef_Jacob_qt_syms(3, 19);
    double coef_Jacob4_qt_syms21 = coef_Jacob_qt_syms(3, 20);
    double coef_Jacob4_qt_syms22 = coef_Jacob_qt_syms(3, 21);
    double coef_Jacob4_qt_syms23 = coef_Jacob_qt_syms(3, 22);
    double coef_Jacob4_qt_syms24 = coef_Jacob_qt_syms(3, 23);
    double coef_Jacob4_qt_syms25 = coef_Jacob_qt_syms(3, 24);
    double coef_Jacob4_qt_syms26 = coef_Jacob_qt_syms(3, 25);
    double coef_Jacob4_qt_syms27 = coef_Jacob_qt_syms(3, 26);
    double coef_Jacob4_qt_syms28 = coef_Jacob_qt_syms(3, 27);
    double coef_Jacob4_qt_syms29 = coef_Jacob_qt_syms(3, 28);
    double coef_Jacob4_qt_syms30 = coef_Jacob_qt_syms(3, 29);
    double coef_Jacob4_qt_syms31 = coef_Jacob_qt_syms(3, 30);
    double coef_Jacob4_qt_syms32 = coef_Jacob_qt_syms(3, 31);
    double coef_Jacob4_qt_syms33 = coef_Jacob_qt_syms(3, 32);
    double coef_Jacob4_qt_syms34 = coef_Jacob_qt_syms(3, 33);
    double coef_Jacob4_qt_syms35 = coef_Jacob_qt_syms(3, 34);
    double coef_Jacob4_qt_syms36 = coef_Jacob_qt_syms(3, 35);

    coef_f_q_sym(0, 0) = coef_Jacob1_qt_syms1 + coefs_tq1_1 * coef_Jacob1_qt_syms11 * pinvG1_1 +
                         coefs_tq1_1 * coef_Jacob1_qt_syms12 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob1_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob1_qt_syms13 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob1_qt_syms12 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob1_qt_syms11 * pinvG1_3 +
                         coefs_tq2_1 * coef_Jacob1_qt_syms13 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob1_qt_syms12 * pinvG2_3 +
                         coefs_tq3_1 * coef_Jacob1_qt_syms13 * pinvG3_3;
    coef_f_q_sym(0, 1) = coef_Jacob1_qt_syms2 + coefs_tq1_2 * coef_Jacob1_qt_syms11 * pinvG1_1 +
                         coefs_tq1_1 * coef_Jacob1_qt_syms21 * pinvG1_1 +
                         coefs_tq1_2 * coef_Jacob1_qt_syms12 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob1_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob1_qt_syms22 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob1_qt_syms21 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob1_qt_syms13 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob1_qt_syms12 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob1_qt_syms11 * pinvG1_3 +
                         coefs_tq1_1 * coef_Jacob1_qt_syms23 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob1_qt_syms22 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob1_qt_syms21 * pinvG1_3 +
                         coefs_tq2_2 * coef_Jacob1_qt_syms13 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob1_qt_syms12 * pinvG2_3 +
                         coefs_tq2_1 * coef_Jacob1_qt_syms23 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob1_qt_syms22 * pinvG2_3 +
                         coefs_tq3_2 * coef_Jacob1_qt_syms13 * pinvG3_3 +
                         coefs_tq3_1 * coef_Jacob1_qt_syms23 * pinvG3_3;
    coef_f_q_sym(0, 2) = coef_Jacob1_qt_syms3 + coefs_tq1_3 * coef_Jacob1_qt_syms11 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob1_qt_syms12 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob1_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob1_qt_syms28 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob1_qt_syms13 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob1_qt_syms12 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob1_qt_syms11 * pinvG1_3 +
                         coefs_tq1_1 * coef_Jacob1_qt_syms29 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob1_qt_syms28 * pinvG1_2 +
                         coefs_tq2_3 * coef_Jacob1_qt_syms13 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob1_qt_syms12 * pinvG2_3 +
                         coefs_tq1_1 * coef_Jacob1_qt_syms30 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob1_qt_syms29 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob1_qt_syms28 * pinvG1_3 +
                         coefs_tq3_3 * coef_Jacob1_qt_syms13 * pinvG3_3 +
                         coefs_tq2_1 * coef_Jacob1_qt_syms30 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob1_qt_syms29 * pinvG2_3 +
                         coefs_tq3_1 * coef_Jacob1_qt_syms30 * pinvG3_3;
    coef_f_q_sym(0, 3) = coef_Jacob1_qt_syms4 + coefs_tq1_4 * coef_Jacob1_qt_syms11 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob1_qt_syms12 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob1_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob1_qt_syms33 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob1_qt_syms13 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob1_qt_syms12 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob1_qt_syms11 * pinvG1_3 +
                         coefs_tq1_1 * coef_Jacob1_qt_syms34 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob1_qt_syms33 * pinvG1_2 +
                         coefs_tq2_4 * coef_Jacob1_qt_syms13 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob1_qt_syms12 * pinvG2_3 +
                         coefs_tq1_1 * coef_Jacob1_qt_syms35 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob1_qt_syms34 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob1_qt_syms33 * pinvG1_3 +
                         coefs_tq3_4 * coef_Jacob1_qt_syms13 * pinvG3_3 +
                         coefs_tq2_1 * coef_Jacob1_qt_syms35 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob1_qt_syms34 * pinvG2_3 +
                         coefs_tq3_1 * coef_Jacob1_qt_syms35 * pinvG3_3;
    coef_f_q_sym(0, 4) = coef_Jacob1_qt_syms5 + coefs_tq1_5 * coef_Jacob1_qt_syms11 * pinvG1_1 +
                         coefs_tq1_2 * coef_Jacob1_qt_syms21 * pinvG1_1 +
                         coefs_tq1_5 * coef_Jacob1_qt_syms12 * pinvG2_1 +
                         coefs_tq2_5 * coef_Jacob1_qt_syms11 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob1_qt_syms22 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob1_qt_syms21 * pinvG1_2 +
                         coefs_tq1_5 * coef_Jacob1_qt_syms13 * pinvG3_1 +
                         coefs_tq2_5 * coef_Jacob1_qt_syms12 * pinvG2_2 +
                         coefs_tq3_5 * coef_Jacob1_qt_syms11 * pinvG1_3 +
                         coefs_tq1_2 * coef_Jacob1_qt_syms23 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob1_qt_syms22 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob1_qt_syms21 * pinvG1_3 +
                         coefs_tq2_5 * coef_Jacob1_qt_syms13 * pinvG3_2 +
                         coefs_tq3_5 * coef_Jacob1_qt_syms12 * pinvG2_3 +
                         coefs_tq2_2 * coef_Jacob1_qt_syms23 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob1_qt_syms22 * pinvG2_3 +
                         coefs_tq3_5 * coef_Jacob1_qt_syms13 * pinvG3_3 +
                         coefs_tq3_2 * coef_Jacob1_qt_syms23 * pinvG3_3;
    coef_f_q_sym(0, 5) = coef_Jacob1_qt_syms6 + coefs_tq1_6 * coef_Jacob1_qt_syms11 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob1_qt_syms21 * pinvG1_1 +
                         coefs_tq1_6 * coef_Jacob1_qt_syms12 * pinvG2_1 +
                         coefs_tq2_6 * coef_Jacob1_qt_syms11 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob1_qt_syms28 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob1_qt_syms22 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob1_qt_syms21 * pinvG1_2 +
                         coefs_tq1_6 * coef_Jacob1_qt_syms13 * pinvG3_1 +
                         coefs_tq2_6 * coef_Jacob1_qt_syms12 * pinvG2_2 +
                         coefs_tq3_6 * coef_Jacob1_qt_syms11 * pinvG1_3 +
                         coefs_tq1_2 * coef_Jacob1_qt_syms29 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob1_qt_syms28 * pinvG1_2 +
                         coefs_tq1_3 * coef_Jacob1_qt_syms23 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob1_qt_syms22 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob1_qt_syms21 * pinvG1_3 +
                         coefs_tq2_6 * coef_Jacob1_qt_syms13 * pinvG3_2 +
                         coefs_tq3_6 * coef_Jacob1_qt_syms12 * pinvG2_3 +
                         coefs_tq1_2 * coef_Jacob1_qt_syms30 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob1_qt_syms29 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob1_qt_syms28 * pinvG1_3 +
                         coefs_tq2_3 * coef_Jacob1_qt_syms23 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob1_qt_syms22 * pinvG2_3 +
                         coefs_tq3_6 * coef_Jacob1_qt_syms13 * pinvG3_3 +
                         coefs_tq2_2 * coef_Jacob1_qt_syms30 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob1_qt_syms29 * pinvG2_3 +
                         coefs_tq3_3 * coef_Jacob1_qt_syms23 * pinvG3_3 +
                         coefs_tq3_2 * coef_Jacob1_qt_syms30 * pinvG3_3;
    coef_f_q_sym(0, 6) = coef_Jacob1_qt_syms7 + coefs_tq1_7 * coef_Jacob1_qt_syms11 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob1_qt_syms21 * pinvG1_1 +
                         coefs_tq1_7 * coef_Jacob1_qt_syms12 * pinvG2_1 +
                         coefs_tq2_7 * coef_Jacob1_qt_syms11 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob1_qt_syms33 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob1_qt_syms22 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob1_qt_syms21 * pinvG1_2 +
                         coefs_tq1_7 * coef_Jacob1_qt_syms13 * pinvG3_1 +
                         coefs_tq2_7 * coef_Jacob1_qt_syms12 * pinvG2_2 +
                         coefs_tq3_7 * coef_Jacob1_qt_syms11 * pinvG1_3 +
                         coefs_tq1_2 * coef_Jacob1_qt_syms34 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob1_qt_syms33 * pinvG1_2 +
                         coefs_tq1_4 * coef_Jacob1_qt_syms23 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob1_qt_syms22 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob1_qt_syms21 * pinvG1_3 +
                         coefs_tq2_7 * coef_Jacob1_qt_syms13 * pinvG3_2 +
                         coefs_tq3_7 * coef_Jacob1_qt_syms12 * pinvG2_3 +
                         coefs_tq1_2 * coef_Jacob1_qt_syms35 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob1_qt_syms34 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob1_qt_syms33 * pinvG1_3 +
                         coefs_tq2_4 * coef_Jacob1_qt_syms23 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob1_qt_syms22 * pinvG2_3 +
                         coefs_tq3_7 * coef_Jacob1_qt_syms13 * pinvG3_3 +
                         coefs_tq2_2 * coef_Jacob1_qt_syms35 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob1_qt_syms34 * pinvG2_3 +
                         coefs_tq3_4 * coef_Jacob1_qt_syms23 * pinvG3_3 +
                         coefs_tq3_2 * coef_Jacob1_qt_syms35 * pinvG3_3;
    coef_f_q_sym(0, 7) = coef_Jacob1_qt_syms8 + coefs_tq1_8 * coef_Jacob1_qt_syms11 * pinvG1_1 +
                         coefs_tq1_8 * coef_Jacob1_qt_syms12 * pinvG2_1 +
                         coefs_tq2_8 * coef_Jacob1_qt_syms11 * pinvG1_2 +
                         coefs_tq1_3 * coef_Jacob1_qt_syms28 * pinvG1_1 +
                         coefs_tq1_8 * coef_Jacob1_qt_syms13 * pinvG3_1 +
                         coefs_tq2_8 * coef_Jacob1_qt_syms12 * pinvG2_2 +
                         coefs_tq3_8 * coef_Jacob1_qt_syms11 * pinvG1_3 +
                         coefs_tq1_3 * coef_Jacob1_qt_syms29 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob1_qt_syms28 * pinvG1_2 +
                         coefs_tq2_8 * coef_Jacob1_qt_syms13 * pinvG3_2 +
                         coefs_tq3_8 * coef_Jacob1_qt_syms12 * pinvG2_3 +
                         coefs_tq1_3 * coef_Jacob1_qt_syms30 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob1_qt_syms29 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob1_qt_syms28 * pinvG1_3 +
                         coefs_tq3_8 * coef_Jacob1_qt_syms13 * pinvG3_3 +
                         coefs_tq2_3 * coef_Jacob1_qt_syms30 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob1_qt_syms29 * pinvG2_3 +
                         coefs_tq3_3 * coef_Jacob1_qt_syms30 * pinvG3_3;
    coef_f_q_sym(0, 8) = coef_Jacob1_qt_syms9 + coefs_tq1_9 * coef_Jacob1_qt_syms11 * pinvG1_1 +
                         coefs_tq1_9 * coef_Jacob1_qt_syms12 * pinvG2_1 +
                         coefs_tq2_9 * coef_Jacob1_qt_syms11 * pinvG1_2 +
                         coefs_tq1_4 * coef_Jacob1_qt_syms28 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob1_qt_syms33 * pinvG1_1 +
                         coefs_tq1_9 * coef_Jacob1_qt_syms13 * pinvG3_1 +
                         coefs_tq2_9 * coef_Jacob1_qt_syms12 * pinvG2_2 +
                         coefs_tq3_9 * coef_Jacob1_qt_syms11 * pinvG1_3 +
                         coefs_tq1_4 * coef_Jacob1_qt_syms29 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob1_qt_syms28 * pinvG1_2 +
                         coefs_tq1_3 * coef_Jacob1_qt_syms34 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob1_qt_syms33 * pinvG1_2 +
                         coefs_tq2_9 * coef_Jacob1_qt_syms13 * pinvG3_2 +
                         coefs_tq3_9 * coef_Jacob1_qt_syms12 * pinvG2_3 +
                         coefs_tq1_4 * coef_Jacob1_qt_syms30 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob1_qt_syms29 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob1_qt_syms28 * pinvG1_3 +
                         coefs_tq1_3 * coef_Jacob1_qt_syms35 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob1_qt_syms34 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob1_qt_syms33 * pinvG1_3 +
                         coefs_tq3_9 * coef_Jacob1_qt_syms13 * pinvG3_3 +
                         coefs_tq2_4 * coef_Jacob1_qt_syms30 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob1_qt_syms29 * pinvG2_3 +
                         coefs_tq2_3 * coef_Jacob1_qt_syms35 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob1_qt_syms34 * pinvG2_3 +
                         coefs_tq3_4 * coef_Jacob1_qt_syms30 * pinvG3_3 +
                         coefs_tq3_3 * coef_Jacob1_qt_syms35 * pinvG3_3;
    coef_f_q_sym(0, 9) = coef_Jacob1_qt_syms10 + coefs_tq1_4 * coef_Jacob1_qt_syms33 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob1_qt_syms34 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob1_qt_syms33 * pinvG1_2 +
                         coefs_tq1_4 * coef_Jacob1_qt_syms35 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob1_qt_syms34 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob1_qt_syms33 * pinvG1_3 +
                         coefs_tq2_4 * coef_Jacob1_qt_syms35 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob1_qt_syms34 * pinvG2_3 +
                         coefs_tq3_4 * coef_Jacob1_qt_syms35 * pinvG3_3 +
                         coefs_tq1_10 * coef_Jacob1_qt_syms11 * pinvG1_1 +
                         coefs_tq1_10 * coef_Jacob1_qt_syms12 * pinvG2_1 +
                         coefs_tq1_10 * coef_Jacob1_qt_syms13 * pinvG3_1 +
                         coefs_tq2_10 * coef_Jacob1_qt_syms11 * pinvG1_2 +
                         coefs_tq2_10 * coef_Jacob1_qt_syms12 * pinvG2_2 +
                         coefs_tq2_10 * coef_Jacob1_qt_syms13 * pinvG3_2 +
                         coefs_tq3_10 * coef_Jacob1_qt_syms11 * pinvG1_3 +
                         coefs_tq3_10 * coef_Jacob1_qt_syms12 * pinvG2_3 +
                         coefs_tq3_10 * coef_Jacob1_qt_syms13 * pinvG3_3;
    coef_f_q_sym(0, 10) = coef_Jacob1_qt_syms14 + coefs_tq1_11 * coef_Jacob1_qt_syms11 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob1_qt_syms12 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob1_qt_syms13 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob1_qt_syms11 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob1_qt_syms12 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob1_qt_syms13 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob1_qt_syms11 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob1_qt_syms12 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob1_qt_syms13 * pinvG3_3;
    coef_f_q_sym(0, 11) = coef_Jacob1_qt_syms15 + coefs_tq1_5 * coef_Jacob1_qt_syms21 * pinvG1_1 +
                          coefs_tq1_5 * coef_Jacob1_qt_syms22 * pinvG2_1 +
                          coefs_tq2_5 * coef_Jacob1_qt_syms21 * pinvG1_2 +
                          coefs_tq1_5 * coef_Jacob1_qt_syms23 * pinvG3_1 +
                          coefs_tq2_5 * coef_Jacob1_qt_syms22 * pinvG2_2 +
                          coefs_tq3_5 * coef_Jacob1_qt_syms21 * pinvG1_3 +
                          coefs_tq2_5 * coef_Jacob1_qt_syms23 * pinvG3_2 +
                          coefs_tq3_5 * coef_Jacob1_qt_syms22 * pinvG2_3 +
                          coefs_tq3_5 * coef_Jacob1_qt_syms23 * pinvG3_3;
    coef_f_q_sym(0, 12) = coef_Jacob1_qt_syms16 + coefs_tq1_6 * coef_Jacob1_qt_syms21 * pinvG1_1 +
                          coefs_tq1_5 * coef_Jacob1_qt_syms28 * pinvG1_1 +
                          coefs_tq1_6 * coef_Jacob1_qt_syms22 * pinvG2_1 +
                          coefs_tq2_6 * coef_Jacob1_qt_syms21 * pinvG1_2 +
                          coefs_tq1_5 * coef_Jacob1_qt_syms29 * pinvG2_1 +
                          coefs_tq2_5 * coef_Jacob1_qt_syms28 * pinvG1_2 +
                          coefs_tq1_6 * coef_Jacob1_qt_syms23 * pinvG3_1 +
                          coefs_tq2_6 * coef_Jacob1_qt_syms22 * pinvG2_2 +
                          coefs_tq3_6 * coef_Jacob1_qt_syms21 * pinvG1_3 +
                          coefs_tq1_5 * coef_Jacob1_qt_syms30 * pinvG3_1 +
                          coefs_tq2_5 * coef_Jacob1_qt_syms29 * pinvG2_2 +
                          coefs_tq3_5 * coef_Jacob1_qt_syms28 * pinvG1_3 +
                          coefs_tq2_6 * coef_Jacob1_qt_syms23 * pinvG3_2 +
                          coefs_tq3_6 * coef_Jacob1_qt_syms22 * pinvG2_3 +
                          coefs_tq2_5 * coef_Jacob1_qt_syms30 * pinvG3_2 +
                          coefs_tq3_5 * coef_Jacob1_qt_syms29 * pinvG2_3 +
                          coefs_tq3_6 * coef_Jacob1_qt_syms23 * pinvG3_3 +
                          coefs_tq3_5 * coef_Jacob1_qt_syms30 * pinvG3_3;
    coef_f_q_sym(0, 13) = coef_Jacob1_qt_syms17 + coefs_tq1_7 * coef_Jacob1_qt_syms21 * pinvG1_1 +
                          coefs_tq1_5 * coef_Jacob1_qt_syms33 * pinvG1_1 +
                          coefs_tq1_7 * coef_Jacob1_qt_syms22 * pinvG2_1 +
                          coefs_tq2_7 * coef_Jacob1_qt_syms21 * pinvG1_2 +
                          coefs_tq1_5 * coef_Jacob1_qt_syms34 * pinvG2_1 +
                          coefs_tq2_5 * coef_Jacob1_qt_syms33 * pinvG1_2 +
                          coefs_tq1_7 * coef_Jacob1_qt_syms23 * pinvG3_1 +
                          coefs_tq2_7 * coef_Jacob1_qt_syms22 * pinvG2_2 +
                          coefs_tq3_7 * coef_Jacob1_qt_syms21 * pinvG1_3 +
                          coefs_tq1_5 * coef_Jacob1_qt_syms35 * pinvG3_1 +
                          coefs_tq2_5 * coef_Jacob1_qt_syms34 * pinvG2_2 +
                          coefs_tq3_5 * coef_Jacob1_qt_syms33 * pinvG1_3 +
                          coefs_tq2_7 * coef_Jacob1_qt_syms23 * pinvG3_2 +
                          coefs_tq3_7 * coef_Jacob1_qt_syms22 * pinvG2_3 +
                          coefs_tq2_5 * coef_Jacob1_qt_syms35 * pinvG3_2 +
                          coefs_tq3_5 * coef_Jacob1_qt_syms34 * pinvG2_3 +
                          coefs_tq3_7 * coef_Jacob1_qt_syms23 * pinvG3_3 +
                          coefs_tq3_5 * coef_Jacob1_qt_syms35 * pinvG3_3;
    coef_f_q_sym(0, 14) = coef_Jacob1_qt_syms18 + coefs_tq1_8 * coef_Jacob1_qt_syms21 * pinvG1_1 +
                          coefs_tq1_6 * coef_Jacob1_qt_syms28 * pinvG1_1 +
                          coefs_tq1_8 * coef_Jacob1_qt_syms22 * pinvG2_1 +
                          coefs_tq2_8 * coef_Jacob1_qt_syms21 * pinvG1_2 +
                          coefs_tq1_6 * coef_Jacob1_qt_syms29 * pinvG2_1 +
                          coefs_tq2_6 * coef_Jacob1_qt_syms28 * pinvG1_2 +
                          coefs_tq1_8 * coef_Jacob1_qt_syms23 * pinvG3_1 +
                          coefs_tq2_8 * coef_Jacob1_qt_syms22 * pinvG2_2 +
                          coefs_tq3_8 * coef_Jacob1_qt_syms21 * pinvG1_3 +
                          coefs_tq1_6 * coef_Jacob1_qt_syms30 * pinvG3_1 +
                          coefs_tq2_6 * coef_Jacob1_qt_syms29 * pinvG2_2 +
                          coefs_tq3_6 * coef_Jacob1_qt_syms28 * pinvG1_3 +
                          coefs_tq2_8 * coef_Jacob1_qt_syms23 * pinvG3_2 +
                          coefs_tq3_8 * coef_Jacob1_qt_syms22 * pinvG2_3 +
                          coefs_tq2_6 * coef_Jacob1_qt_syms30 * pinvG3_2 +
                          coefs_tq3_6 * coef_Jacob1_qt_syms29 * pinvG2_3 +
                          coefs_tq3_8 * coef_Jacob1_qt_syms23 * pinvG3_3 +
                          coefs_tq3_6 * coef_Jacob1_qt_syms30 * pinvG3_3;
    coef_f_q_sym(0, 15) = coef_Jacob1_qt_syms19 + coefs_tq1_9 * coef_Jacob1_qt_syms21 * pinvG1_1 +
                          coefs_tq1_7 * coef_Jacob1_qt_syms28 * pinvG1_1 +
                          coefs_tq1_6 * coef_Jacob1_qt_syms33 * pinvG1_1 +
                          coefs_tq1_9 * coef_Jacob1_qt_syms22 * pinvG2_1 +
                          coefs_tq2_9 * coef_Jacob1_qt_syms21 * pinvG1_2 +
                          coefs_tq1_7 * coef_Jacob1_qt_syms29 * pinvG2_1 +
                          coefs_tq2_7 * coef_Jacob1_qt_syms28 * pinvG1_2 +
                          coefs_tq1_6 * coef_Jacob1_qt_syms34 * pinvG2_1 +
                          coefs_tq2_6 * coef_Jacob1_qt_syms33 * pinvG1_2 +
                          coefs_tq1_9 * coef_Jacob1_qt_syms23 * pinvG3_1 +
                          coefs_tq2_9 * coef_Jacob1_qt_syms22 * pinvG2_2 +
                          coefs_tq3_9 * coef_Jacob1_qt_syms21 * pinvG1_3 +
                          coefs_tq1_7 * coef_Jacob1_qt_syms30 * pinvG3_1 +
                          coefs_tq2_7 * coef_Jacob1_qt_syms29 * pinvG2_2 +
                          coefs_tq3_7 * coef_Jacob1_qt_syms28 * pinvG1_3 +
                          coefs_tq1_6 * coef_Jacob1_qt_syms35 * pinvG3_1 +
                          coefs_tq2_6 * coef_Jacob1_qt_syms34 * pinvG2_2 +
                          coefs_tq3_6 * coef_Jacob1_qt_syms33 * pinvG1_3 +
                          coefs_tq2_9 * coef_Jacob1_qt_syms23 * pinvG3_2 +
                          coefs_tq3_9 * coef_Jacob1_qt_syms22 * pinvG2_3 +
                          coefs_tq2_7 * coef_Jacob1_qt_syms30 * pinvG3_2 +
                          coefs_tq3_7 * coef_Jacob1_qt_syms29 * pinvG2_3 +
                          coefs_tq2_6 * coef_Jacob1_qt_syms35 * pinvG3_2 +
                          coefs_tq3_6 * coef_Jacob1_qt_syms34 * pinvG2_3 +
                          coefs_tq3_9 * coef_Jacob1_qt_syms23 * pinvG3_3 +
                          coefs_tq3_7 * coef_Jacob1_qt_syms30 * pinvG3_3 +
                          coefs_tq3_6 * coef_Jacob1_qt_syms35 * pinvG3_3;
    coef_f_q_sym(0, 16) = coef_Jacob1_qt_syms20 + coefs_tq1_7 * coef_Jacob1_qt_syms33 * pinvG1_1 +
                          coefs_tq1_7 * coef_Jacob1_qt_syms34 * pinvG2_1 +
                          coefs_tq2_7 * coef_Jacob1_qt_syms33 * pinvG1_2 +
                          coefs_tq1_7 * coef_Jacob1_qt_syms35 * pinvG3_1 +
                          coefs_tq2_7 * coef_Jacob1_qt_syms34 * pinvG2_2 +
                          coefs_tq3_7 * coef_Jacob1_qt_syms33 * pinvG1_3 +
                          coefs_tq2_7 * coef_Jacob1_qt_syms35 * pinvG3_2 +
                          coefs_tq3_7 * coef_Jacob1_qt_syms34 * pinvG2_3 +
                          coefs_tq3_7 * coef_Jacob1_qt_syms35 * pinvG3_3 +
                          coefs_tq1_10 * coef_Jacob1_qt_syms21 * pinvG1_1 +
                          coefs_tq1_10 * coef_Jacob1_qt_syms22 * pinvG2_1 +
                          coefs_tq1_10 * coef_Jacob1_qt_syms23 * pinvG3_1 +
                          coefs_tq2_10 * coef_Jacob1_qt_syms21 * pinvG1_2 +
                          coefs_tq2_10 * coef_Jacob1_qt_syms22 * pinvG2_2 +
                          coefs_tq2_10 * coef_Jacob1_qt_syms23 * pinvG3_2 +
                          coefs_tq3_10 * coef_Jacob1_qt_syms21 * pinvG1_3 +
                          coefs_tq3_10 * coef_Jacob1_qt_syms22 * pinvG2_3 +
                          coefs_tq3_10 * coef_Jacob1_qt_syms23 * pinvG3_3;
    coef_f_q_sym(0, 17) = coef_Jacob1_qt_syms24 + coefs_tq1_11 * coef_Jacob1_qt_syms21 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob1_qt_syms22 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob1_qt_syms23 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob1_qt_syms21 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob1_qt_syms22 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob1_qt_syms23 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob1_qt_syms21 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob1_qt_syms22 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob1_qt_syms23 * pinvG3_3;
    coef_f_q_sym(0, 18) = coef_Jacob1_qt_syms25 + coefs_tq1_8 * coef_Jacob1_qt_syms28 * pinvG1_1 +
                          coefs_tq1_8 * coef_Jacob1_qt_syms29 * pinvG2_1 +
                          coefs_tq2_8 * coef_Jacob1_qt_syms28 * pinvG1_2 +
                          coefs_tq1_8 * coef_Jacob1_qt_syms30 * pinvG3_1 +
                          coefs_tq2_8 * coef_Jacob1_qt_syms29 * pinvG2_2 +
                          coefs_tq3_8 * coef_Jacob1_qt_syms28 * pinvG1_3 +
                          coefs_tq2_8 * coef_Jacob1_qt_syms30 * pinvG3_2 +
                          coefs_tq3_8 * coef_Jacob1_qt_syms29 * pinvG2_3 +
                          coefs_tq3_8 * coef_Jacob1_qt_syms30 * pinvG3_3;
    coef_f_q_sym(0, 19) = coef_Jacob1_qt_syms26 + coefs_tq1_9 * coef_Jacob1_qt_syms28 * pinvG1_1 +
                          coefs_tq1_8 * coef_Jacob1_qt_syms33 * pinvG1_1 +
                          coefs_tq1_9 * coef_Jacob1_qt_syms29 * pinvG2_1 +
                          coefs_tq2_9 * coef_Jacob1_qt_syms28 * pinvG1_2 +
                          coefs_tq1_8 * coef_Jacob1_qt_syms34 * pinvG2_1 +
                          coefs_tq2_8 * coef_Jacob1_qt_syms33 * pinvG1_2 +
                          coefs_tq1_9 * coef_Jacob1_qt_syms30 * pinvG3_1 +
                          coefs_tq2_9 * coef_Jacob1_qt_syms29 * pinvG2_2 +
                          coefs_tq3_9 * coef_Jacob1_qt_syms28 * pinvG1_3 +
                          coefs_tq1_8 * coef_Jacob1_qt_syms35 * pinvG3_1 +
                          coefs_tq2_8 * coef_Jacob1_qt_syms34 * pinvG2_2 +
                          coefs_tq3_8 * coef_Jacob1_qt_syms33 * pinvG1_3 +
                          coefs_tq2_9 * coef_Jacob1_qt_syms30 * pinvG3_2 +
                          coefs_tq3_9 * coef_Jacob1_qt_syms29 * pinvG2_3 +
                          coefs_tq2_8 * coef_Jacob1_qt_syms35 * pinvG3_2 +
                          coefs_tq3_8 * coef_Jacob1_qt_syms34 * pinvG2_3 +
                          coefs_tq3_9 * coef_Jacob1_qt_syms30 * pinvG3_3 +
                          coefs_tq3_8 * coef_Jacob1_qt_syms35 * pinvG3_3;
    coef_f_q_sym(0, 20) = coef_Jacob1_qt_syms27 + coefs_tq1_9 * coef_Jacob1_qt_syms33 * pinvG1_1 +
                          coefs_tq1_9 * coef_Jacob1_qt_syms34 * pinvG2_1 +
                          coefs_tq2_9 * coef_Jacob1_qt_syms33 * pinvG1_2 +
                          coefs_tq1_9 * coef_Jacob1_qt_syms35 * pinvG3_1 +
                          coefs_tq2_9 * coef_Jacob1_qt_syms34 * pinvG2_2 +
                          coefs_tq3_9 * coef_Jacob1_qt_syms33 * pinvG1_3 +
                          coefs_tq2_9 * coef_Jacob1_qt_syms35 * pinvG3_2 +
                          coefs_tq3_9 * coef_Jacob1_qt_syms34 * pinvG2_3 +
                          coefs_tq3_9 * coef_Jacob1_qt_syms35 * pinvG3_3 +
                          coefs_tq1_10 * coef_Jacob1_qt_syms28 * pinvG1_1 +
                          coefs_tq1_10 * coef_Jacob1_qt_syms29 * pinvG2_1 +
                          coefs_tq1_10 * coef_Jacob1_qt_syms30 * pinvG3_1 +
                          coefs_tq2_10 * coef_Jacob1_qt_syms28 * pinvG1_2 +
                          coefs_tq2_10 * coef_Jacob1_qt_syms29 * pinvG2_2 +
                          coefs_tq2_10 * coef_Jacob1_qt_syms30 * pinvG3_2 +
                          coefs_tq3_10 * coef_Jacob1_qt_syms28 * pinvG1_3 +
                          coefs_tq3_10 * coef_Jacob1_qt_syms29 * pinvG2_3 +
                          coefs_tq3_10 * coef_Jacob1_qt_syms30 * pinvG3_3;
    coef_f_q_sym(0, 21) = coef_Jacob1_qt_syms31 + coefs_tq1_11 * coef_Jacob1_qt_syms28 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob1_qt_syms29 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob1_qt_syms30 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob1_qt_syms28 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob1_qt_syms29 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob1_qt_syms30 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob1_qt_syms28 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob1_qt_syms29 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob1_qt_syms30 * pinvG3_3;
    coef_f_q_sym(0, 22) = coef_Jacob1_qt_syms32 + coefs_tq1_10 * coef_Jacob1_qt_syms33 * pinvG1_1 +
                          coefs_tq1_10 * coef_Jacob1_qt_syms34 * pinvG2_1 +
                          coefs_tq1_10 * coef_Jacob1_qt_syms35 * pinvG3_1 +
                          coefs_tq2_10 * coef_Jacob1_qt_syms33 * pinvG1_2 +
                          coefs_tq2_10 * coef_Jacob1_qt_syms34 * pinvG2_2 +
                          coefs_tq2_10 * coef_Jacob1_qt_syms35 * pinvG3_2 +
                          coefs_tq3_10 * coef_Jacob1_qt_syms33 * pinvG1_3 +
                          coefs_tq3_10 * coef_Jacob1_qt_syms34 * pinvG2_3 +
                          coefs_tq3_10 * coef_Jacob1_qt_syms35 * pinvG3_3;
    coef_f_q_sym(0, 23) = coef_Jacob1_qt_syms36 + coefs_tq1_11 * coef_Jacob1_qt_syms33 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob1_qt_syms34 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob1_qt_syms35 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob1_qt_syms33 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob1_qt_syms34 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob1_qt_syms35 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob1_qt_syms33 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob1_qt_syms34 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob1_qt_syms35 * pinvG3_3;
    coef_f_q_sym(1, 0) = coef_Jacob2_qt_syms1 + coefs_tq1_1 * coef_Jacob2_qt_syms11 * pinvG1_1 +
                         coefs_tq1_1 * coef_Jacob2_qt_syms12 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob2_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob2_qt_syms13 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob2_qt_syms12 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob2_qt_syms11 * pinvG1_3 +
                         coefs_tq2_1 * coef_Jacob2_qt_syms13 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob2_qt_syms12 * pinvG2_3 +
                         coefs_tq3_1 * coef_Jacob2_qt_syms13 * pinvG3_3;
    coef_f_q_sym(1, 1) = coef_Jacob2_qt_syms2 + coefs_tq1_2 * coef_Jacob2_qt_syms11 * pinvG1_1 +
                         coefs_tq1_1 * coef_Jacob2_qt_syms21 * pinvG1_1 +
                         coefs_tq1_2 * coef_Jacob2_qt_syms12 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob2_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob2_qt_syms22 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob2_qt_syms21 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob2_qt_syms13 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob2_qt_syms12 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob2_qt_syms11 * pinvG1_3 +
                         coefs_tq1_1 * coef_Jacob2_qt_syms23 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob2_qt_syms22 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob2_qt_syms21 * pinvG1_3 +
                         coefs_tq2_2 * coef_Jacob2_qt_syms13 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob2_qt_syms12 * pinvG2_3 +
                         coefs_tq2_1 * coef_Jacob2_qt_syms23 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob2_qt_syms22 * pinvG2_3 +
                         coefs_tq3_2 * coef_Jacob2_qt_syms13 * pinvG3_3 +
                         coefs_tq3_1 * coef_Jacob2_qt_syms23 * pinvG3_3;
    coef_f_q_sym(1, 2) = coef_Jacob2_qt_syms3 + coefs_tq1_3 * coef_Jacob2_qt_syms11 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob2_qt_syms12 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob2_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob2_qt_syms28 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob2_qt_syms13 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob2_qt_syms12 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob2_qt_syms11 * pinvG1_3 +
                         coefs_tq1_1 * coef_Jacob2_qt_syms29 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob2_qt_syms28 * pinvG1_2 +
                         coefs_tq2_3 * coef_Jacob2_qt_syms13 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob2_qt_syms12 * pinvG2_3 +
                         coefs_tq1_1 * coef_Jacob2_qt_syms30 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob2_qt_syms29 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob2_qt_syms28 * pinvG1_3 +
                         coefs_tq3_3 * coef_Jacob2_qt_syms13 * pinvG3_3 +
                         coefs_tq2_1 * coef_Jacob2_qt_syms30 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob2_qt_syms29 * pinvG2_3 +
                         coefs_tq3_1 * coef_Jacob2_qt_syms30 * pinvG3_3;
    coef_f_q_sym(1, 3) = coef_Jacob2_qt_syms4 + coefs_tq1_4 * coef_Jacob2_qt_syms11 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob2_qt_syms12 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob2_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob2_qt_syms33 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob2_qt_syms13 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob2_qt_syms12 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob2_qt_syms11 * pinvG1_3 +
                         coefs_tq1_1 * coef_Jacob2_qt_syms34 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob2_qt_syms33 * pinvG1_2 +
                         coefs_tq2_4 * coef_Jacob2_qt_syms13 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob2_qt_syms12 * pinvG2_3 +
                         coefs_tq1_1 * coef_Jacob2_qt_syms35 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob2_qt_syms34 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob2_qt_syms33 * pinvG1_3 +
                         coefs_tq3_4 * coef_Jacob2_qt_syms13 * pinvG3_3 +
                         coefs_tq2_1 * coef_Jacob2_qt_syms35 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob2_qt_syms34 * pinvG2_3 +
                         coefs_tq3_1 * coef_Jacob2_qt_syms35 * pinvG3_3;
    coef_f_q_sym(1, 4) = coef_Jacob2_qt_syms5 + coefs_tq1_5 * coef_Jacob2_qt_syms11 * pinvG1_1 +
                         coefs_tq1_2 * coef_Jacob2_qt_syms21 * pinvG1_1 +
                         coefs_tq1_5 * coef_Jacob2_qt_syms12 * pinvG2_1 +
                         coefs_tq2_5 * coef_Jacob2_qt_syms11 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob2_qt_syms22 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob2_qt_syms21 * pinvG1_2 +
                         coefs_tq1_5 * coef_Jacob2_qt_syms13 * pinvG3_1 +
                         coefs_tq2_5 * coef_Jacob2_qt_syms12 * pinvG2_2 +
                         coefs_tq3_5 * coef_Jacob2_qt_syms11 * pinvG1_3 +
                         coefs_tq1_2 * coef_Jacob2_qt_syms23 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob2_qt_syms22 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob2_qt_syms21 * pinvG1_3 +
                         coefs_tq2_5 * coef_Jacob2_qt_syms13 * pinvG3_2 +
                         coefs_tq3_5 * coef_Jacob2_qt_syms12 * pinvG2_3 +
                         coefs_tq2_2 * coef_Jacob2_qt_syms23 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob2_qt_syms22 * pinvG2_3 +
                         coefs_tq3_5 * coef_Jacob2_qt_syms13 * pinvG3_3 +
                         coefs_tq3_2 * coef_Jacob2_qt_syms23 * pinvG3_3;
    coef_f_q_sym(1, 5) = coef_Jacob2_qt_syms6 + coefs_tq1_6 * coef_Jacob2_qt_syms11 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob2_qt_syms21 * pinvG1_1 +
                         coefs_tq1_6 * coef_Jacob2_qt_syms12 * pinvG2_1 +
                         coefs_tq2_6 * coef_Jacob2_qt_syms11 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob2_qt_syms28 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob2_qt_syms22 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob2_qt_syms21 * pinvG1_2 +
                         coefs_tq1_6 * coef_Jacob2_qt_syms13 * pinvG3_1 +
                         coefs_tq2_6 * coef_Jacob2_qt_syms12 * pinvG2_2 +
                         coefs_tq3_6 * coef_Jacob2_qt_syms11 * pinvG1_3 +
                         coefs_tq1_2 * coef_Jacob2_qt_syms29 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob2_qt_syms28 * pinvG1_2 +
                         coefs_tq1_3 * coef_Jacob2_qt_syms23 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob2_qt_syms22 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob2_qt_syms21 * pinvG1_3 +
                         coefs_tq2_6 * coef_Jacob2_qt_syms13 * pinvG3_2 +
                         coefs_tq3_6 * coef_Jacob2_qt_syms12 * pinvG2_3 +
                         coefs_tq1_2 * coef_Jacob2_qt_syms30 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob2_qt_syms29 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob2_qt_syms28 * pinvG1_3 +
                         coefs_tq2_3 * coef_Jacob2_qt_syms23 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob2_qt_syms22 * pinvG2_3 +
                         coefs_tq3_6 * coef_Jacob2_qt_syms13 * pinvG3_3 +
                         coefs_tq2_2 * coef_Jacob2_qt_syms30 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob2_qt_syms29 * pinvG2_3 +
                         coefs_tq3_3 * coef_Jacob2_qt_syms23 * pinvG3_3 +
                         coefs_tq3_2 * coef_Jacob2_qt_syms30 * pinvG3_3;
    coef_f_q_sym(1, 6) = coef_Jacob2_qt_syms7 + coefs_tq1_7 * coef_Jacob2_qt_syms11 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob2_qt_syms21 * pinvG1_1 +
                         coefs_tq1_7 * coef_Jacob2_qt_syms12 * pinvG2_1 +
                         coefs_tq2_7 * coef_Jacob2_qt_syms11 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob2_qt_syms33 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob2_qt_syms22 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob2_qt_syms21 * pinvG1_2 +
                         coefs_tq1_7 * coef_Jacob2_qt_syms13 * pinvG3_1 +
                         coefs_tq2_7 * coef_Jacob2_qt_syms12 * pinvG2_2 +
                         coefs_tq3_7 * coef_Jacob2_qt_syms11 * pinvG1_3 +
                         coefs_tq1_2 * coef_Jacob2_qt_syms34 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob2_qt_syms33 * pinvG1_2 +
                         coefs_tq1_4 * coef_Jacob2_qt_syms23 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob2_qt_syms22 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob2_qt_syms21 * pinvG1_3 +
                         coefs_tq2_7 * coef_Jacob2_qt_syms13 * pinvG3_2 +
                         coefs_tq3_7 * coef_Jacob2_qt_syms12 * pinvG2_3 +
                         coefs_tq1_2 * coef_Jacob2_qt_syms35 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob2_qt_syms34 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob2_qt_syms33 * pinvG1_3 +
                         coefs_tq2_4 * coef_Jacob2_qt_syms23 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob2_qt_syms22 * pinvG2_3 +
                         coefs_tq3_7 * coef_Jacob2_qt_syms13 * pinvG3_3 +
                         coefs_tq2_2 * coef_Jacob2_qt_syms35 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob2_qt_syms34 * pinvG2_3 +
                         coefs_tq3_4 * coef_Jacob2_qt_syms23 * pinvG3_3 +
                         coefs_tq3_2 * coef_Jacob2_qt_syms35 * pinvG3_3;
    coef_f_q_sym(1, 7) = coef_Jacob2_qt_syms8 + coefs_tq1_8 * coef_Jacob2_qt_syms11 * pinvG1_1 +
                         coefs_tq1_8 * coef_Jacob2_qt_syms12 * pinvG2_1 +
                         coefs_tq2_8 * coef_Jacob2_qt_syms11 * pinvG1_2 +
                         coefs_tq1_3 * coef_Jacob2_qt_syms28 * pinvG1_1 +
                         coefs_tq1_8 * coef_Jacob2_qt_syms13 * pinvG3_1 +
                         coefs_tq2_8 * coef_Jacob2_qt_syms12 * pinvG2_2 +
                         coefs_tq3_8 * coef_Jacob2_qt_syms11 * pinvG1_3 +
                         coefs_tq1_3 * coef_Jacob2_qt_syms29 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob2_qt_syms28 * pinvG1_2 +
                         coefs_tq2_8 * coef_Jacob2_qt_syms13 * pinvG3_2 +
                         coefs_tq3_8 * coef_Jacob2_qt_syms12 * pinvG2_3 +
                         coefs_tq1_3 * coef_Jacob2_qt_syms30 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob2_qt_syms29 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob2_qt_syms28 * pinvG1_3 +
                         coefs_tq3_8 * coef_Jacob2_qt_syms13 * pinvG3_3 +
                         coefs_tq2_3 * coef_Jacob2_qt_syms30 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob2_qt_syms29 * pinvG2_3 +
                         coefs_tq3_3 * coef_Jacob2_qt_syms30 * pinvG3_3;
    coef_f_q_sym(1, 8) = coef_Jacob2_qt_syms9 + coefs_tq1_9 * coef_Jacob2_qt_syms11 * pinvG1_1 +
                         coefs_tq1_9 * coef_Jacob2_qt_syms12 * pinvG2_1 +
                         coefs_tq2_9 * coef_Jacob2_qt_syms11 * pinvG1_2 +
                         coefs_tq1_4 * coef_Jacob2_qt_syms28 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob2_qt_syms33 * pinvG1_1 +
                         coefs_tq1_9 * coef_Jacob2_qt_syms13 * pinvG3_1 +
                         coefs_tq2_9 * coef_Jacob2_qt_syms12 * pinvG2_2 +
                         coefs_tq3_9 * coef_Jacob2_qt_syms11 * pinvG1_3 +
                         coefs_tq1_4 * coef_Jacob2_qt_syms29 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob2_qt_syms28 * pinvG1_2 +
                         coefs_tq1_3 * coef_Jacob2_qt_syms34 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob2_qt_syms33 * pinvG1_2 +
                         coefs_tq2_9 * coef_Jacob2_qt_syms13 * pinvG3_2 +
                         coefs_tq3_9 * coef_Jacob2_qt_syms12 * pinvG2_3 +
                         coefs_tq1_4 * coef_Jacob2_qt_syms30 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob2_qt_syms29 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob2_qt_syms28 * pinvG1_3 +
                         coefs_tq1_3 * coef_Jacob2_qt_syms35 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob2_qt_syms34 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob2_qt_syms33 * pinvG1_3 +
                         coefs_tq3_9 * coef_Jacob2_qt_syms13 * pinvG3_3 +
                         coefs_tq2_4 * coef_Jacob2_qt_syms30 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob2_qt_syms29 * pinvG2_3 +
                         coefs_tq2_3 * coef_Jacob2_qt_syms35 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob2_qt_syms34 * pinvG2_3 +
                         coefs_tq3_4 * coef_Jacob2_qt_syms30 * pinvG3_3 +
                         coefs_tq3_3 * coef_Jacob2_qt_syms35 * pinvG3_3;
    coef_f_q_sym(1, 9) = coef_Jacob2_qt_syms10 + coefs_tq1_4 * coef_Jacob2_qt_syms33 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob2_qt_syms34 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob2_qt_syms33 * pinvG1_2 +
                         coefs_tq1_4 * coef_Jacob2_qt_syms35 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob2_qt_syms34 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob2_qt_syms33 * pinvG1_3 +
                         coefs_tq2_4 * coef_Jacob2_qt_syms35 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob2_qt_syms34 * pinvG2_3 +
                         coefs_tq3_4 * coef_Jacob2_qt_syms35 * pinvG3_3 +
                         coefs_tq1_10 * coef_Jacob2_qt_syms11 * pinvG1_1 +
                         coefs_tq1_10 * coef_Jacob2_qt_syms12 * pinvG2_1 +
                         coefs_tq1_10 * coef_Jacob2_qt_syms13 * pinvG3_1 +
                         coefs_tq2_10 * coef_Jacob2_qt_syms11 * pinvG1_2 +
                         coefs_tq2_10 * coef_Jacob2_qt_syms12 * pinvG2_2 +
                         coefs_tq2_10 * coef_Jacob2_qt_syms13 * pinvG3_2 +
                         coefs_tq3_10 * coef_Jacob2_qt_syms11 * pinvG1_3 +
                         coefs_tq3_10 * coef_Jacob2_qt_syms12 * pinvG2_3 +
                         coefs_tq3_10 * coef_Jacob2_qt_syms13 * pinvG3_3;
    coef_f_q_sym(1, 10) = coef_Jacob2_qt_syms14 + coefs_tq1_11 * coef_Jacob2_qt_syms11 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob2_qt_syms12 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob2_qt_syms13 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob2_qt_syms11 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob2_qt_syms12 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob2_qt_syms13 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob2_qt_syms11 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob2_qt_syms12 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob2_qt_syms13 * pinvG3_3;
    coef_f_q_sym(1, 11) = coef_Jacob2_qt_syms15 + coefs_tq1_5 * coef_Jacob2_qt_syms21 * pinvG1_1 +
                          coefs_tq1_5 * coef_Jacob2_qt_syms22 * pinvG2_1 +
                          coefs_tq2_5 * coef_Jacob2_qt_syms21 * pinvG1_2 +
                          coefs_tq1_5 * coef_Jacob2_qt_syms23 * pinvG3_1 +
                          coefs_tq2_5 * coef_Jacob2_qt_syms22 * pinvG2_2 +
                          coefs_tq3_5 * coef_Jacob2_qt_syms21 * pinvG1_3 +
                          coefs_tq2_5 * coef_Jacob2_qt_syms23 * pinvG3_2 +
                          coefs_tq3_5 * coef_Jacob2_qt_syms22 * pinvG2_3 +
                          coefs_tq3_5 * coef_Jacob2_qt_syms23 * pinvG3_3;
    coef_f_q_sym(1, 12) = coef_Jacob2_qt_syms16 + coefs_tq1_6 * coef_Jacob2_qt_syms21 * pinvG1_1 +
                          coefs_tq1_5 * coef_Jacob2_qt_syms28 * pinvG1_1 +
                          coefs_tq1_6 * coef_Jacob2_qt_syms22 * pinvG2_1 +
                          coefs_tq2_6 * coef_Jacob2_qt_syms21 * pinvG1_2 +
                          coefs_tq1_5 * coef_Jacob2_qt_syms29 * pinvG2_1 +
                          coefs_tq2_5 * coef_Jacob2_qt_syms28 * pinvG1_2 +
                          coefs_tq1_6 * coef_Jacob2_qt_syms23 * pinvG3_1 +
                          coefs_tq2_6 * coef_Jacob2_qt_syms22 * pinvG2_2 +
                          coefs_tq3_6 * coef_Jacob2_qt_syms21 * pinvG1_3 +
                          coefs_tq1_5 * coef_Jacob2_qt_syms30 * pinvG3_1 +
                          coefs_tq2_5 * coef_Jacob2_qt_syms29 * pinvG2_2 +
                          coefs_tq3_5 * coef_Jacob2_qt_syms28 * pinvG1_3 +
                          coefs_tq2_6 * coef_Jacob2_qt_syms23 * pinvG3_2 +
                          coefs_tq3_6 * coef_Jacob2_qt_syms22 * pinvG2_3 +
                          coefs_tq2_5 * coef_Jacob2_qt_syms30 * pinvG3_2 +
                          coefs_tq3_5 * coef_Jacob2_qt_syms29 * pinvG2_3 +
                          coefs_tq3_6 * coef_Jacob2_qt_syms23 * pinvG3_3 +
                          coefs_tq3_5 * coef_Jacob2_qt_syms30 * pinvG3_3;
    coef_f_q_sym(1, 13) = coef_Jacob2_qt_syms17 + coefs_tq1_7 * coef_Jacob2_qt_syms21 * pinvG1_1 +
                          coefs_tq1_5 * coef_Jacob2_qt_syms33 * pinvG1_1 +
                          coefs_tq1_7 * coef_Jacob2_qt_syms22 * pinvG2_1 +
                          coefs_tq2_7 * coef_Jacob2_qt_syms21 * pinvG1_2 +
                          coefs_tq1_5 * coef_Jacob2_qt_syms34 * pinvG2_1 +
                          coefs_tq2_5 * coef_Jacob2_qt_syms33 * pinvG1_2 +
                          coefs_tq1_7 * coef_Jacob2_qt_syms23 * pinvG3_1 +
                          coefs_tq2_7 * coef_Jacob2_qt_syms22 * pinvG2_2 +
                          coefs_tq3_7 * coef_Jacob2_qt_syms21 * pinvG1_3 +
                          coefs_tq1_5 * coef_Jacob2_qt_syms35 * pinvG3_1 +
                          coefs_tq2_5 * coef_Jacob2_qt_syms34 * pinvG2_2 +
                          coefs_tq3_5 * coef_Jacob2_qt_syms33 * pinvG1_3 +
                          coefs_tq2_7 * coef_Jacob2_qt_syms23 * pinvG3_2 +
                          coefs_tq3_7 * coef_Jacob2_qt_syms22 * pinvG2_3 +
                          coefs_tq2_5 * coef_Jacob2_qt_syms35 * pinvG3_2 +
                          coefs_tq3_5 * coef_Jacob2_qt_syms34 * pinvG2_3 +
                          coefs_tq3_7 * coef_Jacob2_qt_syms23 * pinvG3_3 +
                          coefs_tq3_5 * coef_Jacob2_qt_syms35 * pinvG3_3;
    coef_f_q_sym(1, 14) = coef_Jacob2_qt_syms18 + coefs_tq1_8 * coef_Jacob2_qt_syms21 * pinvG1_1 +
                          coefs_tq1_6 * coef_Jacob2_qt_syms28 * pinvG1_1 +
                          coefs_tq1_8 * coef_Jacob2_qt_syms22 * pinvG2_1 +
                          coefs_tq2_8 * coef_Jacob2_qt_syms21 * pinvG1_2 +
                          coefs_tq1_6 * coef_Jacob2_qt_syms29 * pinvG2_1 +
                          coefs_tq2_6 * coef_Jacob2_qt_syms28 * pinvG1_2 +
                          coefs_tq1_8 * coef_Jacob2_qt_syms23 * pinvG3_1 +
                          coefs_tq2_8 * coef_Jacob2_qt_syms22 * pinvG2_2 +
                          coefs_tq3_8 * coef_Jacob2_qt_syms21 * pinvG1_3 +
                          coefs_tq1_6 * coef_Jacob2_qt_syms30 * pinvG3_1 +
                          coefs_tq2_6 * coef_Jacob2_qt_syms29 * pinvG2_2 +
                          coefs_tq3_6 * coef_Jacob2_qt_syms28 * pinvG1_3 +
                          coefs_tq2_8 * coef_Jacob2_qt_syms23 * pinvG3_2 +
                          coefs_tq3_8 * coef_Jacob2_qt_syms22 * pinvG2_3 +
                          coefs_tq2_6 * coef_Jacob2_qt_syms30 * pinvG3_2 +
                          coefs_tq3_6 * coef_Jacob2_qt_syms29 * pinvG2_3 +
                          coefs_tq3_8 * coef_Jacob2_qt_syms23 * pinvG3_3 +
                          coefs_tq3_6 * coef_Jacob2_qt_syms30 * pinvG3_3;
    coef_f_q_sym(1, 15) = coef_Jacob2_qt_syms19 + coefs_tq1_9 * coef_Jacob2_qt_syms21 * pinvG1_1 +
                          coefs_tq1_7 * coef_Jacob2_qt_syms28 * pinvG1_1 +
                          coefs_tq1_6 * coef_Jacob2_qt_syms33 * pinvG1_1 +
                          coefs_tq1_9 * coef_Jacob2_qt_syms22 * pinvG2_1 +
                          coefs_tq2_9 * coef_Jacob2_qt_syms21 * pinvG1_2 +
                          coefs_tq1_7 * coef_Jacob2_qt_syms29 * pinvG2_1 +
                          coefs_tq2_7 * coef_Jacob2_qt_syms28 * pinvG1_2 +
                          coefs_tq1_6 * coef_Jacob2_qt_syms34 * pinvG2_1 +
                          coefs_tq2_6 * coef_Jacob2_qt_syms33 * pinvG1_2 +
                          coefs_tq1_9 * coef_Jacob2_qt_syms23 * pinvG3_1 +
                          coefs_tq2_9 * coef_Jacob2_qt_syms22 * pinvG2_2 +
                          coefs_tq3_9 * coef_Jacob2_qt_syms21 * pinvG1_3 +
                          coefs_tq1_7 * coef_Jacob2_qt_syms30 * pinvG3_1 +
                          coefs_tq2_7 * coef_Jacob2_qt_syms29 * pinvG2_2 +
                          coefs_tq3_7 * coef_Jacob2_qt_syms28 * pinvG1_3 +
                          coefs_tq1_6 * coef_Jacob2_qt_syms35 * pinvG3_1 +
                          coefs_tq2_6 * coef_Jacob2_qt_syms34 * pinvG2_2 +
                          coefs_tq3_6 * coef_Jacob2_qt_syms33 * pinvG1_3 +
                          coefs_tq2_9 * coef_Jacob2_qt_syms23 * pinvG3_2 +
                          coefs_tq3_9 * coef_Jacob2_qt_syms22 * pinvG2_3 +
                          coefs_tq2_7 * coef_Jacob2_qt_syms30 * pinvG3_2 +
                          coefs_tq3_7 * coef_Jacob2_qt_syms29 * pinvG2_3 +
                          coefs_tq2_6 * coef_Jacob2_qt_syms35 * pinvG3_2 +
                          coefs_tq3_6 * coef_Jacob2_qt_syms34 * pinvG2_3 +
                          coefs_tq3_9 * coef_Jacob2_qt_syms23 * pinvG3_3 +
                          coefs_tq3_7 * coef_Jacob2_qt_syms30 * pinvG3_3 +
                          coefs_tq3_6 * coef_Jacob2_qt_syms35 * pinvG3_3;
    coef_f_q_sym(1, 16) = coef_Jacob2_qt_syms20 + coefs_tq1_7 * coef_Jacob2_qt_syms33 * pinvG1_1 +
                          coefs_tq1_7 * coef_Jacob2_qt_syms34 * pinvG2_1 +
                          coefs_tq2_7 * coef_Jacob2_qt_syms33 * pinvG1_2 +
                          coefs_tq1_7 * coef_Jacob2_qt_syms35 * pinvG3_1 +
                          coefs_tq2_7 * coef_Jacob2_qt_syms34 * pinvG2_2 +
                          coefs_tq3_7 * coef_Jacob2_qt_syms33 * pinvG1_3 +
                          coefs_tq2_7 * coef_Jacob2_qt_syms35 * pinvG3_2 +
                          coefs_tq3_7 * coef_Jacob2_qt_syms34 * pinvG2_3 +
                          coefs_tq3_7 * coef_Jacob2_qt_syms35 * pinvG3_3 +
                          coefs_tq1_10 * coef_Jacob2_qt_syms21 * pinvG1_1 +
                          coefs_tq1_10 * coef_Jacob2_qt_syms22 * pinvG2_1 +
                          coefs_tq1_10 * coef_Jacob2_qt_syms23 * pinvG3_1 +
                          coefs_tq2_10 * coef_Jacob2_qt_syms21 * pinvG1_2 +
                          coefs_tq2_10 * coef_Jacob2_qt_syms22 * pinvG2_2 +
                          coefs_tq2_10 * coef_Jacob2_qt_syms23 * pinvG3_2 +
                          coefs_tq3_10 * coef_Jacob2_qt_syms21 * pinvG1_3 +
                          coefs_tq3_10 * coef_Jacob2_qt_syms22 * pinvG2_3 +
                          coefs_tq3_10 * coef_Jacob2_qt_syms23 * pinvG3_3;
    coef_f_q_sym(1, 17) = coef_Jacob2_qt_syms24 + coefs_tq1_11 * coef_Jacob2_qt_syms21 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob2_qt_syms22 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob2_qt_syms23 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob2_qt_syms21 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob2_qt_syms22 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob2_qt_syms23 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob2_qt_syms21 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob2_qt_syms22 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob2_qt_syms23 * pinvG3_3;
    coef_f_q_sym(1, 18) = coef_Jacob2_qt_syms25 + coefs_tq1_8 * coef_Jacob2_qt_syms28 * pinvG1_1 +
                          coefs_tq1_8 * coef_Jacob2_qt_syms29 * pinvG2_1 +
                          coefs_tq2_8 * coef_Jacob2_qt_syms28 * pinvG1_2 +
                          coefs_tq1_8 * coef_Jacob2_qt_syms30 * pinvG3_1 +
                          coefs_tq2_8 * coef_Jacob2_qt_syms29 * pinvG2_2 +
                          coefs_tq3_8 * coef_Jacob2_qt_syms28 * pinvG1_3 +
                          coefs_tq2_8 * coef_Jacob2_qt_syms30 * pinvG3_2 +
                          coefs_tq3_8 * coef_Jacob2_qt_syms29 * pinvG2_3 +
                          coefs_tq3_8 * coef_Jacob2_qt_syms30 * pinvG3_3;
    coef_f_q_sym(1, 19) = coef_Jacob2_qt_syms26 + coefs_tq1_9 * coef_Jacob2_qt_syms28 * pinvG1_1 +
                          coefs_tq1_8 * coef_Jacob2_qt_syms33 * pinvG1_1 +
                          coefs_tq1_9 * coef_Jacob2_qt_syms29 * pinvG2_1 +
                          coefs_tq2_9 * coef_Jacob2_qt_syms28 * pinvG1_2 +
                          coefs_tq1_8 * coef_Jacob2_qt_syms34 * pinvG2_1 +
                          coefs_tq2_8 * coef_Jacob2_qt_syms33 * pinvG1_2 +
                          coefs_tq1_9 * coef_Jacob2_qt_syms30 * pinvG3_1 +
                          coefs_tq2_9 * coef_Jacob2_qt_syms29 * pinvG2_2 +
                          coefs_tq3_9 * coef_Jacob2_qt_syms28 * pinvG1_3 +
                          coefs_tq1_8 * coef_Jacob2_qt_syms35 * pinvG3_1 +
                          coefs_tq2_8 * coef_Jacob2_qt_syms34 * pinvG2_2 +
                          coefs_tq3_8 * coef_Jacob2_qt_syms33 * pinvG1_3 +
                          coefs_tq2_9 * coef_Jacob2_qt_syms30 * pinvG3_2 +
                          coefs_tq3_9 * coef_Jacob2_qt_syms29 * pinvG2_3 +
                          coefs_tq2_8 * coef_Jacob2_qt_syms35 * pinvG3_2 +
                          coefs_tq3_8 * coef_Jacob2_qt_syms34 * pinvG2_3 +
                          coefs_tq3_9 * coef_Jacob2_qt_syms30 * pinvG3_3 +
                          coefs_tq3_8 * coef_Jacob2_qt_syms35 * pinvG3_3;
    coef_f_q_sym(1, 20) = coef_Jacob2_qt_syms27 + coefs_tq1_9 * coef_Jacob2_qt_syms33 * pinvG1_1 +
                          coefs_tq1_9 * coef_Jacob2_qt_syms34 * pinvG2_1 +
                          coefs_tq2_9 * coef_Jacob2_qt_syms33 * pinvG1_2 +
                          coefs_tq1_9 * coef_Jacob2_qt_syms35 * pinvG3_1 +
                          coefs_tq2_9 * coef_Jacob2_qt_syms34 * pinvG2_2 +
                          coefs_tq3_9 * coef_Jacob2_qt_syms33 * pinvG1_3 +
                          coefs_tq2_9 * coef_Jacob2_qt_syms35 * pinvG3_2 +
                          coefs_tq3_9 * coef_Jacob2_qt_syms34 * pinvG2_3 +
                          coefs_tq3_9 * coef_Jacob2_qt_syms35 * pinvG3_3 +
                          coefs_tq1_10 * coef_Jacob2_qt_syms28 * pinvG1_1 +
                          coefs_tq1_10 * coef_Jacob2_qt_syms29 * pinvG2_1 +
                          coefs_tq1_10 * coef_Jacob2_qt_syms30 * pinvG3_1 +
                          coefs_tq2_10 * coef_Jacob2_qt_syms28 * pinvG1_2 +
                          coefs_tq2_10 * coef_Jacob2_qt_syms29 * pinvG2_2 +
                          coefs_tq2_10 * coef_Jacob2_qt_syms30 * pinvG3_2 +
                          coefs_tq3_10 * coef_Jacob2_qt_syms28 * pinvG1_3 +
                          coefs_tq3_10 * coef_Jacob2_qt_syms29 * pinvG2_3 +
                          coefs_tq3_10 * coef_Jacob2_qt_syms30 * pinvG3_3;
    coef_f_q_sym(1, 21) = coef_Jacob2_qt_syms31 + coefs_tq1_11 * coef_Jacob2_qt_syms28 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob2_qt_syms29 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob2_qt_syms30 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob2_qt_syms28 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob2_qt_syms29 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob2_qt_syms30 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob2_qt_syms28 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob2_qt_syms29 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob2_qt_syms30 * pinvG3_3;
    coef_f_q_sym(1, 22) = coef_Jacob2_qt_syms32 + coefs_tq1_10 * coef_Jacob2_qt_syms33 * pinvG1_1 +
                          coefs_tq1_10 * coef_Jacob2_qt_syms34 * pinvG2_1 +
                          coefs_tq1_10 * coef_Jacob2_qt_syms35 * pinvG3_1 +
                          coefs_tq2_10 * coef_Jacob2_qt_syms33 * pinvG1_2 +
                          coefs_tq2_10 * coef_Jacob2_qt_syms34 * pinvG2_2 +
                          coefs_tq2_10 * coef_Jacob2_qt_syms35 * pinvG3_2 +
                          coefs_tq3_10 * coef_Jacob2_qt_syms33 * pinvG1_3 +
                          coefs_tq3_10 * coef_Jacob2_qt_syms34 * pinvG2_3 +
                          coefs_tq3_10 * coef_Jacob2_qt_syms35 * pinvG3_3;
    coef_f_q_sym(1, 23) = coef_Jacob2_qt_syms36 + coefs_tq1_11 * coef_Jacob2_qt_syms33 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob2_qt_syms34 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob2_qt_syms35 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob2_qt_syms33 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob2_qt_syms34 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob2_qt_syms35 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob2_qt_syms33 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob2_qt_syms34 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob2_qt_syms35 * pinvG3_3;
    coef_f_q_sym(2, 0) = coef_Jacob3_qt_syms1 + coefs_tq1_1 * coef_Jacob3_qt_syms11 * pinvG1_1 +
                         coefs_tq1_1 * coef_Jacob3_qt_syms12 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob3_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob3_qt_syms13 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob3_qt_syms12 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob3_qt_syms11 * pinvG1_3 +
                         coefs_tq2_1 * coef_Jacob3_qt_syms13 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob3_qt_syms12 * pinvG2_3 +
                         coefs_tq3_1 * coef_Jacob3_qt_syms13 * pinvG3_3;
    coef_f_q_sym(2, 1) = coef_Jacob3_qt_syms2 + coefs_tq1_2 * coef_Jacob3_qt_syms11 * pinvG1_1 +
                         coefs_tq1_1 * coef_Jacob3_qt_syms21 * pinvG1_1 +
                         coefs_tq1_2 * coef_Jacob3_qt_syms12 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob3_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob3_qt_syms22 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob3_qt_syms21 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob3_qt_syms13 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob3_qt_syms12 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob3_qt_syms11 * pinvG1_3 +
                         coefs_tq1_1 * coef_Jacob3_qt_syms23 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob3_qt_syms22 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob3_qt_syms21 * pinvG1_3 +
                         coefs_tq2_2 * coef_Jacob3_qt_syms13 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob3_qt_syms12 * pinvG2_3 +
                         coefs_tq2_1 * coef_Jacob3_qt_syms23 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob3_qt_syms22 * pinvG2_3 +
                         coefs_tq3_2 * coef_Jacob3_qt_syms13 * pinvG3_3 +
                         coefs_tq3_1 * coef_Jacob3_qt_syms23 * pinvG3_3;
    coef_f_q_sym(2, 2) = coef_Jacob3_qt_syms3 + coefs_tq1_3 * coef_Jacob3_qt_syms11 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob3_qt_syms12 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob3_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob3_qt_syms28 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob3_qt_syms13 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob3_qt_syms12 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob3_qt_syms11 * pinvG1_3 +
                         coefs_tq1_1 * coef_Jacob3_qt_syms29 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob3_qt_syms28 * pinvG1_2 +
                         coefs_tq2_3 * coef_Jacob3_qt_syms13 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob3_qt_syms12 * pinvG2_3 +
                         coefs_tq1_1 * coef_Jacob3_qt_syms30 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob3_qt_syms29 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob3_qt_syms28 * pinvG1_3 +
                         coefs_tq3_3 * coef_Jacob3_qt_syms13 * pinvG3_3 +
                         coefs_tq2_1 * coef_Jacob3_qt_syms30 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob3_qt_syms29 * pinvG2_3 +
                         coefs_tq3_1 * coef_Jacob3_qt_syms30 * pinvG3_3;
    coef_f_q_sym(2, 3) = coef_Jacob3_qt_syms4 + coefs_tq1_4 * coef_Jacob3_qt_syms11 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob3_qt_syms12 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob3_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob3_qt_syms33 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob3_qt_syms13 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob3_qt_syms12 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob3_qt_syms11 * pinvG1_3 +
                         coefs_tq1_1 * coef_Jacob3_qt_syms34 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob3_qt_syms33 * pinvG1_2 +
                         coefs_tq2_4 * coef_Jacob3_qt_syms13 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob3_qt_syms12 * pinvG2_3 +
                         coefs_tq1_1 * coef_Jacob3_qt_syms35 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob3_qt_syms34 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob3_qt_syms33 * pinvG1_3 +
                         coefs_tq3_4 * coef_Jacob3_qt_syms13 * pinvG3_3 +
                         coefs_tq2_1 * coef_Jacob3_qt_syms35 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob3_qt_syms34 * pinvG2_3 +
                         coefs_tq3_1 * coef_Jacob3_qt_syms35 * pinvG3_3;
    coef_f_q_sym(2, 4) = coef_Jacob3_qt_syms5 + coefs_tq1_5 * coef_Jacob3_qt_syms11 * pinvG1_1 +
                         coefs_tq1_2 * coef_Jacob3_qt_syms21 * pinvG1_1 +
                         coefs_tq1_5 * coef_Jacob3_qt_syms12 * pinvG2_1 +
                         coefs_tq2_5 * coef_Jacob3_qt_syms11 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob3_qt_syms22 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob3_qt_syms21 * pinvG1_2 +
                         coefs_tq1_5 * coef_Jacob3_qt_syms13 * pinvG3_1 +
                         coefs_tq2_5 * coef_Jacob3_qt_syms12 * pinvG2_2 +
                         coefs_tq3_5 * coef_Jacob3_qt_syms11 * pinvG1_3 +
                         coefs_tq1_2 * coef_Jacob3_qt_syms23 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob3_qt_syms22 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob3_qt_syms21 * pinvG1_3 +
                         coefs_tq2_5 * coef_Jacob3_qt_syms13 * pinvG3_2 +
                         coefs_tq3_5 * coef_Jacob3_qt_syms12 * pinvG2_3 +
                         coefs_tq2_2 * coef_Jacob3_qt_syms23 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob3_qt_syms22 * pinvG2_3 +
                         coefs_tq3_5 * coef_Jacob3_qt_syms13 * pinvG3_3 +
                         coefs_tq3_2 * coef_Jacob3_qt_syms23 * pinvG3_3;
    coef_f_q_sym(2, 5) = coef_Jacob3_qt_syms6 + coefs_tq1_6 * coef_Jacob3_qt_syms11 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob3_qt_syms21 * pinvG1_1 +
                         coefs_tq1_6 * coef_Jacob3_qt_syms12 * pinvG2_1 +
                         coefs_tq2_6 * coef_Jacob3_qt_syms11 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob3_qt_syms28 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob3_qt_syms22 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob3_qt_syms21 * pinvG1_2 +
                         coefs_tq1_6 * coef_Jacob3_qt_syms13 * pinvG3_1 +
                         coefs_tq2_6 * coef_Jacob3_qt_syms12 * pinvG2_2 +
                         coefs_tq3_6 * coef_Jacob3_qt_syms11 * pinvG1_3 +
                         coefs_tq1_2 * coef_Jacob3_qt_syms29 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob3_qt_syms28 * pinvG1_2 +
                         coefs_tq1_3 * coef_Jacob3_qt_syms23 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob3_qt_syms22 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob3_qt_syms21 * pinvG1_3 +
                         coefs_tq2_6 * coef_Jacob3_qt_syms13 * pinvG3_2 +
                         coefs_tq3_6 * coef_Jacob3_qt_syms12 * pinvG2_3 +
                         coefs_tq1_2 * coef_Jacob3_qt_syms30 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob3_qt_syms29 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob3_qt_syms28 * pinvG1_3 +
                         coefs_tq2_3 * coef_Jacob3_qt_syms23 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob3_qt_syms22 * pinvG2_3 +
                         coefs_tq3_6 * coef_Jacob3_qt_syms13 * pinvG3_3 +
                         coefs_tq2_2 * coef_Jacob3_qt_syms30 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob3_qt_syms29 * pinvG2_3 +
                         coefs_tq3_3 * coef_Jacob3_qt_syms23 * pinvG3_3 +
                         coefs_tq3_2 * coef_Jacob3_qt_syms30 * pinvG3_3;
    coef_f_q_sym(2, 6) = coef_Jacob3_qt_syms7 + coefs_tq1_7 * coef_Jacob3_qt_syms11 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob3_qt_syms21 * pinvG1_1 +
                         coefs_tq1_7 * coef_Jacob3_qt_syms12 * pinvG2_1 +
                         coefs_tq2_7 * coef_Jacob3_qt_syms11 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob3_qt_syms33 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob3_qt_syms22 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob3_qt_syms21 * pinvG1_2 +
                         coefs_tq1_7 * coef_Jacob3_qt_syms13 * pinvG3_1 +
                         coefs_tq2_7 * coef_Jacob3_qt_syms12 * pinvG2_2 +
                         coefs_tq3_7 * coef_Jacob3_qt_syms11 * pinvG1_3 +
                         coefs_tq1_2 * coef_Jacob3_qt_syms34 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob3_qt_syms33 * pinvG1_2 +
                         coefs_tq1_4 * coef_Jacob3_qt_syms23 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob3_qt_syms22 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob3_qt_syms21 * pinvG1_3 +
                         coefs_tq2_7 * coef_Jacob3_qt_syms13 * pinvG3_2 +
                         coefs_tq3_7 * coef_Jacob3_qt_syms12 * pinvG2_3 +
                         coefs_tq1_2 * coef_Jacob3_qt_syms35 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob3_qt_syms34 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob3_qt_syms33 * pinvG1_3 +
                         coefs_tq2_4 * coef_Jacob3_qt_syms23 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob3_qt_syms22 * pinvG2_3 +
                         coefs_tq3_7 * coef_Jacob3_qt_syms13 * pinvG3_3 +
                         coefs_tq2_2 * coef_Jacob3_qt_syms35 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob3_qt_syms34 * pinvG2_3 +
                         coefs_tq3_4 * coef_Jacob3_qt_syms23 * pinvG3_3 +
                         coefs_tq3_2 * coef_Jacob3_qt_syms35 * pinvG3_3;
    coef_f_q_sym(2, 7) = coef_Jacob3_qt_syms8 + coefs_tq1_8 * coef_Jacob3_qt_syms11 * pinvG1_1 +
                         coefs_tq1_8 * coef_Jacob3_qt_syms12 * pinvG2_1 +
                         coefs_tq2_8 * coef_Jacob3_qt_syms11 * pinvG1_2 +
                         coefs_tq1_3 * coef_Jacob3_qt_syms28 * pinvG1_1 +
                         coefs_tq1_8 * coef_Jacob3_qt_syms13 * pinvG3_1 +
                         coefs_tq2_8 * coef_Jacob3_qt_syms12 * pinvG2_2 +
                         coefs_tq3_8 * coef_Jacob3_qt_syms11 * pinvG1_3 +
                         coefs_tq1_3 * coef_Jacob3_qt_syms29 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob3_qt_syms28 * pinvG1_2 +
                         coefs_tq2_8 * coef_Jacob3_qt_syms13 * pinvG3_2 +
                         coefs_tq3_8 * coef_Jacob3_qt_syms12 * pinvG2_3 +
                         coefs_tq1_3 * coef_Jacob3_qt_syms30 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob3_qt_syms29 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob3_qt_syms28 * pinvG1_3 +
                         coefs_tq3_8 * coef_Jacob3_qt_syms13 * pinvG3_3 +
                         coefs_tq2_3 * coef_Jacob3_qt_syms30 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob3_qt_syms29 * pinvG2_3 +
                         coefs_tq3_3 * coef_Jacob3_qt_syms30 * pinvG3_3;
    coef_f_q_sym(2, 8) = coef_Jacob3_qt_syms9 + coefs_tq1_9 * coef_Jacob3_qt_syms11 * pinvG1_1 +
                         coefs_tq1_9 * coef_Jacob3_qt_syms12 * pinvG2_1 +
                         coefs_tq2_9 * coef_Jacob3_qt_syms11 * pinvG1_2 +
                         coefs_tq1_4 * coef_Jacob3_qt_syms28 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob3_qt_syms33 * pinvG1_1 +
                         coefs_tq1_9 * coef_Jacob3_qt_syms13 * pinvG3_1 +
                         coefs_tq2_9 * coef_Jacob3_qt_syms12 * pinvG2_2 +
                         coefs_tq3_9 * coef_Jacob3_qt_syms11 * pinvG1_3 +
                         coefs_tq1_4 * coef_Jacob3_qt_syms29 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob3_qt_syms28 * pinvG1_2 +
                         coefs_tq1_3 * coef_Jacob3_qt_syms34 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob3_qt_syms33 * pinvG1_2 +
                         coefs_tq2_9 * coef_Jacob3_qt_syms13 * pinvG3_2 +
                         coefs_tq3_9 * coef_Jacob3_qt_syms12 * pinvG2_3 +
                         coefs_tq1_4 * coef_Jacob3_qt_syms30 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob3_qt_syms29 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob3_qt_syms28 * pinvG1_3 +
                         coefs_tq1_3 * coef_Jacob3_qt_syms35 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob3_qt_syms34 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob3_qt_syms33 * pinvG1_3 +
                         coefs_tq3_9 * coef_Jacob3_qt_syms13 * pinvG3_3 +
                         coefs_tq2_4 * coef_Jacob3_qt_syms30 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob3_qt_syms29 * pinvG2_3 +
                         coefs_tq2_3 * coef_Jacob3_qt_syms35 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob3_qt_syms34 * pinvG2_3 +
                         coefs_tq3_4 * coef_Jacob3_qt_syms30 * pinvG3_3 +
                         coefs_tq3_3 * coef_Jacob3_qt_syms35 * pinvG3_3;
    coef_f_q_sym(2, 9) = coef_Jacob3_qt_syms10 + coefs_tq1_4 * coef_Jacob3_qt_syms33 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob3_qt_syms34 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob3_qt_syms33 * pinvG1_2 +
                         coefs_tq1_4 * coef_Jacob3_qt_syms35 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob3_qt_syms34 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob3_qt_syms33 * pinvG1_3 +
                         coefs_tq2_4 * coef_Jacob3_qt_syms35 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob3_qt_syms34 * pinvG2_3 +
                         coefs_tq3_4 * coef_Jacob3_qt_syms35 * pinvG3_3 +
                         coefs_tq1_10 * coef_Jacob3_qt_syms11 * pinvG1_1 +
                         coefs_tq1_10 * coef_Jacob3_qt_syms12 * pinvG2_1 +
                         coefs_tq1_10 * coef_Jacob3_qt_syms13 * pinvG3_1 +
                         coefs_tq2_10 * coef_Jacob3_qt_syms11 * pinvG1_2 +
                         coefs_tq2_10 * coef_Jacob3_qt_syms12 * pinvG2_2 +
                         coefs_tq2_10 * coef_Jacob3_qt_syms13 * pinvG3_2 +
                         coefs_tq3_10 * coef_Jacob3_qt_syms11 * pinvG1_3 +
                         coefs_tq3_10 * coef_Jacob3_qt_syms12 * pinvG2_3 +
                         coefs_tq3_10 * coef_Jacob3_qt_syms13 * pinvG3_3;
    coef_f_q_sym(2, 10) = coef_Jacob3_qt_syms14 + coefs_tq1_11 * coef_Jacob3_qt_syms11 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob3_qt_syms12 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob3_qt_syms13 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob3_qt_syms11 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob3_qt_syms12 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob3_qt_syms13 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob3_qt_syms11 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob3_qt_syms12 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob3_qt_syms13 * pinvG3_3;
    coef_f_q_sym(2, 11) = coef_Jacob3_qt_syms15 + coefs_tq1_5 * coef_Jacob3_qt_syms21 * pinvG1_1 +
                          coefs_tq1_5 * coef_Jacob3_qt_syms22 * pinvG2_1 +
                          coefs_tq2_5 * coef_Jacob3_qt_syms21 * pinvG1_2 +
                          coefs_tq1_5 * coef_Jacob3_qt_syms23 * pinvG3_1 +
                          coefs_tq2_5 * coef_Jacob3_qt_syms22 * pinvG2_2 +
                          coefs_tq3_5 * coef_Jacob3_qt_syms21 * pinvG1_3 +
                          coefs_tq2_5 * coef_Jacob3_qt_syms23 * pinvG3_2 +
                          coefs_tq3_5 * coef_Jacob3_qt_syms22 * pinvG2_3 +
                          coefs_tq3_5 * coef_Jacob3_qt_syms23 * pinvG3_3;
    coef_f_q_sym(2, 12) = coef_Jacob3_qt_syms16 + coefs_tq1_6 * coef_Jacob3_qt_syms21 * pinvG1_1 +
                          coefs_tq1_5 * coef_Jacob3_qt_syms28 * pinvG1_1 +
                          coefs_tq1_6 * coef_Jacob3_qt_syms22 * pinvG2_1 +
                          coefs_tq2_6 * coef_Jacob3_qt_syms21 * pinvG1_2 +
                          coefs_tq1_5 * coef_Jacob3_qt_syms29 * pinvG2_1 +
                          coefs_tq2_5 * coef_Jacob3_qt_syms28 * pinvG1_2 +
                          coefs_tq1_6 * coef_Jacob3_qt_syms23 * pinvG3_1 +
                          coefs_tq2_6 * coef_Jacob3_qt_syms22 * pinvG2_2 +
                          coefs_tq3_6 * coef_Jacob3_qt_syms21 * pinvG1_3 +
                          coefs_tq1_5 * coef_Jacob3_qt_syms30 * pinvG3_1 +
                          coefs_tq2_5 * coef_Jacob3_qt_syms29 * pinvG2_2 +
                          coefs_tq3_5 * coef_Jacob3_qt_syms28 * pinvG1_3 +
                          coefs_tq2_6 * coef_Jacob3_qt_syms23 * pinvG3_2 +
                          coefs_tq3_6 * coef_Jacob3_qt_syms22 * pinvG2_3 +
                          coefs_tq2_5 * coef_Jacob3_qt_syms30 * pinvG3_2 +
                          coefs_tq3_5 * coef_Jacob3_qt_syms29 * pinvG2_3 +
                          coefs_tq3_6 * coef_Jacob3_qt_syms23 * pinvG3_3 +
                          coefs_tq3_5 * coef_Jacob3_qt_syms30 * pinvG3_3;
    coef_f_q_sym(2, 13) = coef_Jacob3_qt_syms17 + coefs_tq1_7 * coef_Jacob3_qt_syms21 * pinvG1_1 +
                          coefs_tq1_5 * coef_Jacob3_qt_syms33 * pinvG1_1 +
                          coefs_tq1_7 * coef_Jacob3_qt_syms22 * pinvG2_1 +
                          coefs_tq2_7 * coef_Jacob3_qt_syms21 * pinvG1_2 +
                          coefs_tq1_5 * coef_Jacob3_qt_syms34 * pinvG2_1 +
                          coefs_tq2_5 * coef_Jacob3_qt_syms33 * pinvG1_2 +
                          coefs_tq1_7 * coef_Jacob3_qt_syms23 * pinvG3_1 +
                          coefs_tq2_7 * coef_Jacob3_qt_syms22 * pinvG2_2 +
                          coefs_tq3_7 * coef_Jacob3_qt_syms21 * pinvG1_3 +
                          coefs_tq1_5 * coef_Jacob3_qt_syms35 * pinvG3_1 +
                          coefs_tq2_5 * coef_Jacob3_qt_syms34 * pinvG2_2 +
                          coefs_tq3_5 * coef_Jacob3_qt_syms33 * pinvG1_3 +
                          coefs_tq2_7 * coef_Jacob3_qt_syms23 * pinvG3_2 +
                          coefs_tq3_7 * coef_Jacob3_qt_syms22 * pinvG2_3 +
                          coefs_tq2_5 * coef_Jacob3_qt_syms35 * pinvG3_2 +
                          coefs_tq3_5 * coef_Jacob3_qt_syms34 * pinvG2_3 +
                          coefs_tq3_7 * coef_Jacob3_qt_syms23 * pinvG3_3 +
                          coefs_tq3_5 * coef_Jacob3_qt_syms35 * pinvG3_3;
    coef_f_q_sym(2, 14) = coef_Jacob3_qt_syms18 + coefs_tq1_8 * coef_Jacob3_qt_syms21 * pinvG1_1 +
                          coefs_tq1_6 * coef_Jacob3_qt_syms28 * pinvG1_1 +
                          coefs_tq1_8 * coef_Jacob3_qt_syms22 * pinvG2_1 +
                          coefs_tq2_8 * coef_Jacob3_qt_syms21 * pinvG1_2 +
                          coefs_tq1_6 * coef_Jacob3_qt_syms29 * pinvG2_1 +
                          coefs_tq2_6 * coef_Jacob3_qt_syms28 * pinvG1_2 +
                          coefs_tq1_8 * coef_Jacob3_qt_syms23 * pinvG3_1 +
                          coefs_tq2_8 * coef_Jacob3_qt_syms22 * pinvG2_2 +
                          coefs_tq3_8 * coef_Jacob3_qt_syms21 * pinvG1_3 +
                          coefs_tq1_6 * coef_Jacob3_qt_syms30 * pinvG3_1 +
                          coefs_tq2_6 * coef_Jacob3_qt_syms29 * pinvG2_2 +
                          coefs_tq3_6 * coef_Jacob3_qt_syms28 * pinvG1_3 +
                          coefs_tq2_8 * coef_Jacob3_qt_syms23 * pinvG3_2 +
                          coefs_tq3_8 * coef_Jacob3_qt_syms22 * pinvG2_3 +
                          coefs_tq2_6 * coef_Jacob3_qt_syms30 * pinvG3_2 +
                          coefs_tq3_6 * coef_Jacob3_qt_syms29 * pinvG2_3 +
                          coefs_tq3_8 * coef_Jacob3_qt_syms23 * pinvG3_3 +
                          coefs_tq3_6 * coef_Jacob3_qt_syms30 * pinvG3_3;
    coef_f_q_sym(2, 15) = coef_Jacob3_qt_syms19 + coefs_tq1_9 * coef_Jacob3_qt_syms21 * pinvG1_1 +
                          coefs_tq1_7 * coef_Jacob3_qt_syms28 * pinvG1_1 +
                          coefs_tq1_6 * coef_Jacob3_qt_syms33 * pinvG1_1 +
                          coefs_tq1_9 * coef_Jacob3_qt_syms22 * pinvG2_1 +
                          coefs_tq2_9 * coef_Jacob3_qt_syms21 * pinvG1_2 +
                          coefs_tq1_7 * coef_Jacob3_qt_syms29 * pinvG2_1 +
                          coefs_tq2_7 * coef_Jacob3_qt_syms28 * pinvG1_2 +
                          coefs_tq1_6 * coef_Jacob3_qt_syms34 * pinvG2_1 +
                          coefs_tq2_6 * coef_Jacob3_qt_syms33 * pinvG1_2 +
                          coefs_tq1_9 * coef_Jacob3_qt_syms23 * pinvG3_1 +
                          coefs_tq2_9 * coef_Jacob3_qt_syms22 * pinvG2_2 +
                          coefs_tq3_9 * coef_Jacob3_qt_syms21 * pinvG1_3 +
                          coefs_tq1_7 * coef_Jacob3_qt_syms30 * pinvG3_1 +
                          coefs_tq2_7 * coef_Jacob3_qt_syms29 * pinvG2_2 +
                          coefs_tq3_7 * coef_Jacob3_qt_syms28 * pinvG1_3 +
                          coefs_tq1_6 * coef_Jacob3_qt_syms35 * pinvG3_1 +
                          coefs_tq2_6 * coef_Jacob3_qt_syms34 * pinvG2_2 +
                          coefs_tq3_6 * coef_Jacob3_qt_syms33 * pinvG1_3 +
                          coefs_tq2_9 * coef_Jacob3_qt_syms23 * pinvG3_2 +
                          coefs_tq3_9 * coef_Jacob3_qt_syms22 * pinvG2_3 +
                          coefs_tq2_7 * coef_Jacob3_qt_syms30 * pinvG3_2 +
                          coefs_tq3_7 * coef_Jacob3_qt_syms29 * pinvG2_3 +
                          coefs_tq2_6 * coef_Jacob3_qt_syms35 * pinvG3_2 +
                          coefs_tq3_6 * coef_Jacob3_qt_syms34 * pinvG2_3 +
                          coefs_tq3_9 * coef_Jacob3_qt_syms23 * pinvG3_3 +
                          coefs_tq3_7 * coef_Jacob3_qt_syms30 * pinvG3_3 +
                          coefs_tq3_6 * coef_Jacob3_qt_syms35 * pinvG3_3;
    coef_f_q_sym(2, 16) = coef_Jacob3_qt_syms20 + coefs_tq1_7 * coef_Jacob3_qt_syms33 * pinvG1_1 +
                          coefs_tq1_7 * coef_Jacob3_qt_syms34 * pinvG2_1 +
                          coefs_tq2_7 * coef_Jacob3_qt_syms33 * pinvG1_2 +
                          coefs_tq1_7 * coef_Jacob3_qt_syms35 * pinvG3_1 +
                          coefs_tq2_7 * coef_Jacob3_qt_syms34 * pinvG2_2 +
                          coefs_tq3_7 * coef_Jacob3_qt_syms33 * pinvG1_3 +
                          coefs_tq2_7 * coef_Jacob3_qt_syms35 * pinvG3_2 +
                          coefs_tq3_7 * coef_Jacob3_qt_syms34 * pinvG2_3 +
                          coefs_tq3_7 * coef_Jacob3_qt_syms35 * pinvG3_3 +
                          coefs_tq1_10 * coef_Jacob3_qt_syms21 * pinvG1_1 +
                          coefs_tq1_10 * coef_Jacob3_qt_syms22 * pinvG2_1 +
                          coefs_tq1_10 * coef_Jacob3_qt_syms23 * pinvG3_1 +
                          coefs_tq2_10 * coef_Jacob3_qt_syms21 * pinvG1_2 +
                          coefs_tq2_10 * coef_Jacob3_qt_syms22 * pinvG2_2 +
                          coefs_tq2_10 * coef_Jacob3_qt_syms23 * pinvG3_2 +
                          coefs_tq3_10 * coef_Jacob3_qt_syms21 * pinvG1_3 +
                          coefs_tq3_10 * coef_Jacob3_qt_syms22 * pinvG2_3 +
                          coefs_tq3_10 * coef_Jacob3_qt_syms23 * pinvG3_3;
    coef_f_q_sym(2, 17) = coef_Jacob3_qt_syms24 + coefs_tq1_11 * coef_Jacob3_qt_syms21 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob3_qt_syms22 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob3_qt_syms23 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob3_qt_syms21 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob3_qt_syms22 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob3_qt_syms23 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob3_qt_syms21 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob3_qt_syms22 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob3_qt_syms23 * pinvG3_3;
    coef_f_q_sym(2, 18) = coef_Jacob3_qt_syms25 + coefs_tq1_8 * coef_Jacob3_qt_syms28 * pinvG1_1 +
                          coefs_tq1_8 * coef_Jacob3_qt_syms29 * pinvG2_1 +
                          coefs_tq2_8 * coef_Jacob3_qt_syms28 * pinvG1_2 +
                          coefs_tq1_8 * coef_Jacob3_qt_syms30 * pinvG3_1 +
                          coefs_tq2_8 * coef_Jacob3_qt_syms29 * pinvG2_2 +
                          coefs_tq3_8 * coef_Jacob3_qt_syms28 * pinvG1_3 +
                          coefs_tq2_8 * coef_Jacob3_qt_syms30 * pinvG3_2 +
                          coefs_tq3_8 * coef_Jacob3_qt_syms29 * pinvG2_3 +
                          coefs_tq3_8 * coef_Jacob3_qt_syms30 * pinvG3_3;
    coef_f_q_sym(2, 19) = coef_Jacob3_qt_syms26 + coefs_tq1_9 * coef_Jacob3_qt_syms28 * pinvG1_1 +
                          coefs_tq1_8 * coef_Jacob3_qt_syms33 * pinvG1_1 +
                          coefs_tq1_9 * coef_Jacob3_qt_syms29 * pinvG2_1 +
                          coefs_tq2_9 * coef_Jacob3_qt_syms28 * pinvG1_2 +
                          coefs_tq1_8 * coef_Jacob3_qt_syms34 * pinvG2_1 +
                          coefs_tq2_8 * coef_Jacob3_qt_syms33 * pinvG1_2 +
                          coefs_tq1_9 * coef_Jacob3_qt_syms30 * pinvG3_1 +
                          coefs_tq2_9 * coef_Jacob3_qt_syms29 * pinvG2_2 +
                          coefs_tq3_9 * coef_Jacob3_qt_syms28 * pinvG1_3 +
                          coefs_tq1_8 * coef_Jacob3_qt_syms35 * pinvG3_1 +
                          coefs_tq2_8 * coef_Jacob3_qt_syms34 * pinvG2_2 +
                          coefs_tq3_8 * coef_Jacob3_qt_syms33 * pinvG1_3 +
                          coefs_tq2_9 * coef_Jacob3_qt_syms30 * pinvG3_2 +
                          coefs_tq3_9 * coef_Jacob3_qt_syms29 * pinvG2_3 +
                          coefs_tq2_8 * coef_Jacob3_qt_syms35 * pinvG3_2 +
                          coefs_tq3_8 * coef_Jacob3_qt_syms34 * pinvG2_3 +
                          coefs_tq3_9 * coef_Jacob3_qt_syms30 * pinvG3_3 +
                          coefs_tq3_8 * coef_Jacob3_qt_syms35 * pinvG3_3;
    coef_f_q_sym(2, 20) = coef_Jacob3_qt_syms27 + coefs_tq1_9 * coef_Jacob3_qt_syms33 * pinvG1_1 +
                          coefs_tq1_9 * coef_Jacob3_qt_syms34 * pinvG2_1 +
                          coefs_tq2_9 * coef_Jacob3_qt_syms33 * pinvG1_2 +
                          coefs_tq1_9 * coef_Jacob3_qt_syms35 * pinvG3_1 +
                          coefs_tq2_9 * coef_Jacob3_qt_syms34 * pinvG2_2 +
                          coefs_tq3_9 * coef_Jacob3_qt_syms33 * pinvG1_3 +
                          coefs_tq2_9 * coef_Jacob3_qt_syms35 * pinvG3_2 +
                          coefs_tq3_9 * coef_Jacob3_qt_syms34 * pinvG2_3 +
                          coefs_tq3_9 * coef_Jacob3_qt_syms35 * pinvG3_3 +
                          coefs_tq1_10 * coef_Jacob3_qt_syms28 * pinvG1_1 +
                          coefs_tq1_10 * coef_Jacob3_qt_syms29 * pinvG2_1 +
                          coefs_tq1_10 * coef_Jacob3_qt_syms30 * pinvG3_1 +
                          coefs_tq2_10 * coef_Jacob3_qt_syms28 * pinvG1_2 +
                          coefs_tq2_10 * coef_Jacob3_qt_syms29 * pinvG2_2 +
                          coefs_tq2_10 * coef_Jacob3_qt_syms30 * pinvG3_2 +
                          coefs_tq3_10 * coef_Jacob3_qt_syms28 * pinvG1_3 +
                          coefs_tq3_10 * coef_Jacob3_qt_syms29 * pinvG2_3 +
                          coefs_tq3_10 * coef_Jacob3_qt_syms30 * pinvG3_3;
    coef_f_q_sym(2, 21) = coef_Jacob3_qt_syms31 + coefs_tq1_11 * coef_Jacob3_qt_syms28 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob3_qt_syms29 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob3_qt_syms30 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob3_qt_syms28 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob3_qt_syms29 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob3_qt_syms30 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob3_qt_syms28 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob3_qt_syms29 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob3_qt_syms30 * pinvG3_3;
    coef_f_q_sym(2, 22) = coef_Jacob3_qt_syms32 + coefs_tq1_10 * coef_Jacob3_qt_syms33 * pinvG1_1 +
                          coefs_tq1_10 * coef_Jacob3_qt_syms34 * pinvG2_1 +
                          coefs_tq1_10 * coef_Jacob3_qt_syms35 * pinvG3_1 +
                          coefs_tq2_10 * coef_Jacob3_qt_syms33 * pinvG1_2 +
                          coefs_tq2_10 * coef_Jacob3_qt_syms34 * pinvG2_2 +
                          coefs_tq2_10 * coef_Jacob3_qt_syms35 * pinvG3_2 +
                          coefs_tq3_10 * coef_Jacob3_qt_syms33 * pinvG1_3 +
                          coefs_tq3_10 * coef_Jacob3_qt_syms34 * pinvG2_3 +
                          coefs_tq3_10 * coef_Jacob3_qt_syms35 * pinvG3_3;
    coef_f_q_sym(2, 23) = coef_Jacob3_qt_syms36 + coefs_tq1_11 * coef_Jacob3_qt_syms33 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob3_qt_syms34 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob3_qt_syms35 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob3_qt_syms33 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob3_qt_syms34 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob3_qt_syms35 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob3_qt_syms33 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob3_qt_syms34 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob3_qt_syms35 * pinvG3_3;
    coef_f_q_sym(3, 0) = coef_Jacob4_qt_syms1 + coefs_tq1_1 * coef_Jacob4_qt_syms11 * pinvG1_1 +
                         coefs_tq1_1 * coef_Jacob4_qt_syms12 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob4_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob4_qt_syms13 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob4_qt_syms12 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob4_qt_syms11 * pinvG1_3 +
                         coefs_tq2_1 * coef_Jacob4_qt_syms13 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob4_qt_syms12 * pinvG2_3 +
                         coefs_tq3_1 * coef_Jacob4_qt_syms13 * pinvG3_3;
    coef_f_q_sym(3, 1) = coef_Jacob4_qt_syms2 + coefs_tq1_2 * coef_Jacob4_qt_syms11 * pinvG1_1 +
                         coefs_tq1_1 * coef_Jacob4_qt_syms21 * pinvG1_1 +
                         coefs_tq1_2 * coef_Jacob4_qt_syms12 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob4_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob4_qt_syms22 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob4_qt_syms21 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob4_qt_syms13 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob4_qt_syms12 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob4_qt_syms11 * pinvG1_3 +
                         coefs_tq1_1 * coef_Jacob4_qt_syms23 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob4_qt_syms22 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob4_qt_syms21 * pinvG1_3 +
                         coefs_tq2_2 * coef_Jacob4_qt_syms13 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob4_qt_syms12 * pinvG2_3 +
                         coefs_tq2_1 * coef_Jacob4_qt_syms23 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob4_qt_syms22 * pinvG2_3 +
                         coefs_tq3_2 * coef_Jacob4_qt_syms13 * pinvG3_3 +
                         coefs_tq3_1 * coef_Jacob4_qt_syms23 * pinvG3_3;
    coef_f_q_sym(3, 2) = coef_Jacob4_qt_syms3 + coefs_tq1_3 * coef_Jacob4_qt_syms11 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob4_qt_syms12 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob4_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob4_qt_syms28 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob4_qt_syms13 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob4_qt_syms12 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob4_qt_syms11 * pinvG1_3 +
                         coefs_tq1_1 * coef_Jacob4_qt_syms29 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob4_qt_syms28 * pinvG1_2 +
                         coefs_tq2_3 * coef_Jacob4_qt_syms13 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob4_qt_syms12 * pinvG2_3 +
                         coefs_tq1_1 * coef_Jacob4_qt_syms30 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob4_qt_syms29 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob4_qt_syms28 * pinvG1_3 +
                         coefs_tq3_3 * coef_Jacob4_qt_syms13 * pinvG3_3 +
                         coefs_tq2_1 * coef_Jacob4_qt_syms30 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob4_qt_syms29 * pinvG2_3 +
                         coefs_tq3_1 * coef_Jacob4_qt_syms30 * pinvG3_3;
    coef_f_q_sym(3, 3) = coef_Jacob4_qt_syms4 + coefs_tq1_4 * coef_Jacob4_qt_syms11 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob4_qt_syms12 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob4_qt_syms11 * pinvG1_2 +
                         coefs_tq1_1 * coef_Jacob4_qt_syms33 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob4_qt_syms13 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob4_qt_syms12 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob4_qt_syms11 * pinvG1_3 +
                         coefs_tq1_1 * coef_Jacob4_qt_syms34 * pinvG2_1 +
                         coefs_tq2_1 * coef_Jacob4_qt_syms33 * pinvG1_2 +
                         coefs_tq2_4 * coef_Jacob4_qt_syms13 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob4_qt_syms12 * pinvG2_3 +
                         coefs_tq1_1 * coef_Jacob4_qt_syms35 * pinvG3_1 +
                         coefs_tq2_1 * coef_Jacob4_qt_syms34 * pinvG2_2 +
                         coefs_tq3_1 * coef_Jacob4_qt_syms33 * pinvG1_3 +
                         coefs_tq3_4 * coef_Jacob4_qt_syms13 * pinvG3_3 +
                         coefs_tq2_1 * coef_Jacob4_qt_syms35 * pinvG3_2 +
                         coefs_tq3_1 * coef_Jacob4_qt_syms34 * pinvG2_3 +
                         coefs_tq3_1 * coef_Jacob4_qt_syms35 * pinvG3_3;
    coef_f_q_sym(3, 4) = coef_Jacob4_qt_syms5 + coefs_tq1_5 * coef_Jacob4_qt_syms11 * pinvG1_1 +
                         coefs_tq1_2 * coef_Jacob4_qt_syms21 * pinvG1_1 +
                         coefs_tq1_5 * coef_Jacob4_qt_syms12 * pinvG2_1 +
                         coefs_tq2_5 * coef_Jacob4_qt_syms11 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob4_qt_syms22 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob4_qt_syms21 * pinvG1_2 +
                         coefs_tq1_5 * coef_Jacob4_qt_syms13 * pinvG3_1 +
                         coefs_tq2_5 * coef_Jacob4_qt_syms12 * pinvG2_2 +
                         coefs_tq3_5 * coef_Jacob4_qt_syms11 * pinvG1_3 +
                         coefs_tq1_2 * coef_Jacob4_qt_syms23 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob4_qt_syms22 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob4_qt_syms21 * pinvG1_3 +
                         coefs_tq2_5 * coef_Jacob4_qt_syms13 * pinvG3_2 +
                         coefs_tq3_5 * coef_Jacob4_qt_syms12 * pinvG2_3 +
                         coefs_tq2_2 * coef_Jacob4_qt_syms23 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob4_qt_syms22 * pinvG2_3 +
                         coefs_tq3_5 * coef_Jacob4_qt_syms13 * pinvG3_3 +
                         coefs_tq3_2 * coef_Jacob4_qt_syms23 * pinvG3_3;
    coef_f_q_sym(3, 5) = coef_Jacob4_qt_syms6 + coefs_tq1_6 * coef_Jacob4_qt_syms11 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob4_qt_syms21 * pinvG1_1 +
                         coefs_tq1_6 * coef_Jacob4_qt_syms12 * pinvG2_1 +
                         coefs_tq2_6 * coef_Jacob4_qt_syms11 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob4_qt_syms28 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob4_qt_syms22 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob4_qt_syms21 * pinvG1_2 +
                         coefs_tq1_6 * coef_Jacob4_qt_syms13 * pinvG3_1 +
                         coefs_tq2_6 * coef_Jacob4_qt_syms12 * pinvG2_2 +
                         coefs_tq3_6 * coef_Jacob4_qt_syms11 * pinvG1_3 +
                         coefs_tq1_2 * coef_Jacob4_qt_syms29 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob4_qt_syms28 * pinvG1_2 +
                         coefs_tq1_3 * coef_Jacob4_qt_syms23 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob4_qt_syms22 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob4_qt_syms21 * pinvG1_3 +
                         coefs_tq2_6 * coef_Jacob4_qt_syms13 * pinvG3_2 +
                         coefs_tq3_6 * coef_Jacob4_qt_syms12 * pinvG2_3 +
                         coefs_tq1_2 * coef_Jacob4_qt_syms30 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob4_qt_syms29 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob4_qt_syms28 * pinvG1_3 +
                         coefs_tq2_3 * coef_Jacob4_qt_syms23 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob4_qt_syms22 * pinvG2_3 +
                         coefs_tq3_6 * coef_Jacob4_qt_syms13 * pinvG3_3 +
                         coefs_tq2_2 * coef_Jacob4_qt_syms30 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob4_qt_syms29 * pinvG2_3 +
                         coefs_tq3_3 * coef_Jacob4_qt_syms23 * pinvG3_3 +
                         coefs_tq3_2 * coef_Jacob4_qt_syms30 * pinvG3_3;
    coef_f_q_sym(3, 6) = coef_Jacob4_qt_syms7 + coefs_tq1_7 * coef_Jacob4_qt_syms11 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob4_qt_syms21 * pinvG1_1 +
                         coefs_tq1_7 * coef_Jacob4_qt_syms12 * pinvG2_1 +
                         coefs_tq2_7 * coef_Jacob4_qt_syms11 * pinvG1_2 +
                         coefs_tq1_2 * coef_Jacob4_qt_syms33 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob4_qt_syms22 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob4_qt_syms21 * pinvG1_2 +
                         coefs_tq1_7 * coef_Jacob4_qt_syms13 * pinvG3_1 +
                         coefs_tq2_7 * coef_Jacob4_qt_syms12 * pinvG2_2 +
                         coefs_tq3_7 * coef_Jacob4_qt_syms11 * pinvG1_3 +
                         coefs_tq1_2 * coef_Jacob4_qt_syms34 * pinvG2_1 +
                         coefs_tq2_2 * coef_Jacob4_qt_syms33 * pinvG1_2 +
                         coefs_tq1_4 * coef_Jacob4_qt_syms23 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob4_qt_syms22 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob4_qt_syms21 * pinvG1_3 +
                         coefs_tq2_7 * coef_Jacob4_qt_syms13 * pinvG3_2 +
                         coefs_tq3_7 * coef_Jacob4_qt_syms12 * pinvG2_3 +
                         coefs_tq1_2 * coef_Jacob4_qt_syms35 * pinvG3_1 +
                         coefs_tq2_2 * coef_Jacob4_qt_syms34 * pinvG2_2 +
                         coefs_tq3_2 * coef_Jacob4_qt_syms33 * pinvG1_3 +
                         coefs_tq2_4 * coef_Jacob4_qt_syms23 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob4_qt_syms22 * pinvG2_3 +
                         coefs_tq3_7 * coef_Jacob4_qt_syms13 * pinvG3_3 +
                         coefs_tq2_2 * coef_Jacob4_qt_syms35 * pinvG3_2 +
                         coefs_tq3_2 * coef_Jacob4_qt_syms34 * pinvG2_3 +
                         coefs_tq3_4 * coef_Jacob4_qt_syms23 * pinvG3_3 +
                         coefs_tq3_2 * coef_Jacob4_qt_syms35 * pinvG3_3;
    coef_f_q_sym(3, 7) = coef_Jacob4_qt_syms8 + coefs_tq1_8 * coef_Jacob4_qt_syms11 * pinvG1_1 +
                         coefs_tq1_8 * coef_Jacob4_qt_syms12 * pinvG2_1 +
                         coefs_tq2_8 * coef_Jacob4_qt_syms11 * pinvG1_2 +
                         coefs_tq1_3 * coef_Jacob4_qt_syms28 * pinvG1_1 +
                         coefs_tq1_8 * coef_Jacob4_qt_syms13 * pinvG3_1 +
                         coefs_tq2_8 * coef_Jacob4_qt_syms12 * pinvG2_2 +
                         coefs_tq3_8 * coef_Jacob4_qt_syms11 * pinvG1_3 +
                         coefs_tq1_3 * coef_Jacob4_qt_syms29 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob4_qt_syms28 * pinvG1_2 +
                         coefs_tq2_8 * coef_Jacob4_qt_syms13 * pinvG3_2 +
                         coefs_tq3_8 * coef_Jacob4_qt_syms12 * pinvG2_3 +
                         coefs_tq1_3 * coef_Jacob4_qt_syms30 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob4_qt_syms29 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob4_qt_syms28 * pinvG1_3 +
                         coefs_tq3_8 * coef_Jacob4_qt_syms13 * pinvG3_3 +
                         coefs_tq2_3 * coef_Jacob4_qt_syms30 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob4_qt_syms29 * pinvG2_3 +
                         coefs_tq3_3 * coef_Jacob4_qt_syms30 * pinvG3_3;
    coef_f_q_sym(3, 8) = coef_Jacob4_qt_syms9 + coefs_tq1_9 * coef_Jacob4_qt_syms11 * pinvG1_1 +
                         coefs_tq1_9 * coef_Jacob4_qt_syms12 * pinvG2_1 +
                         coefs_tq2_9 * coef_Jacob4_qt_syms11 * pinvG1_2 +
                         coefs_tq1_4 * coef_Jacob4_qt_syms28 * pinvG1_1 +
                         coefs_tq1_3 * coef_Jacob4_qt_syms33 * pinvG1_1 +
                         coefs_tq1_9 * coef_Jacob4_qt_syms13 * pinvG3_1 +
                         coefs_tq2_9 * coef_Jacob4_qt_syms12 * pinvG2_2 +
                         coefs_tq3_9 * coef_Jacob4_qt_syms11 * pinvG1_3 +
                         coefs_tq1_4 * coef_Jacob4_qt_syms29 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob4_qt_syms28 * pinvG1_2 +
                         coefs_tq1_3 * coef_Jacob4_qt_syms34 * pinvG2_1 +
                         coefs_tq2_3 * coef_Jacob4_qt_syms33 * pinvG1_2 +
                         coefs_tq2_9 * coef_Jacob4_qt_syms13 * pinvG3_2 +
                         coefs_tq3_9 * coef_Jacob4_qt_syms12 * pinvG2_3 +
                         coefs_tq1_4 * coef_Jacob4_qt_syms30 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob4_qt_syms29 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob4_qt_syms28 * pinvG1_3 +
                         coefs_tq1_3 * coef_Jacob4_qt_syms35 * pinvG3_1 +
                         coefs_tq2_3 * coef_Jacob4_qt_syms34 * pinvG2_2 +
                         coefs_tq3_3 * coef_Jacob4_qt_syms33 * pinvG1_3 +
                         coefs_tq3_9 * coef_Jacob4_qt_syms13 * pinvG3_3 +
                         coefs_tq2_4 * coef_Jacob4_qt_syms30 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob4_qt_syms29 * pinvG2_3 +
                         coefs_tq2_3 * coef_Jacob4_qt_syms35 * pinvG3_2 +
                         coefs_tq3_3 * coef_Jacob4_qt_syms34 * pinvG2_3 +
                         coefs_tq3_4 * coef_Jacob4_qt_syms30 * pinvG3_3 +
                         coefs_tq3_3 * coef_Jacob4_qt_syms35 * pinvG3_3;
    coef_f_q_sym(3, 9) = coef_Jacob4_qt_syms10 + coefs_tq1_4 * coef_Jacob4_qt_syms33 * pinvG1_1 +
                         coefs_tq1_4 * coef_Jacob4_qt_syms34 * pinvG2_1 +
                         coefs_tq2_4 * coef_Jacob4_qt_syms33 * pinvG1_2 +
                         coefs_tq1_4 * coef_Jacob4_qt_syms35 * pinvG3_1 +
                         coefs_tq2_4 * coef_Jacob4_qt_syms34 * pinvG2_2 +
                         coefs_tq3_4 * coef_Jacob4_qt_syms33 * pinvG1_3 +
                         coefs_tq2_4 * coef_Jacob4_qt_syms35 * pinvG3_2 +
                         coefs_tq3_4 * coef_Jacob4_qt_syms34 * pinvG2_3 +
                         coefs_tq3_4 * coef_Jacob4_qt_syms35 * pinvG3_3 +
                         coefs_tq1_10 * coef_Jacob4_qt_syms11 * pinvG1_1 +
                         coefs_tq1_10 * coef_Jacob4_qt_syms12 * pinvG2_1 +
                         coefs_tq1_10 * coef_Jacob4_qt_syms13 * pinvG3_1 +
                         coefs_tq2_10 * coef_Jacob4_qt_syms11 * pinvG1_2 +
                         coefs_tq2_10 * coef_Jacob4_qt_syms12 * pinvG2_2 +
                         coefs_tq2_10 * coef_Jacob4_qt_syms13 * pinvG3_2 +
                         coefs_tq3_10 * coef_Jacob4_qt_syms11 * pinvG1_3 +
                         coefs_tq3_10 * coef_Jacob4_qt_syms12 * pinvG2_3 +
                         coefs_tq3_10 * coef_Jacob4_qt_syms13 * pinvG3_3;
    coef_f_q_sym(3, 10) = coef_Jacob4_qt_syms14 + coefs_tq1_11 * coef_Jacob4_qt_syms11 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob4_qt_syms12 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob4_qt_syms13 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob4_qt_syms11 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob4_qt_syms12 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob4_qt_syms13 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob4_qt_syms11 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob4_qt_syms12 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob4_qt_syms13 * pinvG3_3;
    coef_f_q_sym(3, 11) = coef_Jacob4_qt_syms15 + coefs_tq1_5 * coef_Jacob4_qt_syms21 * pinvG1_1 +
                          coefs_tq1_5 * coef_Jacob4_qt_syms22 * pinvG2_1 +
                          coefs_tq2_5 * coef_Jacob4_qt_syms21 * pinvG1_2 +
                          coefs_tq1_5 * coef_Jacob4_qt_syms23 * pinvG3_1 +
                          coefs_tq2_5 * coef_Jacob4_qt_syms22 * pinvG2_2 +
                          coefs_tq3_5 * coef_Jacob4_qt_syms21 * pinvG1_3 +
                          coefs_tq2_5 * coef_Jacob4_qt_syms23 * pinvG3_2 +
                          coefs_tq3_5 * coef_Jacob4_qt_syms22 * pinvG2_3 +
                          coefs_tq3_5 * coef_Jacob4_qt_syms23 * pinvG3_3;
    coef_f_q_sym(3, 12) = coef_Jacob4_qt_syms16 + coefs_tq1_6 * coef_Jacob4_qt_syms21 * pinvG1_1 +
                          coefs_tq1_5 * coef_Jacob4_qt_syms28 * pinvG1_1 +
                          coefs_tq1_6 * coef_Jacob4_qt_syms22 * pinvG2_1 +
                          coefs_tq2_6 * coef_Jacob4_qt_syms21 * pinvG1_2 +
                          coefs_tq1_5 * coef_Jacob4_qt_syms29 * pinvG2_1 +
                          coefs_tq2_5 * coef_Jacob4_qt_syms28 * pinvG1_2 +
                          coefs_tq1_6 * coef_Jacob4_qt_syms23 * pinvG3_1 +
                          coefs_tq2_6 * coef_Jacob4_qt_syms22 * pinvG2_2 +
                          coefs_tq3_6 * coef_Jacob4_qt_syms21 * pinvG1_3 +
                          coefs_tq1_5 * coef_Jacob4_qt_syms30 * pinvG3_1 +
                          coefs_tq2_5 * coef_Jacob4_qt_syms29 * pinvG2_2 +
                          coefs_tq3_5 * coef_Jacob4_qt_syms28 * pinvG1_3 +
                          coefs_tq2_6 * coef_Jacob4_qt_syms23 * pinvG3_2 +
                          coefs_tq3_6 * coef_Jacob4_qt_syms22 * pinvG2_3 +
                          coefs_tq2_5 * coef_Jacob4_qt_syms30 * pinvG3_2 +
                          coefs_tq3_5 * coef_Jacob4_qt_syms29 * pinvG2_3 +
                          coefs_tq3_6 * coef_Jacob4_qt_syms23 * pinvG3_3 +
                          coefs_tq3_5 * coef_Jacob4_qt_syms30 * pinvG3_3;
    coef_f_q_sym(3, 13) = coef_Jacob4_qt_syms17 + coefs_tq1_7 * coef_Jacob4_qt_syms21 * pinvG1_1 +
                          coefs_tq1_5 * coef_Jacob4_qt_syms33 * pinvG1_1 +
                          coefs_tq1_7 * coef_Jacob4_qt_syms22 * pinvG2_1 +
                          coefs_tq2_7 * coef_Jacob4_qt_syms21 * pinvG1_2 +
                          coefs_tq1_5 * coef_Jacob4_qt_syms34 * pinvG2_1 +
                          coefs_tq2_5 * coef_Jacob4_qt_syms33 * pinvG1_2 +
                          coefs_tq1_7 * coef_Jacob4_qt_syms23 * pinvG3_1 +
                          coefs_tq2_7 * coef_Jacob4_qt_syms22 * pinvG2_2 +
                          coefs_tq3_7 * coef_Jacob4_qt_syms21 * pinvG1_3 +
                          coefs_tq1_5 * coef_Jacob4_qt_syms35 * pinvG3_1 +
                          coefs_tq2_5 * coef_Jacob4_qt_syms34 * pinvG2_2 +
                          coefs_tq3_5 * coef_Jacob4_qt_syms33 * pinvG1_3 +
                          coefs_tq2_7 * coef_Jacob4_qt_syms23 * pinvG3_2 +
                          coefs_tq3_7 * coef_Jacob4_qt_syms22 * pinvG2_3 +
                          coefs_tq2_5 * coef_Jacob4_qt_syms35 * pinvG3_2 +
                          coefs_tq3_5 * coef_Jacob4_qt_syms34 * pinvG2_3 +
                          coefs_tq3_7 * coef_Jacob4_qt_syms23 * pinvG3_3 +
                          coefs_tq3_5 * coef_Jacob4_qt_syms35 * pinvG3_3;
    coef_f_q_sym(3, 14) = coef_Jacob4_qt_syms18 + coefs_tq1_8 * coef_Jacob4_qt_syms21 * pinvG1_1 +
                          coefs_tq1_6 * coef_Jacob4_qt_syms28 * pinvG1_1 +
                          coefs_tq1_8 * coef_Jacob4_qt_syms22 * pinvG2_1 +
                          coefs_tq2_8 * coef_Jacob4_qt_syms21 * pinvG1_2 +
                          coefs_tq1_6 * coef_Jacob4_qt_syms29 * pinvG2_1 +
                          coefs_tq2_6 * coef_Jacob4_qt_syms28 * pinvG1_2 +
                          coefs_tq1_8 * coef_Jacob4_qt_syms23 * pinvG3_1 +
                          coefs_tq2_8 * coef_Jacob4_qt_syms22 * pinvG2_2 +
                          coefs_tq3_8 * coef_Jacob4_qt_syms21 * pinvG1_3 +
                          coefs_tq1_6 * coef_Jacob4_qt_syms30 * pinvG3_1 +
                          coefs_tq2_6 * coef_Jacob4_qt_syms29 * pinvG2_2 +
                          coefs_tq3_6 * coef_Jacob4_qt_syms28 * pinvG1_3 +
                          coefs_tq2_8 * coef_Jacob4_qt_syms23 * pinvG3_2 +
                          coefs_tq3_8 * coef_Jacob4_qt_syms22 * pinvG2_3 +
                          coefs_tq2_6 * coef_Jacob4_qt_syms30 * pinvG3_2 +
                          coefs_tq3_6 * coef_Jacob4_qt_syms29 * pinvG2_3 +
                          coefs_tq3_8 * coef_Jacob4_qt_syms23 * pinvG3_3 +
                          coefs_tq3_6 * coef_Jacob4_qt_syms30 * pinvG3_3;
    coef_f_q_sym(3, 15) = coef_Jacob4_qt_syms19 + coefs_tq1_9 * coef_Jacob4_qt_syms21 * pinvG1_1 +
                          coefs_tq1_7 * coef_Jacob4_qt_syms28 * pinvG1_1 +
                          coefs_tq1_6 * coef_Jacob4_qt_syms33 * pinvG1_1 +
                          coefs_tq1_9 * coef_Jacob4_qt_syms22 * pinvG2_1 +
                          coefs_tq2_9 * coef_Jacob4_qt_syms21 * pinvG1_2 +
                          coefs_tq1_7 * coef_Jacob4_qt_syms29 * pinvG2_1 +
                          coefs_tq2_7 * coef_Jacob4_qt_syms28 * pinvG1_2 +
                          coefs_tq1_6 * coef_Jacob4_qt_syms34 * pinvG2_1 +
                          coefs_tq2_6 * coef_Jacob4_qt_syms33 * pinvG1_2 +
                          coefs_tq1_9 * coef_Jacob4_qt_syms23 * pinvG3_1 +
                          coefs_tq2_9 * coef_Jacob4_qt_syms22 * pinvG2_2 +
                          coefs_tq3_9 * coef_Jacob4_qt_syms21 * pinvG1_3 +
                          coefs_tq1_7 * coef_Jacob4_qt_syms30 * pinvG3_1 +
                          coefs_tq2_7 * coef_Jacob4_qt_syms29 * pinvG2_2 +
                          coefs_tq3_7 * coef_Jacob4_qt_syms28 * pinvG1_3 +
                          coefs_tq1_6 * coef_Jacob4_qt_syms35 * pinvG3_1 +
                          coefs_tq2_6 * coef_Jacob4_qt_syms34 * pinvG2_2 +
                          coefs_tq3_6 * coef_Jacob4_qt_syms33 * pinvG1_3 +
                          coefs_tq2_9 * coef_Jacob4_qt_syms23 * pinvG3_2 +
                          coefs_tq3_9 * coef_Jacob4_qt_syms22 * pinvG2_3 +
                          coefs_tq2_7 * coef_Jacob4_qt_syms30 * pinvG3_2 +
                          coefs_tq3_7 * coef_Jacob4_qt_syms29 * pinvG2_3 +
                          coefs_tq2_6 * coef_Jacob4_qt_syms35 * pinvG3_2 +
                          coefs_tq3_6 * coef_Jacob4_qt_syms34 * pinvG2_3 +
                          coefs_tq3_9 * coef_Jacob4_qt_syms23 * pinvG3_3 +
                          coefs_tq3_7 * coef_Jacob4_qt_syms30 * pinvG3_3 +
                          coefs_tq3_6 * coef_Jacob4_qt_syms35 * pinvG3_3;
    coef_f_q_sym(3, 16) = coef_Jacob4_qt_syms20 + coefs_tq1_7 * coef_Jacob4_qt_syms33 * pinvG1_1 +
                          coefs_tq1_7 * coef_Jacob4_qt_syms34 * pinvG2_1 +
                          coefs_tq2_7 * coef_Jacob4_qt_syms33 * pinvG1_2 +
                          coefs_tq1_7 * coef_Jacob4_qt_syms35 * pinvG3_1 +
                          coefs_tq2_7 * coef_Jacob4_qt_syms34 * pinvG2_2 +
                          coefs_tq3_7 * coef_Jacob4_qt_syms33 * pinvG1_3 +
                          coefs_tq2_7 * coef_Jacob4_qt_syms35 * pinvG3_2 +
                          coefs_tq3_7 * coef_Jacob4_qt_syms34 * pinvG2_3 +
                          coefs_tq3_7 * coef_Jacob4_qt_syms35 * pinvG3_3 +
                          coefs_tq1_10 * coef_Jacob4_qt_syms21 * pinvG1_1 +
                          coefs_tq1_10 * coef_Jacob4_qt_syms22 * pinvG2_1 +
                          coefs_tq1_10 * coef_Jacob4_qt_syms23 * pinvG3_1 +
                          coefs_tq2_10 * coef_Jacob4_qt_syms21 * pinvG1_2 +
                          coefs_tq2_10 * coef_Jacob4_qt_syms22 * pinvG2_2 +
                          coefs_tq2_10 * coef_Jacob4_qt_syms23 * pinvG3_2 +
                          coefs_tq3_10 * coef_Jacob4_qt_syms21 * pinvG1_3 +
                          coefs_tq3_10 * coef_Jacob4_qt_syms22 * pinvG2_3 +
                          coefs_tq3_10 * coef_Jacob4_qt_syms23 * pinvG3_3;
    coef_f_q_sym(3, 17) = coef_Jacob4_qt_syms24 + coefs_tq1_11 * coef_Jacob4_qt_syms21 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob4_qt_syms22 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob4_qt_syms23 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob4_qt_syms21 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob4_qt_syms22 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob4_qt_syms23 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob4_qt_syms21 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob4_qt_syms22 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob4_qt_syms23 * pinvG3_3;
    coef_f_q_sym(3, 18) = coef_Jacob4_qt_syms25 + coefs_tq1_8 * coef_Jacob4_qt_syms28 * pinvG1_1 +
                          coefs_tq1_8 * coef_Jacob4_qt_syms29 * pinvG2_1 +
                          coefs_tq2_8 * coef_Jacob4_qt_syms28 * pinvG1_2 +
                          coefs_tq1_8 * coef_Jacob4_qt_syms30 * pinvG3_1 +
                          coefs_tq2_8 * coef_Jacob4_qt_syms29 * pinvG2_2 +
                          coefs_tq3_8 * coef_Jacob4_qt_syms28 * pinvG1_3 +
                          coefs_tq2_8 * coef_Jacob4_qt_syms30 * pinvG3_2 +
                          coefs_tq3_8 * coef_Jacob4_qt_syms29 * pinvG2_3 +
                          coefs_tq3_8 * coef_Jacob4_qt_syms30 * pinvG3_3;
    coef_f_q_sym(3, 19) = coef_Jacob4_qt_syms26 + coefs_tq1_9 * coef_Jacob4_qt_syms28 * pinvG1_1 +
                          coefs_tq1_8 * coef_Jacob4_qt_syms33 * pinvG1_1 +
                          coefs_tq1_9 * coef_Jacob4_qt_syms29 * pinvG2_1 +
                          coefs_tq2_9 * coef_Jacob4_qt_syms28 * pinvG1_2 +
                          coefs_tq1_8 * coef_Jacob4_qt_syms34 * pinvG2_1 +
                          coefs_tq2_8 * coef_Jacob4_qt_syms33 * pinvG1_2 +
                          coefs_tq1_9 * coef_Jacob4_qt_syms30 * pinvG3_1 +
                          coefs_tq2_9 * coef_Jacob4_qt_syms29 * pinvG2_2 +
                          coefs_tq3_9 * coef_Jacob4_qt_syms28 * pinvG1_3 +
                          coefs_tq1_8 * coef_Jacob4_qt_syms35 * pinvG3_1 +
                          coefs_tq2_8 * coef_Jacob4_qt_syms34 * pinvG2_2 +
                          coefs_tq3_8 * coef_Jacob4_qt_syms33 * pinvG1_3 +
                          coefs_tq2_9 * coef_Jacob4_qt_syms30 * pinvG3_2 +
                          coefs_tq3_9 * coef_Jacob4_qt_syms29 * pinvG2_3 +
                          coefs_tq2_8 * coef_Jacob4_qt_syms35 * pinvG3_2 +
                          coefs_tq3_8 * coef_Jacob4_qt_syms34 * pinvG2_3 +
                          coefs_tq3_9 * coef_Jacob4_qt_syms30 * pinvG3_3 +
                          coefs_tq3_8 * coef_Jacob4_qt_syms35 * pinvG3_3;
    coef_f_q_sym(3, 20) = coef_Jacob4_qt_syms27 + coefs_tq1_9 * coef_Jacob4_qt_syms33 * pinvG1_1 +
                          coefs_tq1_9 * coef_Jacob4_qt_syms34 * pinvG2_1 +
                          coefs_tq2_9 * coef_Jacob4_qt_syms33 * pinvG1_2 +
                          coefs_tq1_9 * coef_Jacob4_qt_syms35 * pinvG3_1 +
                          coefs_tq2_9 * coef_Jacob4_qt_syms34 * pinvG2_2 +
                          coefs_tq3_9 * coef_Jacob4_qt_syms33 * pinvG1_3 +
                          coefs_tq2_9 * coef_Jacob4_qt_syms35 * pinvG3_2 +
                          coefs_tq3_9 * coef_Jacob4_qt_syms34 * pinvG2_3 +
                          coefs_tq3_9 * coef_Jacob4_qt_syms35 * pinvG3_3 +
                          coefs_tq1_10 * coef_Jacob4_qt_syms28 * pinvG1_1 +
                          coefs_tq1_10 * coef_Jacob4_qt_syms29 * pinvG2_1 +
                          coefs_tq1_10 * coef_Jacob4_qt_syms30 * pinvG3_1 +
                          coefs_tq2_10 * coef_Jacob4_qt_syms28 * pinvG1_2 +
                          coefs_tq2_10 * coef_Jacob4_qt_syms29 * pinvG2_2 +
                          coefs_tq2_10 * coef_Jacob4_qt_syms30 * pinvG3_2 +
                          coefs_tq3_10 * coef_Jacob4_qt_syms28 * pinvG1_3 +
                          coefs_tq3_10 * coef_Jacob4_qt_syms29 * pinvG2_3 +
                          coefs_tq3_10 * coef_Jacob4_qt_syms30 * pinvG3_3;
    coef_f_q_sym(3, 21) = coef_Jacob4_qt_syms31 + coefs_tq1_11 * coef_Jacob4_qt_syms28 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob4_qt_syms29 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob4_qt_syms30 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob4_qt_syms28 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob4_qt_syms29 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob4_qt_syms30 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob4_qt_syms28 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob4_qt_syms29 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob4_qt_syms30 * pinvG3_3;
    coef_f_q_sym(3, 22) = coef_Jacob4_qt_syms32 + coefs_tq1_10 * coef_Jacob4_qt_syms33 * pinvG1_1 +
                          coefs_tq1_10 * coef_Jacob4_qt_syms34 * pinvG2_1 +
                          coefs_tq1_10 * coef_Jacob4_qt_syms35 * pinvG3_1 +
                          coefs_tq2_10 * coef_Jacob4_qt_syms33 * pinvG1_2 +
                          coefs_tq2_10 * coef_Jacob4_qt_syms34 * pinvG2_2 +
                          coefs_tq2_10 * coef_Jacob4_qt_syms35 * pinvG3_2 +
                          coefs_tq3_10 * coef_Jacob4_qt_syms33 * pinvG1_3 +
                          coefs_tq3_10 * coef_Jacob4_qt_syms34 * pinvG2_3 +
                          coefs_tq3_10 * coef_Jacob4_qt_syms35 * pinvG3_3;
    coef_f_q_sym(3, 23) = coef_Jacob4_qt_syms36 + coefs_tq1_11 * coef_Jacob4_qt_syms33 * pinvG1_1 +
                          coefs_tq1_11 * coef_Jacob4_qt_syms34 * pinvG2_1 +
                          coefs_tq1_11 * coef_Jacob4_qt_syms35 * pinvG3_1 +
                          coefs_tq2_11 * coef_Jacob4_qt_syms33 * pinvG1_2 +
                          coefs_tq2_11 * coef_Jacob4_qt_syms34 * pinvG2_2 +
                          coefs_tq2_11 * coef_Jacob4_qt_syms35 * pinvG3_2 +
                          coefs_tq3_11 * coef_Jacob4_qt_syms33 * pinvG1_3 +
                          coefs_tq3_11 * coef_Jacob4_qt_syms34 * pinvG2_3 +
                          coefs_tq3_11 * coef_Jacob4_qt_syms35 * pinvG3_3;


    W(0, 0) = -coef_Jacob1_qt_syms1 - coefs_tq1_1 * coef_Jacob1_qt_syms11 * pinvG1_1 -
              coefs_tq1_1 * coef_Jacob1_qt_syms12 * pinvG2_1 - coefs_tq2_1 * coef_Jacob1_qt_syms11 * pinvG1_2 -
              coefs_tq1_1 * coef_Jacob1_qt_syms13 * pinvG3_1 - coefs_tq2_1 * coef_Jacob1_qt_syms12 * pinvG2_2 -
              coefs_tq3_1 * coef_Jacob1_qt_syms11 * pinvG1_3 - coefs_tq2_1 * coef_Jacob1_qt_syms13 * pinvG3_2 -
              coefs_tq3_1 * coef_Jacob1_qt_syms12 * pinvG2_3 - coefs_tq3_1 * coef_Jacob1_qt_syms13 * pinvG3_3;
    W(0, 1) = coef_Jacob1_qt_syms2 * (-1.0 * 2.0) - (coefs_tq1_2 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_1 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
              (coefs_tq1_2 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
              (coefs_tq1_2 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
              (coefs_tq2_2 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
              (coefs_tq3_2 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0;
    W(0, 2) = coef_Jacob1_qt_syms3 * (-1.0 * 2.0) - (coefs_tq1_3 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
              (coefs_tq2_3 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
              (coefs_tq3_3 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0;
    W(0, 3) = coef_Jacob1_qt_syms4 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_4 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_4 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
              (coefs_tq1_4 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_4 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_4 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
              (coefs_tq2_4 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_4 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
              (coefs_tq3_4 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0;
    W(0, 4) = coef_Jacob1_qt_syms2 * (-1.0 * 2.0) - (coefs_tq1_2 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_1 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
              (coefs_tq1_2 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
              (coefs_tq1_2 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
              (coefs_tq2_2 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
              (coefs_tq3_2 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0;
    W(0, 5) = coef_Jacob1_qt_syms5 * (-1.0 * 2.0) - (coefs_tq1_5 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_2 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
              (coefs_tq1_5 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_5 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_2 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
              (coefs_tq1_5 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_5 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_5 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_2 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
              (coefs_tq2_5 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_5 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq2_2 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
              (coefs_tq3_5 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq3_2 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0;
    W(0, 6) = coef_Jacob1_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
              (coefs_tq1_6 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
              (coefs_tq2_6 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
              (coefs_tq1_2 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
              (coefs_tq2_3 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
              (coefs_tq1_6 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
              (coefs_tq2_6 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
              (coefs_tq3_6 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
              (coefs_tq1_2 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
              (coefs_tq2_2 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
              (coefs_tq1_3 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
              (coefs_tq2_3 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
              (coefs_tq3_3 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
              (coefs_tq2_6 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
              (coefs_tq3_6 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
              (coefs_tq1_2 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
              (coefs_tq2_2 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
              (coefs_tq3_2 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
              (coefs_tq2_3 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
              (coefs_tq3_3 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
              (coefs_tq3_6 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
              (coefs_tq2_2 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
              (coefs_tq3_2 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
              (coefs_tq3_3 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
              (coefs_tq3_2 * coef_Jacob1_qt_syms30 * pinvG3_3) ;
    W(0, 7) = coef_Jacob1_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
              (coefs_tq1_4 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
              (coefs_tq1_7 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
              (coefs_tq2_7 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
              (coefs_tq1_2 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
              (coefs_tq1_4 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
              (coefs_tq2_4 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
              (coefs_tq1_7 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
              (coefs_tq2_7 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
              (coefs_tq3_7 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
              (coefs_tq1_2 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
              (coefs_tq2_2 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
              (coefs_tq1_4 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
              (coefs_tq2_4 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
              (coefs_tq3_4 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
              (coefs_tq2_7 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
              (coefs_tq3_7 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
              (coefs_tq1_2 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
              (coefs_tq2_2 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
              (coefs_tq3_2 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
              (coefs_tq2_4 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
              (coefs_tq3_4 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
              (coefs_tq3_7 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
              (coefs_tq2_2 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
              (coefs_tq3_2 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
              (coefs_tq3_4 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
              (coefs_tq3_2 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 8) = coef_Jacob1_qt_syms3 * (-1.0 * 2.0) - (coefs_tq1_3 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
              (coefs_tq2_3 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
              (coefs_tq3_3 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0;
    W(0, 9) = coef_Jacob1_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
              (coefs_tq1_6 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
              (coefs_tq2_6 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
              (coefs_tq1_2 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
              (coefs_tq2_3 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
              (coefs_tq1_6 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
              (coefs_tq2_6 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
              (coefs_tq3_6 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
              (coefs_tq1_2 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
              (coefs_tq2_2 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
              (coefs_tq1_3 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
              (coefs_tq2_3 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
              (coefs_tq3_3 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
              (coefs_tq2_6 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
              (coefs_tq3_6 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
              (coefs_tq1_2 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
              (coefs_tq2_2 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
              (coefs_tq3_2 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
              (coefs_tq2_3 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
              (coefs_tq3_3 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
              (coefs_tq3_6 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
              (coefs_tq2_2 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
              (coefs_tq3_2 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
              (coefs_tq3_3 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
              (coefs_tq3_2 * coef_Jacob1_qt_syms30 * pinvG3_3) ;
    W(0, 10) = coef_Jacob1_qt_syms8 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_3 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_3 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0;
    W(0, 11) = coef_Jacob1_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 12) = coef_Jacob1_qt_syms4 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq2_4 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0;
    W(0, 13) = coef_Jacob1_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 14) = coef_Jacob1_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 15) = coef_Jacob1_qt_syms10 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_4 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_4 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0;
    W(0, 16) = coef_Jacob1_qt_syms2 * (-1.0 * 2.0) - (coefs_tq1_2 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_1 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_2 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_2 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq2_2 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_2 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0;
    W(0, 17) = coef_Jacob1_qt_syms5 * (-1.0 * 2.0) - (coefs_tq1_5 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_2 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_2 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_2 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq2_2 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq3_2 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0;
    W(0, 18) = coef_Jacob1_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms30 * pinvG3_3) ;
    W(0, 19) = coef_Jacob1_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 20) = coef_Jacob1_qt_syms5 * (-1.0 * 2.0) - (coefs_tq1_5 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_2 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_2 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_2 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq2_2 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq3_2 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0;
    W(0, 21) = -coef_Jacob1_qt_syms15 - coefs_tq1_5 * coef_Jacob1_qt_syms21 * pinvG1_1 -
               coefs_tq1_5 * coef_Jacob1_qt_syms22 * pinvG2_1 - coefs_tq2_5 * coef_Jacob1_qt_syms21 * pinvG1_2 -
               coefs_tq1_5 * coef_Jacob1_qt_syms23 * pinvG3_1 - coefs_tq2_5 * coef_Jacob1_qt_syms22 * pinvG2_2 -
               coefs_tq3_5 * coef_Jacob1_qt_syms21 * pinvG1_3 - coefs_tq2_5 * coef_Jacob1_qt_syms23 * pinvG3_2 -
               coefs_tq3_5 * coef_Jacob1_qt_syms22 * pinvG2_3 - coefs_tq3_5 * coef_Jacob1_qt_syms23 * pinvG3_3;
    W(0, 22) = coef_Jacob1_qt_syms16 * (-1.0 * 2.0) - (coefs_tq1_6 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0;
    W(0, 23) = coef_Jacob1_qt_syms17 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0;
    W(0, 24) = coef_Jacob1_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms30 * pinvG3_3) ;
    W(0, 25) = coef_Jacob1_qt_syms16 * (-1.0 * 2.0) - (coefs_tq1_6 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0;
    W(0, 26) = coef_Jacob1_qt_syms18 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_6 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0;
    W(0, 27) = coef_Jacob1_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 28) = coef_Jacob1_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 29) = coef_Jacob1_qt_syms17 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0;
    W(0, 30) = coef_Jacob1_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 31) = coef_Jacob1_qt_syms20 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0;
    W(0, 32) = coef_Jacob1_qt_syms3 * (-1.0 * 2.0) - (coefs_tq1_3 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_3 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_3 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_3 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0;
    W(0, 33) = coef_Jacob1_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms30 * pinvG3_3) ;
    W(0, 34) = coef_Jacob1_qt_syms8 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_3 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_3 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0;
    W(0, 35) = coef_Jacob1_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 36) = coef_Jacob1_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms30 * pinvG3_3) ;
    W(0, 37) = coef_Jacob1_qt_syms16 * (-1.0 * 2.0) - (coefs_tq1_6 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0;
    W(0, 38) = coef_Jacob1_qt_syms18 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_6 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0;
    W(0, 39) = coef_Jacob1_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 40) = coef_Jacob1_qt_syms8 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_3 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_3 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0;
    W(0, 41) = coef_Jacob1_qt_syms18 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_6 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0;
    W(0, 42) = -coef_Jacob1_qt_syms25 - coefs_tq1_8 * coef_Jacob1_qt_syms28 * pinvG1_1 -
               coefs_tq1_8 * coef_Jacob1_qt_syms29 * pinvG2_1 - coefs_tq2_8 * coef_Jacob1_qt_syms28 * pinvG1_2 -
               coefs_tq1_8 * coef_Jacob1_qt_syms30 * pinvG3_1 - coefs_tq2_8 * coef_Jacob1_qt_syms29 * pinvG2_2 -
               coefs_tq3_8 * coef_Jacob1_qt_syms28 * pinvG1_3 - coefs_tq2_8 * coef_Jacob1_qt_syms30 * pinvG3_2 -
               coefs_tq3_8 * coef_Jacob1_qt_syms29 * pinvG2_3 - coefs_tq3_8 * coef_Jacob1_qt_syms30 * pinvG3_3;
    W(0, 43) = coef_Jacob1_qt_syms26 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0;
    W(0, 44) = coef_Jacob1_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 45) = coef_Jacob1_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 46) = coef_Jacob1_qt_syms26 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0;
    W(0, 47) = coef_Jacob1_qt_syms27 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0;
    W(0, 48) = coef_Jacob1_qt_syms4 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq2_4 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0;
    W(0, 49) = coef_Jacob1_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 50) = coef_Jacob1_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 51) = coef_Jacob1_qt_syms10 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_4 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_4 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0;
    W(0, 52) = coef_Jacob1_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 53) = coef_Jacob1_qt_syms17 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0;
    W(0, 54) = coef_Jacob1_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 55) = coef_Jacob1_qt_syms20 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0;
    W(0, 56) = coef_Jacob1_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob1_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob1_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 57) = coef_Jacob1_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob1_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob1_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob1_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob1_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob1_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob1_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob1_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob1_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob1_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob1_qt_syms35 * pinvG3_3) ;
    W(0, 58) = coef_Jacob1_qt_syms26 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq1_8 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0;
    W(0, 59) = coef_Jacob1_qt_syms27 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0;
    W(0, 60) = coef_Jacob1_qt_syms10 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_4 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_4 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms13 * pinvG3_3) * 2.0;
    W(0, 61) = coef_Jacob1_qt_syms20 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms23 * pinvG3_3) * 2.0;
    W(0, 62) = coef_Jacob1_qt_syms27 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob1_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob1_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob1_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob1_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob1_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob1_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob1_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob1_qt_syms30 * pinvG3_3) * 2.0;
    W(0, 63) = -coef_Jacob1_qt_syms32 - coefs_tq1_10 * coef_Jacob1_qt_syms33 * pinvG1_1 -
               coefs_tq1_10 * coef_Jacob1_qt_syms34 * pinvG2_1 - coefs_tq1_10 * coef_Jacob1_qt_syms35 * pinvG3_1 -
               coefs_tq2_10 * coef_Jacob1_qt_syms33 * pinvG1_2 - coefs_tq2_10 * coef_Jacob1_qt_syms34 * pinvG2_2 -
               coefs_tq2_10 * coef_Jacob1_qt_syms35 * pinvG3_2 - coefs_tq3_10 * coef_Jacob1_qt_syms33 * pinvG1_3 -
               coefs_tq3_10 * coef_Jacob1_qt_syms34 * pinvG2_3 - coefs_tq3_10 * coef_Jacob1_qt_syms35 * pinvG3_3;
    W(1, 0) = -coef_Jacob2_qt_syms1 - coefs_tq1_1 * coef_Jacob2_qt_syms11 * pinvG1_1 -
              coefs_tq1_1 * coef_Jacob2_qt_syms12 * pinvG2_1 - coefs_tq2_1 * coef_Jacob2_qt_syms11 * pinvG1_2 -
              coefs_tq1_1 * coef_Jacob2_qt_syms13 * pinvG3_1 - coefs_tq2_1 * coef_Jacob2_qt_syms12 * pinvG2_2 -
              coefs_tq3_1 * coef_Jacob2_qt_syms11 * pinvG1_3 - coefs_tq2_1 * coef_Jacob2_qt_syms13 * pinvG3_2 -
              coefs_tq3_1 * coef_Jacob2_qt_syms12 * pinvG2_3 - coefs_tq3_1 * coef_Jacob2_qt_syms13 * pinvG3_3;
    W(1, 1) = coef_Jacob2_qt_syms2 * (-1.0 * 2.0) - (coefs_tq1_2 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_1 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
              (coefs_tq1_2 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
              (coefs_tq1_2 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
              (coefs_tq2_2 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
              (coefs_tq3_2 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0;
    W(1, 2) = coef_Jacob2_qt_syms3 * (-1.0 * 2.0) - (coefs_tq1_3 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
              (coefs_tq2_3 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
              (coefs_tq3_3 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0;
    W(1, 3) = coef_Jacob2_qt_syms4 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_4 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_4 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
              (coefs_tq1_4 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_4 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_4 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
              (coefs_tq2_4 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_4 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
              (coefs_tq3_4 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0;
    W(1, 4) = coef_Jacob2_qt_syms2 * (-1.0 * 2.0) - (coefs_tq1_2 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_1 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
              (coefs_tq1_2 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
              (coefs_tq1_2 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
              (coefs_tq2_2 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
              (coefs_tq3_2 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0;
    W(1, 5) = coef_Jacob2_qt_syms5 * (-1.0 * 2.0) - (coefs_tq1_5 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_2 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
              (coefs_tq1_5 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_5 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_2 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
              (coefs_tq1_5 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_5 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_5 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_2 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
              (coefs_tq2_5 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_5 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq2_2 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
              (coefs_tq3_5 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq3_2 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0;
    W(1, 6) = coef_Jacob2_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
              (coefs_tq1_6 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
              (coefs_tq2_6 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
              (coefs_tq1_2 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
              (coefs_tq2_3 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
              (coefs_tq1_6 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
              (coefs_tq2_6 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
              (coefs_tq3_6 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
              (coefs_tq1_2 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
              (coefs_tq2_2 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
              (coefs_tq1_3 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
              (coefs_tq2_3 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
              (coefs_tq3_3 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
              (coefs_tq2_6 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
              (coefs_tq3_6 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
              (coefs_tq1_2 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
              (coefs_tq2_2 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
              (coefs_tq3_2 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
              (coefs_tq2_3 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
              (coefs_tq3_3 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
              (coefs_tq3_6 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
              (coefs_tq2_2 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
              (coefs_tq3_2 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
              (coefs_tq3_3 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
              (coefs_tq3_2 * coef_Jacob2_qt_syms30 * pinvG3_3) ;
    W(1, 7) = coef_Jacob2_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
              (coefs_tq1_4 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
              (coefs_tq1_7 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
              (coefs_tq2_7 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
              (coefs_tq1_2 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
              (coefs_tq1_4 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
              (coefs_tq2_4 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
              (coefs_tq1_7 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
              (coefs_tq2_7 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
              (coefs_tq3_7 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
              (coefs_tq1_2 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
              (coefs_tq2_2 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
              (coefs_tq1_4 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
              (coefs_tq2_4 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
              (coefs_tq3_4 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
              (coefs_tq2_7 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
              (coefs_tq3_7 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
              (coefs_tq1_2 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
              (coefs_tq2_2 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
              (coefs_tq3_2 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
              (coefs_tq2_4 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
              (coefs_tq3_4 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
              (coefs_tq3_7 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
              (coefs_tq2_2 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
              (coefs_tq3_2 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
              (coefs_tq3_4 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
              (coefs_tq3_2 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 8) = coef_Jacob2_qt_syms3 * (-1.0 * 2.0) - (coefs_tq1_3 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
              (coefs_tq2_3 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
              (coefs_tq3_3 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0;
    W(1, 9) = coef_Jacob2_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
              (coefs_tq1_6 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
              (coefs_tq2_6 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
              (coefs_tq1_2 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
              (coefs_tq2_3 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
              (coefs_tq1_6 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
              (coefs_tq2_6 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
              (coefs_tq3_6 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
              (coefs_tq1_2 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
              (coefs_tq2_2 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
              (coefs_tq1_3 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
              (coefs_tq2_3 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
              (coefs_tq3_3 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
              (coefs_tq2_6 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
              (coefs_tq3_6 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
              (coefs_tq1_2 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
              (coefs_tq2_2 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
              (coefs_tq3_2 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
              (coefs_tq2_3 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
              (coefs_tq3_3 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
              (coefs_tq3_6 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
              (coefs_tq2_2 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
              (coefs_tq3_2 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
              (coefs_tq3_3 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
              (coefs_tq3_2 * coef_Jacob2_qt_syms30 * pinvG3_3) ;
    W(1, 10) = coef_Jacob2_qt_syms8 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_3 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_3 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0;
    W(1, 11) = coef_Jacob2_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 12) = coef_Jacob2_qt_syms4 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq2_4 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0;
    W(1, 13) = coef_Jacob2_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 14) = coef_Jacob2_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 15) = coef_Jacob2_qt_syms10 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_4 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_4 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0;
    W(1, 16) = coef_Jacob2_qt_syms2 * (-1.0 * 2.0) - (coefs_tq1_2 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_1 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_2 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_2 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq2_2 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_2 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0;
    W(1, 17) = coef_Jacob2_qt_syms5 * (-1.0 * 2.0) - (coefs_tq1_5 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_2 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_2 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_2 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq2_2 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq3_2 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0;
    W(1, 18) = coef_Jacob2_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms30 * pinvG3_3) ;
    W(1, 19) = coef_Jacob2_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 20) = coef_Jacob2_qt_syms5 * (-1.0 * 2.0) - (coefs_tq1_5 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_2 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_2 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_2 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq2_2 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq3_2 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0;
    W(1, 21) = -coef_Jacob2_qt_syms15 - coefs_tq1_5 * coef_Jacob2_qt_syms21 * pinvG1_1 -
               coefs_tq1_5 * coef_Jacob2_qt_syms22 * pinvG2_1 - coefs_tq2_5 * coef_Jacob2_qt_syms21 * pinvG1_2 -
               coefs_tq1_5 * coef_Jacob2_qt_syms23 * pinvG3_1 - coefs_tq2_5 * coef_Jacob2_qt_syms22 * pinvG2_2 -
               coefs_tq3_5 * coef_Jacob2_qt_syms21 * pinvG1_3 - coefs_tq2_5 * coef_Jacob2_qt_syms23 * pinvG3_2 -
               coefs_tq3_5 * coef_Jacob2_qt_syms22 * pinvG2_3 - coefs_tq3_5 * coef_Jacob2_qt_syms23 * pinvG3_3;
    W(1, 22) = coef_Jacob2_qt_syms16 * (-1.0 * 2.0) - (coefs_tq1_6 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0;
    W(1, 23) = coef_Jacob2_qt_syms17 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0;
    W(1, 24) = coef_Jacob2_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms30 * pinvG3_3) ;
    W(1, 25) = coef_Jacob2_qt_syms16 * (-1.0 * 2.0) - (coefs_tq1_6 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0;
    W(1, 26) = coef_Jacob2_qt_syms18 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_6 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0;
    W(1, 27) = coef_Jacob2_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 28) = coef_Jacob2_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 29) = coef_Jacob2_qt_syms17 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0;
    W(1, 30) = coef_Jacob2_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 31) = coef_Jacob2_qt_syms20 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0;
    W(1, 32) = coef_Jacob2_qt_syms3 * (-1.0 * 2.0) - (coefs_tq1_3 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_3 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_3 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_3 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0;
    W(1, 33) = coef_Jacob2_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms30 * pinvG3_3) ;
    W(1, 34) = coef_Jacob2_qt_syms8 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_3 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_3 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0;
    W(1, 35) = coef_Jacob2_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 36) = coef_Jacob2_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms30 * pinvG3_3) ;
    W(1, 37) = coef_Jacob2_qt_syms16 * (-1.0 * 2.0) - (coefs_tq1_6 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0;
    W(1, 38) = coef_Jacob2_qt_syms18 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_6 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0;
    W(1, 39) = coef_Jacob2_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 40) = coef_Jacob2_qt_syms8 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_3 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_3 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0;
    W(1, 41) = coef_Jacob2_qt_syms18 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_6 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0;
    W(1, 42) = -coef_Jacob2_qt_syms25 - coefs_tq1_8 * coef_Jacob2_qt_syms28 * pinvG1_1 -
               coefs_tq1_8 * coef_Jacob2_qt_syms29 * pinvG2_1 - coefs_tq2_8 * coef_Jacob2_qt_syms28 * pinvG1_2 -
               coefs_tq1_8 * coef_Jacob2_qt_syms30 * pinvG3_1 - coefs_tq2_8 * coef_Jacob2_qt_syms29 * pinvG2_2 -
               coefs_tq3_8 * coef_Jacob2_qt_syms28 * pinvG1_3 - coefs_tq2_8 * coef_Jacob2_qt_syms30 * pinvG3_2 -
               coefs_tq3_8 * coef_Jacob2_qt_syms29 * pinvG2_3 - coefs_tq3_8 * coef_Jacob2_qt_syms30 * pinvG3_3;
    W(1, 43) = coef_Jacob2_qt_syms26 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0;
    W(1, 44) = coef_Jacob2_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 45) = coef_Jacob2_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 46) = coef_Jacob2_qt_syms26 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0;
    W(1, 47) = coef_Jacob2_qt_syms27 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0;
    W(1, 48) = coef_Jacob2_qt_syms4 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq2_4 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0;
    W(1, 49) = coef_Jacob2_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 50) = coef_Jacob2_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 51) = coef_Jacob2_qt_syms10 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_4 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_4 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0;
    W(1, 52) = coef_Jacob2_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 53) = coef_Jacob2_qt_syms17 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0;
    W(1, 54) = coef_Jacob2_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 55) = coef_Jacob2_qt_syms20 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0;
    W(1, 56) = coef_Jacob2_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob2_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob2_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 57) = coef_Jacob2_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob2_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob2_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob2_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob2_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob2_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob2_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob2_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob2_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob2_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob2_qt_syms35 * pinvG3_3) ;
    W(1, 58) = coef_Jacob2_qt_syms26 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq1_8 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0;
    W(1, 59) = coef_Jacob2_qt_syms27 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0;
    W(1, 60) = coef_Jacob2_qt_syms10 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_4 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_4 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms13 * pinvG3_3) * 2.0;
    W(1, 61) = coef_Jacob2_qt_syms20 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms23 * pinvG3_3) * 2.0;
    W(1, 62) = coef_Jacob2_qt_syms27 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob2_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob2_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob2_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob2_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob2_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob2_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob2_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob2_qt_syms30 * pinvG3_3) * 2.0;
    W(1, 63) = -coef_Jacob2_qt_syms32 - coefs_tq1_10 * coef_Jacob2_qt_syms33 * pinvG1_1 -
               coefs_tq1_10 * coef_Jacob2_qt_syms34 * pinvG2_1 - coefs_tq1_10 * coef_Jacob2_qt_syms35 * pinvG3_1 -
               coefs_tq2_10 * coef_Jacob2_qt_syms33 * pinvG1_2 - coefs_tq2_10 * coef_Jacob2_qt_syms34 * pinvG2_2 -
               coefs_tq2_10 * coef_Jacob2_qt_syms35 * pinvG3_2 - coefs_tq3_10 * coef_Jacob2_qt_syms33 * pinvG1_3 -
               coefs_tq3_10 * coef_Jacob2_qt_syms34 * pinvG2_3 - coefs_tq3_10 * coef_Jacob2_qt_syms35 * pinvG3_3;
    W(2, 0) = -coef_Jacob3_qt_syms1 - coefs_tq1_1 * coef_Jacob3_qt_syms11 * pinvG1_1 -
              coefs_tq1_1 * coef_Jacob3_qt_syms12 * pinvG2_1 - coefs_tq2_1 * coef_Jacob3_qt_syms11 * pinvG1_2 -
              coefs_tq1_1 * coef_Jacob3_qt_syms13 * pinvG3_1 - coefs_tq2_1 * coef_Jacob3_qt_syms12 * pinvG2_2 -
              coefs_tq3_1 * coef_Jacob3_qt_syms11 * pinvG1_3 - coefs_tq2_1 * coef_Jacob3_qt_syms13 * pinvG3_2 -
              coefs_tq3_1 * coef_Jacob3_qt_syms12 * pinvG2_3 - coefs_tq3_1 * coef_Jacob3_qt_syms13 * pinvG3_3;
    W(2, 1) = coef_Jacob3_qt_syms2 * (-1.0 * 2.0) - (coefs_tq1_2 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_1 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
              (coefs_tq1_2 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
              (coefs_tq1_2 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
              (coefs_tq2_2 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
              (coefs_tq3_2 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0;
    W(2, 2) = coef_Jacob3_qt_syms3 * (-1.0 * 2.0) - (coefs_tq1_3 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
              (coefs_tq2_3 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
              (coefs_tq3_3 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0;
    W(2, 3) = coef_Jacob3_qt_syms4 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_4 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_4 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
              (coefs_tq1_4 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_4 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_4 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
              (coefs_tq2_4 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_4 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
              (coefs_tq3_4 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0;
    W(2, 4) = coef_Jacob3_qt_syms2 * (-1.0 * 2.0) - (coefs_tq1_2 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_1 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
              (coefs_tq1_2 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
              (coefs_tq1_2 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
              (coefs_tq2_2 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
              (coefs_tq3_2 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0;
    W(2, 5) = coef_Jacob3_qt_syms5 * (-1.0 * 2.0) - (coefs_tq1_5 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_2 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
              (coefs_tq1_5 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_5 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_2 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
              (coefs_tq1_5 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_5 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_5 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_2 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
              (coefs_tq2_5 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_5 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq2_2 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
              (coefs_tq3_5 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq3_2 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0;
    W(2, 6) = coef_Jacob3_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
              (coefs_tq1_6 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
              (coefs_tq2_6 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
              (coefs_tq1_2 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
              (coefs_tq2_3 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
              (coefs_tq1_6 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
              (coefs_tq2_6 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
              (coefs_tq3_6 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
              (coefs_tq1_2 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
              (coefs_tq2_2 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
              (coefs_tq1_3 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
              (coefs_tq2_3 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
              (coefs_tq3_3 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
              (coefs_tq2_6 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
              (coefs_tq3_6 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
              (coefs_tq1_2 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
              (coefs_tq2_2 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
              (coefs_tq3_2 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
              (coefs_tq2_3 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
              (coefs_tq3_3 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
              (coefs_tq3_6 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
              (coefs_tq2_2 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
              (coefs_tq3_2 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
              (coefs_tq3_3 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
              (coefs_tq3_2 * coef_Jacob3_qt_syms30 * pinvG3_3) ;
    W(2, 7) = coef_Jacob3_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
              (coefs_tq1_4 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
              (coefs_tq1_7 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
              (coefs_tq2_7 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
              (coefs_tq1_2 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
              (coefs_tq1_4 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
              (coefs_tq2_4 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
              (coefs_tq1_7 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
              (coefs_tq2_7 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
              (coefs_tq3_7 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
              (coefs_tq1_2 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
              (coefs_tq2_2 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
              (coefs_tq1_4 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
              (coefs_tq2_4 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
              (coefs_tq3_4 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
              (coefs_tq2_7 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
              (coefs_tq3_7 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
              (coefs_tq1_2 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
              (coefs_tq2_2 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
              (coefs_tq3_2 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
              (coefs_tq2_4 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
              (coefs_tq3_4 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
              (coefs_tq3_7 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
              (coefs_tq2_2 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
              (coefs_tq3_2 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
              (coefs_tq3_4 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
              (coefs_tq3_2 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 8) = coef_Jacob3_qt_syms3 * (-1.0 * 2.0) - (coefs_tq1_3 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
              (coefs_tq2_3 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
              (coefs_tq3_3 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0;
    W(2, 9) = coef_Jacob3_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
              (coefs_tq1_6 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
              (coefs_tq2_6 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
              (coefs_tq1_2 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
              (coefs_tq2_3 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
              (coefs_tq1_6 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
              (coefs_tq2_6 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
              (coefs_tq3_6 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
              (coefs_tq1_2 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
              (coefs_tq2_2 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
              (coefs_tq1_3 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
              (coefs_tq2_3 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
              (coefs_tq3_3 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
              (coefs_tq2_6 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
              (coefs_tq3_6 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
              (coefs_tq1_2 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
              (coefs_tq2_2 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
              (coefs_tq3_2 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
              (coefs_tq2_3 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
              (coefs_tq3_3 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
              (coefs_tq3_6 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
              (coefs_tq2_2 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
              (coefs_tq3_2 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
              (coefs_tq3_3 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
              (coefs_tq3_2 * coef_Jacob3_qt_syms30 * pinvG3_3) ;
    W(2, 10) = coef_Jacob3_qt_syms8 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_3 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_3 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0;
    W(2, 11) = coef_Jacob3_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 12) = coef_Jacob3_qt_syms4 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq2_4 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0;
    W(2, 13) = coef_Jacob3_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 14) = coef_Jacob3_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 15) = coef_Jacob3_qt_syms10 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_4 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_4 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0;
    W(2, 16) = coef_Jacob3_qt_syms2 * (-1.0 * 2.0) - (coefs_tq1_2 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_1 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_2 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_2 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq2_2 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_2 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0;
    W(2, 17) = coef_Jacob3_qt_syms5 * (-1.0 * 2.0) - (coefs_tq1_5 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_2 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_2 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_2 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq2_2 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq3_2 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0;
    W(2, 18) = coef_Jacob3_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms30 * pinvG3_3) ;
    W(2, 19) = coef_Jacob3_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 20) = coef_Jacob3_qt_syms5 * (-1.0 * 2.0) - (coefs_tq1_5 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_2 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_2 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_2 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq2_2 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq3_2 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0;
    W(2, 21) = -coef_Jacob3_qt_syms15 - coefs_tq1_5 * coef_Jacob3_qt_syms21 * pinvG1_1 -
               coefs_tq1_5 * coef_Jacob3_qt_syms22 * pinvG2_1 - coefs_tq2_5 * coef_Jacob3_qt_syms21 * pinvG1_2 -
               coefs_tq1_5 * coef_Jacob3_qt_syms23 * pinvG3_1 - coefs_tq2_5 * coef_Jacob3_qt_syms22 * pinvG2_2 -
               coefs_tq3_5 * coef_Jacob3_qt_syms21 * pinvG1_3 - coefs_tq2_5 * coef_Jacob3_qt_syms23 * pinvG3_2 -
               coefs_tq3_5 * coef_Jacob3_qt_syms22 * pinvG2_3 - coefs_tq3_5 * coef_Jacob3_qt_syms23 * pinvG3_3;
    W(2, 22) = coef_Jacob3_qt_syms16 * (-1.0 * 2.0) - (coefs_tq1_6 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0;
    W(2, 23) = coef_Jacob3_qt_syms17 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0;
    W(2, 24) = coef_Jacob3_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms30 * pinvG3_3) ;
    W(2, 25) = coef_Jacob3_qt_syms16 * (-1.0 * 2.0) - (coefs_tq1_6 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0;
    W(2, 26) = coef_Jacob3_qt_syms18 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_6 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0;
    W(2, 27) = coef_Jacob3_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 28) = coef_Jacob3_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 29) = coef_Jacob3_qt_syms17 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0;
    W(2, 30) = coef_Jacob3_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 31) = coef_Jacob3_qt_syms20 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0;
    W(2, 32) = coef_Jacob3_qt_syms3 * (-1.0 * 2.0) - (coefs_tq1_3 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_3 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_3 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_3 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0;
    W(2, 33) = coef_Jacob3_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms30 * pinvG3_3) ;
    W(2, 34) = coef_Jacob3_qt_syms8 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_3 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_3 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0;
    W(2, 35) = coef_Jacob3_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 36) = coef_Jacob3_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms30 * pinvG3_3) ;
    W(2, 37) = coef_Jacob3_qt_syms16 * (-1.0 * 2.0) - (coefs_tq1_6 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0;
    W(2, 38) = coef_Jacob3_qt_syms18 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_6 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0;
    W(2, 39) = coef_Jacob3_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 40) = coef_Jacob3_qt_syms8 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_3 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_3 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0;
    W(2, 41) = coef_Jacob3_qt_syms18 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_6 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0;
    W(2, 42) = -coef_Jacob3_qt_syms25 - coefs_tq1_8 * coef_Jacob3_qt_syms28 * pinvG1_1 -
               coefs_tq1_8 * coef_Jacob3_qt_syms29 * pinvG2_1 - coefs_tq2_8 * coef_Jacob3_qt_syms28 * pinvG1_2 -
               coefs_tq1_8 * coef_Jacob3_qt_syms30 * pinvG3_1 - coefs_tq2_8 * coef_Jacob3_qt_syms29 * pinvG2_2 -
               coefs_tq3_8 * coef_Jacob3_qt_syms28 * pinvG1_3 - coefs_tq2_8 * coef_Jacob3_qt_syms30 * pinvG3_2 -
               coefs_tq3_8 * coef_Jacob3_qt_syms29 * pinvG2_3 - coefs_tq3_8 * coef_Jacob3_qt_syms30 * pinvG3_3;
    W(2, 43) = coef_Jacob3_qt_syms26 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0;
    W(2, 44) = coef_Jacob3_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 45) = coef_Jacob3_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 46) = coef_Jacob3_qt_syms26 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0;
    W(2, 47) = coef_Jacob3_qt_syms27 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0;
    W(2, 48) = coef_Jacob3_qt_syms4 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq2_4 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0;
    W(2, 49) = coef_Jacob3_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 50) = coef_Jacob3_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 51) = coef_Jacob3_qt_syms10 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_4 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_4 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0;
    W(2, 52) = coef_Jacob3_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 53) = coef_Jacob3_qt_syms17 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0;
    W(2, 54) = coef_Jacob3_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 55) = coef_Jacob3_qt_syms20 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0;
    W(2, 56) = coef_Jacob3_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob3_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob3_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 57) = coef_Jacob3_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob3_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob3_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob3_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob3_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob3_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob3_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob3_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob3_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob3_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob3_qt_syms35 * pinvG3_3) ;
    W(2, 58) = coef_Jacob3_qt_syms26 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq1_8 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0;
    W(2, 59) = coef_Jacob3_qt_syms27 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0;
    W(2, 60) = coef_Jacob3_qt_syms10 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_4 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_4 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms13 * pinvG3_3) * 2.0;
    W(2, 61) = coef_Jacob3_qt_syms20 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms23 * pinvG3_3) * 2.0;
    W(2, 62) = coef_Jacob3_qt_syms27 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob3_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob3_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob3_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob3_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob3_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob3_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob3_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob3_qt_syms30 * pinvG3_3) * 2.0;
    W(2, 63) = -coef_Jacob3_qt_syms32 - coefs_tq1_10 * coef_Jacob3_qt_syms33 * pinvG1_1 -
               coefs_tq1_10 * coef_Jacob3_qt_syms34 * pinvG2_1 - coefs_tq1_10 * coef_Jacob3_qt_syms35 * pinvG3_1 -
               coefs_tq2_10 * coef_Jacob3_qt_syms33 * pinvG1_2 - coefs_tq2_10 * coef_Jacob3_qt_syms34 * pinvG2_2 -
               coefs_tq2_10 * coef_Jacob3_qt_syms35 * pinvG3_2 - coefs_tq3_10 * coef_Jacob3_qt_syms33 * pinvG1_3 -
               coefs_tq3_10 * coef_Jacob3_qt_syms34 * pinvG2_3 - coefs_tq3_10 * coef_Jacob3_qt_syms35 * pinvG3_3;
    W(3, 0) = -coef_Jacob4_qt_syms1 - coefs_tq1_1 * coef_Jacob4_qt_syms11 * pinvG1_1 -
              coefs_tq1_1 * coef_Jacob4_qt_syms12 * pinvG2_1 - coefs_tq2_1 * coef_Jacob4_qt_syms11 * pinvG1_2 -
              coefs_tq1_1 * coef_Jacob4_qt_syms13 * pinvG3_1 - coefs_tq2_1 * coef_Jacob4_qt_syms12 * pinvG2_2 -
              coefs_tq3_1 * coef_Jacob4_qt_syms11 * pinvG1_3 - coefs_tq2_1 * coef_Jacob4_qt_syms13 * pinvG3_2 -
              coefs_tq3_1 * coef_Jacob4_qt_syms12 * pinvG2_3 - coefs_tq3_1 * coef_Jacob4_qt_syms13 * pinvG3_3;
    W(3, 1) = coef_Jacob4_qt_syms2 * (-1.0 * 2.0) - (coefs_tq1_2 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_1 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
              (coefs_tq1_2 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
              (coefs_tq1_2 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
              (coefs_tq2_2 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
              (coefs_tq3_2 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0;
    W(3, 2) = coef_Jacob4_qt_syms3 * (-1.0 * 2.0) - (coefs_tq1_3 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
              (coefs_tq2_3 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
              (coefs_tq3_3 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0;
    W(3, 3) = coef_Jacob4_qt_syms4 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_4 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_4 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
              (coefs_tq1_4 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_4 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_4 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
              (coefs_tq2_4 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_4 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
              (coefs_tq3_4 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0;
    W(3, 4) = coef_Jacob4_qt_syms2 * (-1.0 * 2.0) - (coefs_tq1_2 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_1 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
              (coefs_tq1_2 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
              (coefs_tq1_2 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
              (coefs_tq2_2 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
              (coefs_tq3_2 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0;
    W(3, 5) = coef_Jacob4_qt_syms5 * (-1.0 * 2.0) - (coefs_tq1_5 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_2 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
              (coefs_tq1_5 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_5 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_2 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
              (coefs_tq1_5 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_5 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_5 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_2 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
              (coefs_tq2_2 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
              (coefs_tq2_5 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_5 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq2_2 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
              (coefs_tq3_2 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
              (coefs_tq3_5 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq3_2 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0;
    W(3, 6) = coef_Jacob4_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
              (coefs_tq1_6 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
              (coefs_tq2_6 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
              (coefs_tq1_2 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
              (coefs_tq2_3 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
              (coefs_tq1_6 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
              (coefs_tq2_6 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
              (coefs_tq3_6 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
              (coefs_tq1_2 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
              (coefs_tq2_2 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
              (coefs_tq1_3 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
              (coefs_tq2_3 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
              (coefs_tq3_3 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
              (coefs_tq2_6 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
              (coefs_tq3_6 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
              (coefs_tq1_2 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
              (coefs_tq2_2 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
              (coefs_tq3_2 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
              (coefs_tq2_3 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
              (coefs_tq3_3 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
              (coefs_tq3_6 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
              (coefs_tq2_2 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
              (coefs_tq3_2 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
              (coefs_tq3_3 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
              (coefs_tq3_2 * coef_Jacob4_qt_syms30 * pinvG3_3) ;
    W(3, 7) = coef_Jacob4_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
              (coefs_tq1_4 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
              (coefs_tq1_7 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
              (coefs_tq2_7 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
              (coefs_tq1_2 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
              (coefs_tq1_4 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
              (coefs_tq2_4 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
              (coefs_tq1_7 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
              (coefs_tq2_7 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
              (coefs_tq3_7 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
              (coefs_tq1_2 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
              (coefs_tq2_2 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
              (coefs_tq1_4 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
              (coefs_tq2_4 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
              (coefs_tq3_4 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
              (coefs_tq2_7 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
              (coefs_tq3_7 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
              (coefs_tq1_2 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
              (coefs_tq2_2 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
              (coefs_tq3_2 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
              (coefs_tq2_4 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
              (coefs_tq3_4 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
              (coefs_tq3_7 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
              (coefs_tq2_2 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
              (coefs_tq3_2 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
              (coefs_tq3_4 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
              (coefs_tq3_2 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 8) = coef_Jacob4_qt_syms3 * (-1.0 * 2.0) - (coefs_tq1_3 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
              (coefs_tq1_1 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
              (coefs_tq1_3 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
              (coefs_tq2_3 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
              (coefs_tq2_3 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
              (coefs_tq3_3 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
              (coefs_tq1_1 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
              (coefs_tq2_1 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
              (coefs_tq3_3 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0 -
              (coefs_tq2_1 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
              (coefs_tq3_1 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
              (coefs_tq3_1 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0;
    W(3, 9) = coef_Jacob4_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
              (coefs_tq1_6 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
              (coefs_tq2_6 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
              (coefs_tq1_2 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
              (coefs_tq1_3 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
              (coefs_tq2_3 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
              (coefs_tq1_6 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
              (coefs_tq2_6 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
              (coefs_tq3_6 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
              (coefs_tq1_2 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
              (coefs_tq2_2 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
              (coefs_tq1_3 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
              (coefs_tq2_3 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
              (coefs_tq3_3 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
              (coefs_tq2_6 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
              (coefs_tq3_6 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
              (coefs_tq1_2 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
              (coefs_tq2_2 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
              (coefs_tq3_2 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
              (coefs_tq2_3 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
              (coefs_tq3_3 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
              (coefs_tq3_6 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
              (coefs_tq2_2 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
              (coefs_tq3_2 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
              (coefs_tq3_3 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
              (coefs_tq3_2 * coef_Jacob4_qt_syms30 * pinvG3_3) ;
    W(3, 10) = coef_Jacob4_qt_syms8 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_3 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_3 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0;
    W(3, 11) = coef_Jacob4_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 12) = coef_Jacob4_qt_syms4 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq2_4 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0;
    W(3, 13) = coef_Jacob4_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 14) = coef_Jacob4_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 15) = coef_Jacob4_qt_syms10 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_4 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_4 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0;
    W(3, 16) = coef_Jacob4_qt_syms2 * (-1.0 * 2.0) - (coefs_tq1_2 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_1 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_2 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_2 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq2_2 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_2 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0;
    W(3, 17) = coef_Jacob4_qt_syms5 * (-1.0 * 2.0) - (coefs_tq1_5 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_2 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_2 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_2 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq2_2 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq3_2 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0;
    W(3, 18) = coef_Jacob4_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms30 * pinvG3_3) ;
    W(3, 19) = coef_Jacob4_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 20) = coef_Jacob4_qt_syms5 * (-1.0 * 2.0) - (coefs_tq1_5 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_2 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_2 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_2 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_2 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq2_2 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_2 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq3_2 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0;
    W(3, 21) = -coef_Jacob4_qt_syms15 - coefs_tq1_5 * coef_Jacob4_qt_syms21 * pinvG1_1 -
               coefs_tq1_5 * coef_Jacob4_qt_syms22 * pinvG2_1 - coefs_tq2_5 * coef_Jacob4_qt_syms21 * pinvG1_2 -
               coefs_tq1_5 * coef_Jacob4_qt_syms23 * pinvG3_1 - coefs_tq2_5 * coef_Jacob4_qt_syms22 * pinvG2_2 -
               coefs_tq3_5 * coef_Jacob4_qt_syms21 * pinvG1_3 - coefs_tq2_5 * coef_Jacob4_qt_syms23 * pinvG3_2 -
               coefs_tq3_5 * coef_Jacob4_qt_syms22 * pinvG2_3 - coefs_tq3_5 * coef_Jacob4_qt_syms23 * pinvG3_3;
    W(3, 22) = coef_Jacob4_qt_syms16 * (-1.0 * 2.0) - (coefs_tq1_6 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0;
    W(3, 23) = coef_Jacob4_qt_syms17 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0;
    W(3, 24) = coef_Jacob4_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms30 * pinvG3_3) ;
    W(3, 25) = coef_Jacob4_qt_syms16 * (-1.0 * 2.0) - (coefs_tq1_6 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0;
    W(3, 26) = coef_Jacob4_qt_syms18 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_6 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0;
    W(3, 27) = coef_Jacob4_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 28) = coef_Jacob4_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 29) = coef_Jacob4_qt_syms17 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0;
    W(3, 30) = coef_Jacob4_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 31) = coef_Jacob4_qt_syms20 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0;
    W(3, 32) = coef_Jacob4_qt_syms3 * (-1.0 * 2.0) - (coefs_tq1_3 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_3 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_3 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_3 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0;
    W(3, 33) = coef_Jacob4_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms30 * pinvG3_3) ;
    W(3, 34) = coef_Jacob4_qt_syms8 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_3 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_3 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0;
    W(3, 35) = coef_Jacob4_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 36) = coef_Jacob4_qt_syms6 * (-1.0 ) - (coefs_tq1_6 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms30 * pinvG3_3) ;
    W(3, 37) = coef_Jacob4_qt_syms16 * (-1.0 * 2.0) - (coefs_tq1_6 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0;
    W(3, 38) = coef_Jacob4_qt_syms18 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_6 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0;
    W(3, 39) = coef_Jacob4_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 40) = coef_Jacob4_qt_syms8 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_3 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_3 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_3 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_3 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_3 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_3 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0;
    W(3, 41) = coef_Jacob4_qt_syms18 * (-1.0 * 2.0) - (coefs_tq1_8 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_6 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_6 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_6 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_6 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_6 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0;
    W(3, 42) = -coef_Jacob4_qt_syms25 - coefs_tq1_8 * coef_Jacob4_qt_syms28 * pinvG1_1 -
               coefs_tq1_8 * coef_Jacob4_qt_syms29 * pinvG2_1 - coefs_tq2_8 * coef_Jacob4_qt_syms28 * pinvG1_2 -
               coefs_tq1_8 * coef_Jacob4_qt_syms30 * pinvG3_1 - coefs_tq2_8 * coef_Jacob4_qt_syms29 * pinvG2_2 -
               coefs_tq3_8 * coef_Jacob4_qt_syms28 * pinvG1_3 - coefs_tq2_8 * coef_Jacob4_qt_syms30 * pinvG3_2 -
               coefs_tq3_8 * coef_Jacob4_qt_syms29 * pinvG2_3 - coefs_tq3_8 * coef_Jacob4_qt_syms30 * pinvG3_3;
    W(3, 43) = coef_Jacob4_qt_syms26 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0;
    W(3, 44) = coef_Jacob4_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 45) = coef_Jacob4_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 46) = coef_Jacob4_qt_syms26 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0;
    W(3, 47) = coef_Jacob4_qt_syms27 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0;
    W(3, 48) = coef_Jacob4_qt_syms4 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq1_1 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq2_4 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq1_1 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_1 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0 -
               (coefs_tq2_1 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_1 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_1 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0;
    W(3, 49) = coef_Jacob4_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 50) = coef_Jacob4_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 51) = coef_Jacob4_qt_syms10 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_4 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_4 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0;
    W(3, 52) = coef_Jacob4_qt_syms7 * (-1.0 ) - (coefs_tq1_7 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_2 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_2 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_2 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 53) = coef_Jacob4_qt_syms17 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq1_5 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq2_5 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0 -
               (coefs_tq3_5 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0;
    W(3, 54) = coef_Jacob4_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 55) = coef_Jacob4_qt_syms20 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0;
    W(3, 56) = coef_Jacob4_qt_syms9 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob4_qt_syms11 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms12 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms11 * pinvG1_2)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms13 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms12 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms11 * pinvG1_3)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms13 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms12 * pinvG2_3)  -
               (coefs_tq1_4 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_3 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms13 * pinvG3_3)  -
               (coefs_tq2_4 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_3 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_4 * coef_Jacob4_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_3 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 57) = coef_Jacob4_qt_syms19 * (-1.0 ) - (coefs_tq1_9 * coef_Jacob4_qt_syms21 * pinvG1_1)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms28 * pinvG1_1)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms33 * pinvG1_1)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms22 * pinvG2_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms21 * pinvG1_2)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms29 * pinvG2_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms28 * pinvG1_2)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms34 * pinvG2_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms33 * pinvG1_2)  -
               (coefs_tq1_9 * coef_Jacob4_qt_syms23 * pinvG3_1)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms22 * pinvG2_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms21 * pinvG1_3)  -
               (coefs_tq1_7 * coef_Jacob4_qt_syms30 * pinvG3_1)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms29 * pinvG2_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms28 * pinvG1_3)  -
               (coefs_tq1_6 * coef_Jacob4_qt_syms35 * pinvG3_1)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms34 * pinvG2_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms33 * pinvG1_3)  -
               (coefs_tq2_9 * coef_Jacob4_qt_syms23 * pinvG3_2)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms22 * pinvG2_3)  -
               (coefs_tq2_7 * coef_Jacob4_qt_syms30 * pinvG3_2)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms29 * pinvG2_3)  -
               (coefs_tq2_6 * coef_Jacob4_qt_syms35 * pinvG3_2)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms34 * pinvG2_3)  -
               (coefs_tq3_9 * coef_Jacob4_qt_syms23 * pinvG3_3)  -
               (coefs_tq3_7 * coef_Jacob4_qt_syms30 * pinvG3_3)  -
               (coefs_tq3_6 * coef_Jacob4_qt_syms35 * pinvG3_3) ;
    W(3, 58) = coef_Jacob4_qt_syms26 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq1_8 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq2_8 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0 -
               (coefs_tq3_8 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0;
    W(3, 59) = coef_Jacob4_qt_syms27 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0;
    W(3, 60) = coef_Jacob4_qt_syms10 * (-1.0 * 2.0) - (coefs_tq1_4 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_4 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_4 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_4 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_4 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_4 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_4 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms11 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms12 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms13 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms11 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms12 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms13 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms11 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms12 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms13 * pinvG3_3) * 2.0;
    W(3, 61) = coef_Jacob4_qt_syms20 * (-1.0 * 2.0) - (coefs_tq1_7 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_7 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_7 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_7 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_7 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms21 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms22 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms23 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms21 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms22 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms23 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms21 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms22 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms23 * pinvG3_3) * 2.0;
    W(3, 62) = coef_Jacob4_qt_syms27 * (-1.0 * 2.0) - (coefs_tq1_9 * coef_Jacob4_qt_syms33 * pinvG1_1) * 2.0 -
               (coefs_tq1_9 * coef_Jacob4_qt_syms34 * pinvG2_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms33 * pinvG1_2) * 2.0 -
               (coefs_tq1_9 * coef_Jacob4_qt_syms35 * pinvG3_1) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms34 * pinvG2_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms33 * pinvG1_3) * 2.0 -
               (coefs_tq2_9 * coef_Jacob4_qt_syms35 * pinvG3_2) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms34 * pinvG2_3) * 2.0 -
               (coefs_tq3_9 * coef_Jacob4_qt_syms35 * pinvG3_3) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms28 * pinvG1_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms29 * pinvG2_1) * 2.0 -
               (coefs_tq1_10 * coef_Jacob4_qt_syms30 * pinvG3_1) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms28 * pinvG1_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms29 * pinvG2_2) * 2.0 -
               (coefs_tq2_10 * coef_Jacob4_qt_syms30 * pinvG3_2) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms28 * pinvG1_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms29 * pinvG2_3) * 2.0 -
               (coefs_tq3_10 * coef_Jacob4_qt_syms30 * pinvG3_3) * 2.0;
    W(3, 63) = -coef_Jacob4_qt_syms32 - coefs_tq1_10 * coef_Jacob4_qt_syms33 * pinvG1_1 -
               coefs_tq1_10 * coef_Jacob4_qt_syms34 * pinvG2_1 - coefs_tq1_10 * coef_Jacob4_qt_syms35 * pinvG3_1 -
               coefs_tq2_10 * coef_Jacob4_qt_syms33 * pinvG1_2 - coefs_tq2_10 * coef_Jacob4_qt_syms34 * pinvG2_2 -
               coefs_tq2_10 * coef_Jacob4_qt_syms35 * pinvG3_2 - coefs_tq3_10 * coef_Jacob4_qt_syms33 * pinvG1_3 -
               coefs_tq3_10 * coef_Jacob4_qt_syms34 * pinvG2_3 - coefs_tq3_10 * coef_Jacob4_qt_syms35 * pinvG3_3;

    W = W * (1.0 / 6.0);


    Q(0, 0) = coef_Jacob1_qt_syms14 + coefs_tq1_11 * coef_Jacob1_qt_syms11 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob1_qt_syms12 * pinvG2_1 + coefs_tq1_11 * coef_Jacob1_qt_syms13 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob1_qt_syms11 * pinvG1_2 + coefs_tq2_11 * coef_Jacob1_qt_syms12 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob1_qt_syms13 * pinvG3_2 + coefs_tq3_11 * coef_Jacob1_qt_syms11 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob1_qt_syms12 * pinvG2_3 + coefs_tq3_11 * coef_Jacob1_qt_syms13 * pinvG3_3;
    Q(0, 1) = coef_Jacob1_qt_syms24 + coefs_tq1_11 * coef_Jacob1_qt_syms21 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob1_qt_syms22 * pinvG2_1 + coefs_tq1_11 * coef_Jacob1_qt_syms23 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob1_qt_syms21 * pinvG1_2 + coefs_tq2_11 * coef_Jacob1_qt_syms22 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob1_qt_syms23 * pinvG3_2 + coefs_tq3_11 * coef_Jacob1_qt_syms21 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob1_qt_syms22 * pinvG2_3 + coefs_tq3_11 * coef_Jacob1_qt_syms23 * pinvG3_3;
    Q(0, 2) = coef_Jacob1_qt_syms31 + coefs_tq1_11 * coef_Jacob1_qt_syms28 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob1_qt_syms29 * pinvG2_1 + coefs_tq1_11 * coef_Jacob1_qt_syms30 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob1_qt_syms28 * pinvG1_2 + coefs_tq2_11 * coef_Jacob1_qt_syms29 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob1_qt_syms30 * pinvG3_2 + coefs_tq3_11 * coef_Jacob1_qt_syms28 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob1_qt_syms29 * pinvG2_3 + coefs_tq3_11 * coef_Jacob1_qt_syms30 * pinvG3_3;
    Q(0, 3) = coef_Jacob1_qt_syms36 + coefs_tq1_11 * coef_Jacob1_qt_syms33 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob1_qt_syms34 * pinvG2_1 + coefs_tq1_11 * coef_Jacob1_qt_syms35 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob1_qt_syms33 * pinvG1_2 + coefs_tq2_11 * coef_Jacob1_qt_syms34 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob1_qt_syms35 * pinvG3_2 + coefs_tq3_11 * coef_Jacob1_qt_syms33 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob1_qt_syms34 * pinvG2_3 + coefs_tq3_11 * coef_Jacob1_qt_syms35 * pinvG3_3;
    Q(1, 0) = coef_Jacob2_qt_syms14 + coefs_tq1_11 * coef_Jacob2_qt_syms11 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob2_qt_syms12 * pinvG2_1 + coefs_tq1_11 * coef_Jacob2_qt_syms13 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob2_qt_syms11 * pinvG1_2 + coefs_tq2_11 * coef_Jacob2_qt_syms12 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob2_qt_syms13 * pinvG3_2 + coefs_tq3_11 * coef_Jacob2_qt_syms11 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob2_qt_syms12 * pinvG2_3 + coefs_tq3_11 * coef_Jacob2_qt_syms13 * pinvG3_3;
    Q(1, 1) = coef_Jacob2_qt_syms24 + coefs_tq1_11 * coef_Jacob2_qt_syms21 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob2_qt_syms22 * pinvG2_1 + coefs_tq1_11 * coef_Jacob2_qt_syms23 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob2_qt_syms21 * pinvG1_2 + coefs_tq2_11 * coef_Jacob2_qt_syms22 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob2_qt_syms23 * pinvG3_2 + coefs_tq3_11 * coef_Jacob2_qt_syms21 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob2_qt_syms22 * pinvG2_3 + coefs_tq3_11 * coef_Jacob2_qt_syms23 * pinvG3_3;
    Q(1, 2) = coef_Jacob2_qt_syms31 + coefs_tq1_11 * coef_Jacob2_qt_syms28 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob2_qt_syms29 * pinvG2_1 + coefs_tq1_11 * coef_Jacob2_qt_syms30 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob2_qt_syms28 * pinvG1_2 + coefs_tq2_11 * coef_Jacob2_qt_syms29 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob2_qt_syms30 * pinvG3_2 + coefs_tq3_11 * coef_Jacob2_qt_syms28 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob2_qt_syms29 * pinvG2_3 + coefs_tq3_11 * coef_Jacob2_qt_syms30 * pinvG3_3;
    Q(1, 3) = coef_Jacob2_qt_syms36 + coefs_tq1_11 * coef_Jacob2_qt_syms33 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob2_qt_syms34 * pinvG2_1 + coefs_tq1_11 * coef_Jacob2_qt_syms35 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob2_qt_syms33 * pinvG1_2 + coefs_tq2_11 * coef_Jacob2_qt_syms34 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob2_qt_syms35 * pinvG3_2 + coefs_tq3_11 * coef_Jacob2_qt_syms33 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob2_qt_syms34 * pinvG2_3 + coefs_tq3_11 * coef_Jacob2_qt_syms35 * pinvG3_3;
    Q(2, 0) = coef_Jacob3_qt_syms14 + coefs_tq1_11 * coef_Jacob3_qt_syms11 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob3_qt_syms12 * pinvG2_1 + coefs_tq1_11 * coef_Jacob3_qt_syms13 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob3_qt_syms11 * pinvG1_2 + coefs_tq2_11 * coef_Jacob3_qt_syms12 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob3_qt_syms13 * pinvG3_2 + coefs_tq3_11 * coef_Jacob3_qt_syms11 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob3_qt_syms12 * pinvG2_3 + coefs_tq3_11 * coef_Jacob3_qt_syms13 * pinvG3_3;
    Q(2, 1) = coef_Jacob3_qt_syms24 + coefs_tq1_11 * coef_Jacob3_qt_syms21 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob3_qt_syms22 * pinvG2_1 + coefs_tq1_11 * coef_Jacob3_qt_syms23 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob3_qt_syms21 * pinvG1_2 + coefs_tq2_11 * coef_Jacob3_qt_syms22 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob3_qt_syms23 * pinvG3_2 + coefs_tq3_11 * coef_Jacob3_qt_syms21 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob3_qt_syms22 * pinvG2_3 + coefs_tq3_11 * coef_Jacob3_qt_syms23 * pinvG3_3;
    Q(2, 2) = coef_Jacob3_qt_syms31 + coefs_tq1_11 * coef_Jacob3_qt_syms28 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob3_qt_syms29 * pinvG2_1 + coefs_tq1_11 * coef_Jacob3_qt_syms30 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob3_qt_syms28 * pinvG1_2 + coefs_tq2_11 * coef_Jacob3_qt_syms29 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob3_qt_syms30 * pinvG3_2 + coefs_tq3_11 * coef_Jacob3_qt_syms28 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob3_qt_syms29 * pinvG2_3 + coefs_tq3_11 * coef_Jacob3_qt_syms30 * pinvG3_3;
    Q(2, 3) = coef_Jacob3_qt_syms36 + coefs_tq1_11 * coef_Jacob3_qt_syms33 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob3_qt_syms34 * pinvG2_1 + coefs_tq1_11 * coef_Jacob3_qt_syms35 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob3_qt_syms33 * pinvG1_2 + coefs_tq2_11 * coef_Jacob3_qt_syms34 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob3_qt_syms35 * pinvG3_2 + coefs_tq3_11 * coef_Jacob3_qt_syms33 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob3_qt_syms34 * pinvG2_3 + coefs_tq3_11 * coef_Jacob3_qt_syms35 * pinvG3_3;
    Q(3, 0) = coef_Jacob4_qt_syms14 + coefs_tq1_11 * coef_Jacob4_qt_syms11 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob4_qt_syms12 * pinvG2_1 + coefs_tq1_11 * coef_Jacob4_qt_syms13 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob4_qt_syms11 * pinvG1_2 + coefs_tq2_11 * coef_Jacob4_qt_syms12 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob4_qt_syms13 * pinvG3_2 + coefs_tq3_11 * coef_Jacob4_qt_syms11 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob4_qt_syms12 * pinvG2_3 + coefs_tq3_11 * coef_Jacob4_qt_syms13 * pinvG3_3;
    Q(3, 1) = coef_Jacob4_qt_syms24 + coefs_tq1_11 * coef_Jacob4_qt_syms21 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob4_qt_syms22 * pinvG2_1 + coefs_tq1_11 * coef_Jacob4_qt_syms23 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob4_qt_syms21 * pinvG1_2 + coefs_tq2_11 * coef_Jacob4_qt_syms22 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob4_qt_syms23 * pinvG3_2 + coefs_tq3_11 * coef_Jacob4_qt_syms21 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob4_qt_syms22 * pinvG2_3 + coefs_tq3_11 * coef_Jacob4_qt_syms23 * pinvG3_3;
    Q(3, 2) = coef_Jacob4_qt_syms31 + coefs_tq1_11 * coef_Jacob4_qt_syms28 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob4_qt_syms29 * pinvG2_1 + coefs_tq1_11 * coef_Jacob4_qt_syms30 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob4_qt_syms28 * pinvG1_2 + coefs_tq2_11 * coef_Jacob4_qt_syms29 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob4_qt_syms30 * pinvG3_2 + coefs_tq3_11 * coef_Jacob4_qt_syms28 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob4_qt_syms29 * pinvG2_3 + coefs_tq3_11 * coef_Jacob4_qt_syms30 * pinvG3_3;
    Q(3, 3) = coef_Jacob4_qt_syms36 + coefs_tq1_11 * coef_Jacob4_qt_syms33 * pinvG1_1 +
              coefs_tq1_11 * coef_Jacob4_qt_syms34 * pinvG2_1 + coefs_tq1_11 * coef_Jacob4_qt_syms35 * pinvG3_1 +
              coefs_tq2_11 * coef_Jacob4_qt_syms33 * pinvG1_2 + coefs_tq2_11 * coef_Jacob4_qt_syms34 * pinvG2_2 +
              coefs_tq2_11 * coef_Jacob4_qt_syms35 * pinvG3_2 + coefs_tq3_11 * coef_Jacob4_qt_syms33 * pinvG1_3 +
              coefs_tq3_11 * coef_Jacob4_qt_syms34 * pinvG2_3 + coefs_tq3_11 * coef_Jacob4_qt_syms35 * pinvG3_3;
}


void D_pTop_func(Eigen::Matrix<double, 3, 28> &D,
                 Eigen::Matrix<double, 3, 9> &G,
                 Eigen::Vector3d &c,
                 const Eigen::Matrix<double, 4, 24> &coef_f_q_sym) {
    double coef_f0_q_sym1 = coef_f_q_sym(0, 0);
    double coef_f0_q_sym2 = coef_f_q_sym(0, 1);
    double coef_f0_q_sym3 = coef_f_q_sym(0, 2);
    double coef_f0_q_sym4 = coef_f_q_sym(0, 3);
    double coef_f0_q_sym5 = coef_f_q_sym(0, 4);
    double coef_f0_q_sym6 = coef_f_q_sym(0, 5);
    double coef_f0_q_sym7 = coef_f_q_sym(0, 6);
    double coef_f0_q_sym8 = coef_f_q_sym(0, 7);
    double coef_f0_q_sym9 = coef_f_q_sym(0, 8);
    double coef_f0_q_sym10 = coef_f_q_sym(0, 9);
    double coef_f0_q_sym11 = coef_f_q_sym(0, 10);
    double coef_f0_q_sym12 = coef_f_q_sym(0, 11);
    double coef_f0_q_sym13 = coef_f_q_sym(0, 12);
    double coef_f0_q_sym14 = coef_f_q_sym(0, 13);
    double coef_f0_q_sym15 = coef_f_q_sym(0, 14);
    double coef_f0_q_sym16 = coef_f_q_sym(0, 15);
    double coef_f0_q_sym17 = coef_f_q_sym(0, 16);
    double coef_f0_q_sym18 = coef_f_q_sym(0, 17);
    double coef_f0_q_sym19 = coef_f_q_sym(0, 18);
    double coef_f0_q_sym20 = coef_f_q_sym(0, 19);
    double coef_f0_q_sym21 = coef_f_q_sym(0, 20);
    double coef_f0_q_sym22 = coef_f_q_sym(0, 21);
    double coef_f0_q_sym23 = coef_f_q_sym(0, 22);
    double coef_f0_q_sym24 = coef_f_q_sym(0, 23);
    double coef_f1_q_sym1 = coef_f_q_sym(1, 0);
    double coef_f1_q_sym2 = coef_f_q_sym(1, 1);
    double coef_f1_q_sym3 = coef_f_q_sym(1, 2);
    double coef_f1_q_sym4 = coef_f_q_sym(1, 3);
    double coef_f1_q_sym5 = coef_f_q_sym(1, 4);
    double coef_f1_q_sym6 = coef_f_q_sym(1, 5);
    double coef_f1_q_sym7 = coef_f_q_sym(1, 6);
    double coef_f1_q_sym8 = coef_f_q_sym(1, 7);
    double coef_f1_q_sym9 = coef_f_q_sym(1, 8);
    double coef_f1_q_sym10 = coef_f_q_sym(1, 9);
    double coef_f1_q_sym11 = coef_f_q_sym(1, 10);
    double coef_f1_q_sym12 = coef_f_q_sym(1, 11);
    double coef_f1_q_sym13 = coef_f_q_sym(1, 12);
    double coef_f1_q_sym14 = coef_f_q_sym(1, 13);
    double coef_f1_q_sym15 = coef_f_q_sym(1, 14);
    double coef_f1_q_sym16 = coef_f_q_sym(1, 15);
    double coef_f1_q_sym17 = coef_f_q_sym(1, 16);
    double coef_f1_q_sym18 = coef_f_q_sym(1, 17);
    double coef_f1_q_sym19 = coef_f_q_sym(1, 18);
    double coef_f1_q_sym20 = coef_f_q_sym(1, 19);
    double coef_f1_q_sym21 = coef_f_q_sym(1, 20);
    double coef_f1_q_sym22 = coef_f_q_sym(1, 21);
    double coef_f1_q_sym23 = coef_f_q_sym(1, 22);
    double coef_f1_q_sym24 = coef_f_q_sym(1, 23);
    double coef_f2_q_sym1 = coef_f_q_sym(2, 0);
    double coef_f2_q_sym2 = coef_f_q_sym(2, 1);
    double coef_f2_q_sym3 = coef_f_q_sym(2, 2);
    double coef_f2_q_sym4 = coef_f_q_sym(2, 3);
    double coef_f2_q_sym5 = coef_f_q_sym(2, 4);
    double coef_f2_q_sym6 = coef_f_q_sym(2, 5);
    double coef_f2_q_sym7 = coef_f_q_sym(2, 6);
    double coef_f2_q_sym8 = coef_f_q_sym(2, 7);
    double coef_f2_q_sym9 = coef_f_q_sym(2, 8);
    double coef_f2_q_sym10 = coef_f_q_sym(2, 9);
    double coef_f2_q_sym11 = coef_f_q_sym(2, 10);
    double coef_f2_q_sym12 = coef_f_q_sym(2, 11);
    double coef_f2_q_sym13 = coef_f_q_sym(2, 12);
    double coef_f2_q_sym14 = coef_f_q_sym(2, 13);
    double coef_f2_q_sym15 = coef_f_q_sym(2, 14);
    double coef_f2_q_sym16 = coef_f_q_sym(2, 15);
    double coef_f2_q_sym17 = coef_f_q_sym(2, 16);
    double coef_f2_q_sym18 = coef_f_q_sym(2, 17);
    double coef_f2_q_sym19 = coef_f_q_sym(2, 18);
    double coef_f2_q_sym20 = coef_f_q_sym(2, 19);
    double coef_f2_q_sym21 = coef_f_q_sym(2, 20);
    double coef_f2_q_sym22 = coef_f_q_sym(2, 21);
    double coef_f2_q_sym23 = coef_f_q_sym(2, 22);
    double coef_f2_q_sym24 = coef_f_q_sym(2, 23);
    double coef_f3_q_sym1 = coef_f_q_sym(3, 0);
    double coef_f3_q_sym2 = coef_f_q_sym(3, 1);
    double coef_f3_q_sym3 = coef_f_q_sym(3, 2);
    double coef_f3_q_sym4 = coef_f_q_sym(3, 3);
    double coef_f3_q_sym5 = coef_f_q_sym(3, 4);
    double coef_f3_q_sym6 = coef_f_q_sym(3, 5);
    double coef_f3_q_sym7 = coef_f_q_sym(3, 6);
    double coef_f3_q_sym8 = coef_f_q_sym(3, 7);
    double coef_f3_q_sym9 = coef_f_q_sym(3, 8);
    double coef_f3_q_sym10 = coef_f_q_sym(3, 9);
    double coef_f3_q_sym11 = coef_f_q_sym(3, 10);
    double coef_f3_q_sym12 = coef_f_q_sym(3, 11);
    double coef_f3_q_sym13 = coef_f_q_sym(3, 12);
    double coef_f3_q_sym14 = coef_f_q_sym(3, 13);
    double coef_f3_q_sym15 = coef_f_q_sym(3, 14);
    double coef_f3_q_sym16 = coef_f_q_sym(3, 15);
    double coef_f3_q_sym17 = coef_f_q_sym(3, 16);
    double coef_f3_q_sym18 = coef_f_q_sym(3, 17);
    double coef_f3_q_sym19 = coef_f_q_sym(3, 18);
    double coef_f3_q_sym20 = coef_f_q_sym(3, 19);
    double coef_f3_q_sym21 = coef_f_q_sym(3, 20);
    double coef_f3_q_sym22 = coef_f_q_sym(3, 21);
    double coef_f3_q_sym23 = coef_f_q_sym(3, 22);
    double coef_f3_q_sym24 = coef_f_q_sym(3, 23);

    D(0, 0) = coef_f0_q_sym1 - coef_f1_q_sym2;
    D(0, 1) = -coef_f1_q_sym3;
    D(0, 2) = -coef_f1_q_sym4;
    D(0, 3) = coef_f0_q_sym5 - coef_f1_q_sym12;
    D(0, 4) = coef_f0_q_sym6 - coef_f1_q_sym13;
    D(0, 5) = coef_f0_q_sym7 - coef_f1_q_sym14;
    D(0, 6) = coef_f0_q_sym8 - coef_f1_q_sym15;
    D(0, 7) = coef_f0_q_sym9 - coef_f1_q_sym16;
    D(0, 8) = coef_f0_q_sym10 - coef_f1_q_sym17;
    D(0, 9) = -coef_f1_q_sym19;
    D(0, 10) = -coef_f1_q_sym20;
    D(0, 11) = -coef_f1_q_sym21;
    D(0, 12) = -coef_f1_q_sym23;
    D(0, 13) = -coef_f0_q_sym2 - coef_f1_q_sym1 + coef_f0_q_sym12 + coef_f1_q_sym5;
    D(0, 14) = -coef_f0_q_sym3 + coef_f0_q_sym13 + coef_f1_q_sym6;
    D(0, 15) = -coef_f0_q_sym4 + coef_f0_q_sym14 + coef_f1_q_sym7;
    D(0, 16) = -coef_f0_q_sym2 - coef_f1_q_sym1 * 2.0 + coef_f0_q_sym15 + coef_f1_q_sym5 + coef_f1_q_sym8;
    D(0, 17) = coef_f0_q_sym16 + coef_f1_q_sym9;
    D(0, 18) = -coef_f0_q_sym2 - coef_f1_q_sym1 * 2.0 + coef_f1_q_sym5 + coef_f0_q_sym17 + coef_f1_q_sym10;
    D(0, 19) = -coef_f0_q_sym3 + coef_f1_q_sym6 + coef_f0_q_sym19;
    D(0, 20) = -coef_f0_q_sym4 + coef_f1_q_sym7 + coef_f0_q_sym20;
    D(0, 21) = -coef_f0_q_sym3 + coef_f1_q_sym6 + coef_f0_q_sym21;
    D(0, 22) = -coef_f0_q_sym4 + coef_f1_q_sym7 + coef_f0_q_sym23;
    D(0, 23) = -coef_f1_q_sym1 + coef_f1_q_sym8;
    D(0, 24) = coef_f1_q_sym9;
    D(0, 25) = coef_f1_q_sym1 * -2.0 + coef_f1_q_sym8 + coef_f1_q_sym10;
    D(0, 26) = coef_f1_q_sym9;
    D(0, 27) = -coef_f1_q_sym1 + coef_f1_q_sym10;
    D(1, 0) = -coef_f2_q_sym2;
    D(1, 1) = coef_f0_q_sym1 - coef_f2_q_sym3;
    D(1, 2) = -coef_f2_q_sym4;
    D(1, 3) = -coef_f2_q_sym12;
    D(1, 4) = coef_f0_q_sym5 - coef_f2_q_sym13;
    D(1, 5) = -coef_f2_q_sym14;
    D(1, 6) = coef_f0_q_sym6 - coef_f2_q_sym15;
    D(1, 7) = coef_f0_q_sym7 - coef_f2_q_sym16;
    D(1, 8) = -coef_f2_q_sym17;
    D(1, 9) = coef_f0_q_sym8 - coef_f2_q_sym19;
    D(1, 10) = coef_f0_q_sym9 - coef_f2_q_sym20;
    D(1, 11) = coef_f0_q_sym10 - coef_f2_q_sym21;
    D(1, 12) = -coef_f2_q_sym23;
    D(1, 13) = -coef_f2_q_sym1 + coef_f2_q_sym5;
    D(1, 14) = -coef_f0_q_sym2 + coef_f0_q_sym12 + coef_f2_q_sym6;
    D(1, 15) = coef_f2_q_sym7;
    D(1, 16) = -coef_f0_q_sym3 + coef_f0_q_sym13 - coef_f2_q_sym1 * 2.0 + coef_f2_q_sym5 + coef_f2_q_sym8;
    D(1, 17) = -coef_f0_q_sym4 + coef_f0_q_sym14 + coef_f2_q_sym9;
    D(1, 18) = coef_f2_q_sym1 * -2.0 + coef_f2_q_sym5 + coef_f2_q_sym10;
    D(1, 19) = -coef_f0_q_sym2 + coef_f0_q_sym15 + coef_f2_q_sym6;
    D(1, 20) = coef_f0_q_sym16 + coef_f2_q_sym7;
    D(1, 21) = -coef_f0_q_sym2 + coef_f0_q_sym17 + coef_f2_q_sym6;
    D(1, 22) = coef_f2_q_sym7;
    D(1, 23) = -coef_f0_q_sym3 + coef_f0_q_sym19 - coef_f2_q_sym1 + coef_f2_q_sym8;
    D(1, 24) = -coef_f0_q_sym4 + coef_f0_q_sym20 + coef_f2_q_sym9;
    D(1, 25) = -coef_f0_q_sym3 + coef_f0_q_sym21 - coef_f2_q_sym1 * 2.0 + coef_f2_q_sym8 + coef_f2_q_sym10;
    D(1, 26) = -coef_f0_q_sym4 + coef_f0_q_sym23 + coef_f2_q_sym9;
    D(1, 27) = -coef_f2_q_sym1 + coef_f2_q_sym10;
    D(2, 0) = -coef_f3_q_sym2;
    D(2, 1) = -coef_f3_q_sym3;
    D(2, 2) = coef_f0_q_sym1 - coef_f3_q_sym4;
    D(2, 3) = -coef_f3_q_sym12;
    D(2, 4) = -coef_f3_q_sym13;
    D(2, 5) = coef_f0_q_sym5 - coef_f3_q_sym14;
    D(2, 6) = -coef_f3_q_sym15;
    D(2, 7) = coef_f0_q_sym6 - coef_f3_q_sym16;
    D(2, 8) = coef_f0_q_sym7 - coef_f3_q_sym17;
    D(2, 9) = -coef_f3_q_sym19;
    D(2, 10) = coef_f0_q_sym8 - coef_f3_q_sym20;
    D(2, 11) = coef_f0_q_sym9 - coef_f3_q_sym21;
    D(2, 12) = coef_f0_q_sym10 - coef_f3_q_sym23;
    D(2, 13) = -coef_f3_q_sym1 + coef_f3_q_sym5;
    D(2, 14) = coef_f3_q_sym6;
    D(2, 15) = -coef_f0_q_sym2 + coef_f0_q_sym12 + coef_f3_q_sym7;
    D(2, 16) = coef_f3_q_sym1 * -2.0 + coef_f3_q_sym5 + coef_f3_q_sym8;
    D(2, 17) = -coef_f0_q_sym3 + coef_f0_q_sym13 + coef_f3_q_sym9;
    D(2, 18) = -coef_f0_q_sym4 + coef_f0_q_sym14 - coef_f3_q_sym1 * 2.0 + coef_f3_q_sym5 + coef_f3_q_sym10;
    D(2, 19) = coef_f3_q_sym6;
    D(2, 20) = -coef_f0_q_sym2 + coef_f0_q_sym15 + coef_f3_q_sym7;
    D(2, 21) = coef_f0_q_sym16 + coef_f3_q_sym6;
    D(2, 22) = -coef_f0_q_sym2 + coef_f0_q_sym17 + coef_f3_q_sym7;
    D(2, 23) = -coef_f3_q_sym1 + coef_f3_q_sym8;
    D(2, 24) = -coef_f0_q_sym3 + coef_f0_q_sym19 + coef_f3_q_sym9;
    D(2, 25) = -coef_f0_q_sym4 + coef_f0_q_sym20 - coef_f3_q_sym1 * 2.0 + coef_f3_q_sym8 + coef_f3_q_sym10;
    D(2, 26) = -coef_f0_q_sym3 + coef_f0_q_sym21 + coef_f3_q_sym9;
    D(2, 27) = -coef_f0_q_sym4 + coef_f0_q_sym23 - coef_f3_q_sym1 + coef_f3_q_sym10;

    G(0, 0) = coef_f0_q_sym11 - coef_f1_q_sym18;
    G(0, 1) = -coef_f1_q_sym22;
    G(0, 2) = -coef_f1_q_sym24;
    G(0, 3) = coef_f0_q_sym2 + coef_f1_q_sym1 * 2.0 - coef_f1_q_sym5 + coef_f0_q_sym18 + coef_f1_q_sym11;
    G(0, 4) = coef_f0_q_sym3 - coef_f1_q_sym6 + coef_f0_q_sym22;
    G(0, 5) = coef_f0_q_sym4 - coef_f1_q_sym7 + coef_f0_q_sym24;
    G(0, 6) = coef_f1_q_sym1 * 2.0 - coef_f1_q_sym8 + coef_f1_q_sym11;
    G(0, 7) = -coef_f1_q_sym9;
    G(0, 8) = coef_f1_q_sym1 * 2.0 - coef_f1_q_sym10 + coef_f1_q_sym11;
    G(1, 0) = -coef_f2_q_sym18;
    G(1, 1) = coef_f0_q_sym11 - coef_f2_q_sym22;
    G(1, 2) = -coef_f2_q_sym24;
    G(1, 3) = coef_f2_q_sym1 * 2.0 - coef_f2_q_sym5 + coef_f2_q_sym11;
    G(1, 4) = coef_f0_q_sym2 + coef_f0_q_sym18 - coef_f2_q_sym6;
    G(1, 5) = -coef_f2_q_sym7;
    G(1, 6) = coef_f0_q_sym3 + coef_f2_q_sym1 * 2.0 + coef_f0_q_sym22 - coef_f2_q_sym8 + coef_f2_q_sym11;
    G(1, 7) = coef_f0_q_sym4 + coef_f0_q_sym24 - coef_f2_q_sym9;
    G(1, 8) = coef_f2_q_sym1 * 2.0 - coef_f2_q_sym10 + coef_f2_q_sym11;
    G(2, 0) = -coef_f3_q_sym18;
    G(2, 1) = -coef_f3_q_sym22;
    G(2, 2) = coef_f0_q_sym11 - coef_f3_q_sym24;
    G(2, 3) = coef_f3_q_sym1 * 2.0 - coef_f3_q_sym5 + coef_f3_q_sym11;
    G(2, 4) = -coef_f3_q_sym6;
    G(2, 5) = coef_f0_q_sym2 + coef_f0_q_sym18 - coef_f3_q_sym7;
    G(2, 6) = coef_f3_q_sym1 * 2.0 - coef_f3_q_sym8 + coef_f3_q_sym11;
    G(2, 7) = coef_f0_q_sym3 + coef_f0_q_sym22 - coef_f3_q_sym9;
    G(2, 8) = coef_f0_q_sym4 + coef_f0_q_sym24 + coef_f3_q_sym1 * 2.0 - coef_f3_q_sym10 + coef_f3_q_sym11;

    c(0) = -coef_f1_q_sym1 - coef_f1_q_sym11;
    c(1) = -coef_f2_q_sym1 - coef_f2_q_sym11;
    c(2) = -coef_f3_q_sym1 - coef_f3_q_sym11;
}

void pTop_WQD(Eigen::Matrix<double, 4, 64> &W,
              Eigen::Matrix<double, 4, 4> &Q,
              Eigen::Matrix<double, 3, 28> &D,
              Eigen::Matrix<double, 3, 9> &GG,
              Eigen::Vector3d &c,
              Eigen::Matrix<double, 4, 24> &coef_f_q_sym,
              Eigen::Matrix<double, 1, 85> &coef_J_pure,
              Eigen::Matrix<double, 3, 11> &coefs_tq,
              Eigen::Matrix<double, 3, 3> &pinvG,
              const std::vector<Eigen::Vector3d> &rr,
              const std::vector<Eigen::Vector3d> &bb,
              const std::vector<Eigen::Vector3d> &nv) {

    int len = rr.size();

    double factor = 1.0 / ((double) len);

    {
        std::vector<Eigen::Matrix<double, 1, 85> > coef_J_pures(len, Eigen::Matrix<double, 1, 85>::Zero());
        std::vector<Eigen::Matrix<double, 9, 1> > pack(len);

        for (int i = 0; i < len; ++i) {
            pack[i] << rr[i], bb[i], nv[i];
        }

        for (int ii = 0; ii < len; ++ii) {
            mixed_pTop_func(coef_J_pures[ii], pack[ii]);
        }

        coef_J_pure = factor * std::accumulate(coef_J_pures.begin(), coef_J_pures.end(),
                                               Eigen::Matrix<double, 1, 85>::Zero().eval());
    }

    Eigen::Matrix<double, 1, 11> coeftq1;
    Eigen::Matrix<double, 1, 11> coeftq2;
    Eigen::Matrix<double, 1, 11> coeftq3;
    Eigen::Matrix<double, 1, 36> coef_Jacob1_qt;
    Eigen::Matrix<double, 1, 36> coef_Jacob2_qt;
    Eigen::Matrix<double, 1, 36> coef_Jacob3_qt;
    Eigen::Matrix<double, 1, 36> coef_Jacob4_qt;
    Eigen::Matrix3d G;
    coeftq1.setZero();
    coeftq2.setZero();
    coeftq3.setZero();
    coef_Jacob1_qt.setZero();
    coef_Jacob2_qt.setZero();
    coef_Jacob3_qt.setZero();
    coef_Jacob4_qt.setZero();
    G.setZero();
    mixed2_pTop_func(G, coeftq1, coeftq2, coeftq3, coef_Jacob1_qt, coef_Jacob2_qt, coef_Jacob3_qt, coef_Jacob4_qt,
                     coef_J_pure);


    auto svd = G.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
    const auto &singularValues = svd.singularValues();
    Eigen::Matrix3d singularValuesInv;
    singularValuesInv.setZero();
    double pinvtoler = 1.e-9; // choose your tolerance wisely
    for (int i = 0; i < 3; ++i) {
        if (singularValues(i) > pinvtoler)
            singularValuesInv(i, i) = 1.0 / singularValues(i);
        else
            singularValuesInv(i, i) = 0.0;
    }

    pinvG = svd.matrixV() * singularValuesInv * svd.matrixU().transpose();

    coefs_tq << coeftq1, coeftq2, coeftq3;

    Eigen::Matrix<double, 4, 36> coef_Jacob_qt_syms;
    coef_Jacob_qt_syms << coef_Jacob1_qt, coef_Jacob2_qt, coef_Jacob3_qt, coef_Jacob4_qt;
    mixed3_pTop_func(coef_f_q_sym, W, Q, pinvG, coefs_tq, coef_Jacob_qt_syms);
    D_pTop_func(D, GG, c, coef_f_q_sym);
}
