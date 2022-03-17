// LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems,
//          It also gives highly accurate uncertainty description of the solutions.
//
// Authors: Jin Wu and Ming Liu
// Affiliation: Hong Kong University of Science and Technology (HKUST)
// Emails: jin_wu_uestc@hotmail.com; eelium@ust.hk
// Reference: Wu, J., et al. (2022) Quadratic Pose Estimation Problems: 
//                                  Globally Optimal Solutions, 
//                                  Solvability/Observability Analysis,
//                                  and Uncertainty Description.
//                                  IEEE Transactions on Robotics.
//                                  https://doi.org/10.1109/TRO.2022.3155880
//
//
// misc_pTop_funcs.cpp: Functions for computation of Point-to-Plane registration problems


#include "misc_pTop_funcs.h"



Eigen::Vector3d t_pTop_func(const Eigen::MatrixXd& pinvG, const Eigen::MatrixXd& coefs_tq, const Eigen::Quaterniond& q)
{
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
    double q0 = q.x();
    double q1 = q.y();
    double q2 = q.z();
    double q3 = q.w();

    Eigen::Vector3d t;
    t(0) = (q0*q0)*(coefs_tq1_1*pinvG1_1+coefs_tq2_1*pinvG1_2+coefs_tq3_1*pinvG1_3)+(q1*q1)*(coefs_tq1_5*pinvG1_1+coefs_tq2_5*pinvG1_2+coefs_tq3_5*pinvG1_3)+(q2*q2)*(coefs_tq1_8*pinvG1_1+coefs_tq2_8*pinvG1_2+coefs_tq3_8*pinvG1_3)+(q3*q3)*(coefs_tq1_10*pinvG1_1+coefs_tq2_10*pinvG1_2+coefs_tq3_10*pinvG1_3)+coefs_tq1_11*pinvG1_1+coefs_tq2_11*pinvG1_2+coefs_tq3_11*pinvG1_3+q0*q1*(coefs_tq1_2*pinvG1_1+coefs_tq2_2*pinvG1_2+coefs_tq3_2*pinvG1_3)+q0*q2*(coefs_tq1_3*pinvG1_1+coefs_tq2_3*pinvG1_2+coefs_tq3_3*pinvG1_3)+q0*q3*(coefs_tq1_4*pinvG1_1+coefs_tq2_4*pinvG1_2+coefs_tq3_4*pinvG1_3)+q1*q2*(coefs_tq1_6*pinvG1_1+coefs_tq2_6*pinvG1_2+coefs_tq3_6*pinvG1_3)+q1*q3*(coefs_tq1_7*pinvG1_1+coefs_tq2_7*pinvG1_2+coefs_tq3_7*pinvG1_3)+q2*q3*(coefs_tq1_9*pinvG1_1+coefs_tq2_9*pinvG1_2+coefs_tq3_9*pinvG1_3);
    t(1) = (q0*q0)*(coefs_tq1_1*pinvG2_1+coefs_tq2_1*pinvG2_2+coefs_tq3_1*pinvG2_3)+(q1*q1)*(coefs_tq1_5*pinvG2_1+coefs_tq2_5*pinvG2_2+coefs_tq3_5*pinvG2_3)+(q2*q2)*(coefs_tq1_8*pinvG2_1+coefs_tq2_8*pinvG2_2+coefs_tq3_8*pinvG2_3)+(q3*q3)*(coefs_tq1_10*pinvG2_1+coefs_tq2_10*pinvG2_2+coefs_tq3_10*pinvG2_3)+coefs_tq1_11*pinvG2_1+coefs_tq2_11*pinvG2_2+coefs_tq3_11*pinvG2_3+q0*q1*(coefs_tq1_2*pinvG2_1+coefs_tq2_2*pinvG2_2+coefs_tq3_2*pinvG2_3)+q0*q2*(coefs_tq1_3*pinvG2_1+coefs_tq2_3*pinvG2_2+coefs_tq3_3*pinvG2_3)+q0*q3*(coefs_tq1_4*pinvG2_1+coefs_tq2_4*pinvG2_2+coefs_tq3_4*pinvG2_3)+q1*q2*(coefs_tq1_6*pinvG2_1+coefs_tq2_6*pinvG2_2+coefs_tq3_6*pinvG2_3)+q1*q3*(coefs_tq1_7*pinvG2_1+coefs_tq2_7*pinvG2_2+coefs_tq3_7*pinvG2_3)+q2*q3*(coefs_tq1_9*pinvG2_1+coefs_tq2_9*pinvG2_2+coefs_tq3_9*pinvG2_3);
    t(2) = (q0*q0)*(coefs_tq1_1*pinvG3_1+coefs_tq2_1*pinvG3_2+coefs_tq3_1*pinvG3_3)+(q1*q1)*(coefs_tq1_5*pinvG3_1+coefs_tq2_5*pinvG3_2+coefs_tq3_5*pinvG3_3)+(q2*q2)*(coefs_tq1_8*pinvG3_1+coefs_tq2_8*pinvG3_2+coefs_tq3_8*pinvG3_3)+(q3*q3)*(coefs_tq1_10*pinvG3_1+coefs_tq2_10*pinvG3_2+coefs_tq3_10*pinvG3_3)+coefs_tq1_11*pinvG3_1+coefs_tq2_11*pinvG3_2+coefs_tq3_11*pinvG3_3+q0*q1*(coefs_tq1_2*pinvG3_1+coefs_tq2_2*pinvG3_2+coefs_tq3_2*pinvG3_3)+q0*q2*(coefs_tq1_3*pinvG3_1+coefs_tq2_3*pinvG3_2+coefs_tq3_3*pinvG3_3)+q0*q3*(coefs_tq1_4*pinvG3_1+coefs_tq2_4*pinvG3_2+coefs_tq3_4*pinvG3_3)+q1*q2*(coefs_tq1_6*pinvG3_1+coefs_tq2_6*pinvG3_2+coefs_tq3_6*pinvG3_3)+q1*q3*(coefs_tq1_7*pinvG3_1+coefs_tq2_7*pinvG3_2+coefs_tq3_7*pinvG3_3)+q2*q3*(coefs_tq1_9*pinvG3_1+coefs_tq2_9*pinvG3_2+coefs_tq3_9*pinvG3_3);
    return t;
}



std::vector<double> mon_J_pure_pTop_func(const Eigen::Quaterniond& q, const Eigen::Vector3d& t)
{
    double q0 = q.x();
    double q1 = q.y();
    double q2 = q.z();
    double q3 = q.w();
    double t1 = t(0);
    double t2 = t(1);
    double t3 = t(2);

    double symobj[1][85];

    std::memset(symobj, 0, 1 * 85 * sizeof(double));
    symobj[0][0] = q0*q0*q0*q0;
    symobj[0][1] = (q0*q0*q0)*q1;
    symobj[0][2] = (q0*q0*q0)*q2;
    symobj[0][3] = (q0*q0*q0)*q3;
    symobj[0][4] = (q0*q0)*(q1*q1);
    symobj[0][5] = (q0*q0)*q1*q2;
    symobj[0][6] = (q0*q0)*q1*q3;
    symobj[0][7] = (q0*q0)*(q2*q2);
    symobj[0][8] = (q0*q0)*q2*q3;
    symobj[0][9] = (q0*q0)*(q3*q3);
    symobj[0][10] = (q0*q0)*t1;
    symobj[0][11] = (q0*q0)*t2;
    symobj[0][12] = (q0*q0)*t3;
    symobj[0][13] = q0*q0;
    symobj[0][14] = q0*(q1*q1*q1);
    symobj[0][15] = q0*(q1*q1)*q2;
    symobj[0][16] = q0*(q1*q1)*q3;
    symobj[0][17] = q0*q1*(q2*q2);
    symobj[0][18] = q0*q1*q2*q3;
    symobj[0][19] = q0*q1*(q3*q3);
    symobj[0][20] = q0*q1*t1;
    symobj[0][21] = q0*q1*t2;
    symobj[0][22] = q0*q1*t3;
    symobj[0][23] = q0*q1;
    symobj[0][24] = q0*(q2*q2*q2);
    symobj[0][25] = q0*(q2*q2)*q3;
    symobj[0][26] = q0*q2*(q3*q3);
    symobj[0][27] = q0*q2*t1;
    symobj[0][28] = q0*q2*t2;
    symobj[0][29] = q0*q2*t3;
    symobj[0][30] = q0*q2;
    symobj[0][31] = q0*(q3*q3*q3);
    symobj[0][32] = q0*q3*t1;
    symobj[0][33] = q0*q3*t2;
    symobj[0][34] = q0*q3*t3;
    symobj[0][35] = q0*q3;
    symobj[0][36] = q1*q1*q1*q1;
    symobj[0][37] = (q1*q1*q1)*q2;
    symobj[0][38] = (q1*q1*q1)*q3;
    symobj[0][39] = (q1*q1)*(q2*q2);
    symobj[0][40] = (q1*q1)*q2*q3;
    symobj[0][41] = (q1*q1)*(q3*q3);
    symobj[0][42] = (q1*q1)*t1;
    symobj[0][43] = (q1*q1)*t2;
    symobj[0][44] = (q1*q1)*t3;
    symobj[0][45] = q1*q1;
    symobj[0][46] = q1*(q2*q2*q2);
    symobj[0][47] = q1*(q2*q2)*q3;
    symobj[0][48] = q1*q2*(q3*q3);
    symobj[0][49] = q1*q2*t1;
    symobj[0][50] = q1*q2*t2;
    symobj[0][51] = q1*q2*t3;
    symobj[0][52] = q1*q2;
    symobj[0][53] = q1*(q3*q3*q3);
    symobj[0][54] = q1*q3*t1;
    symobj[0][55] = q1*q3*t2;
    symobj[0][56] = q1*q3*t3;
    symobj[0][57] = q1*q3;
    symobj[0][58] = q2*q2*q2*q2;
    symobj[0][59] = (q2*q2*q2)*q3;
    symobj[0][60] = (q2*q2)*(q3*q3);
    symobj[0][61] = (q2*q2)*t1;
    symobj[0][62] = (q2*q2)*t2;
    symobj[0][63] = (q2*q2)*t3;
    symobj[0][64] = q2*q2;
    symobj[0][65] = q2*(q3*q3*q3);
    symobj[0][66] = q2*q3*t1;
    symobj[0][67] = q2*q3*t2;
    symobj[0][68] = q2*q3*t3;
    symobj[0][69] = q2*q3;
    symobj[0][70] = q3*q3*q3*q3;
    symobj[0][71] = (q3*q3)*t1;
    symobj[0][72] = (q3*q3)*t2;
    symobj[0][73] = (q3*q3)*t3;
    symobj[0][74] = q3*q3;
    symobj[0][75] = t1*t1;
    symobj[0][76] = t1*t2;
    symobj[0][77] = t1*t3;
    symobj[0][78] = t1;
    symobj[0][79] = t2*t2;
    symobj[0][80] = t2*t3;
    symobj[0][81] = t2;
    symobj[0][82] = t3*t3;
    symobj[0][83] = t3;
    symobj[0][84] = 1.0;

    std::vector<double> tmp;
    for(int i = 0; i < 1; ++i)
        for(int j = 0; j < 85; ++j)
            tmp.push_back(symobj[i][j]);
    return tmp;
}



void eq_Jacob_pTop_func(Eigen::Matrix<double, 4, 1>& eq,
                                         Eigen::Matrix<double, 4, 4>& Jacob,
                                         const Eigen::MatrixXd& coef_f_q_sym,
                                         const Eigen::Vector4d& q)
{
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
    double q0 = q.x();
    double q1 = q.y();
    double q2 = q.z();
    double q3 = q.w();

    double symobj[4][4];

    std::memset(symobj, 0, 4 * 1 * sizeof(double));
    symobj[0][0] = -coef_f1_q_sym1*(q0*q0*q0*q0)+coef_f0_q_sym12*(q1*q1*q1*q1)+coef_f0_q_sym18*(q1*q1)-coef_f1_q_sym11*(q0*q0)+coef_f0_q_sym1*(q0*q0*q0)*q1+coef_f0_q_sym5*q0*(q1*q1*q1)-coef_f1_q_sym2*(q0*q0*q0)*q1-coef_f1_q_sym3*(q0*q0*q0)*q2+coef_f0_q_sym13*(q1*q1*q1)*q2-coef_f1_q_sym4*(q0*q0*q0)*q3+coef_f0_q_sym14*(q1*q1*q1)*q3+coef_f0_q_sym19*q1*(q2*q2*q2)+coef_f0_q_sym23*q1*(q3*q3*q3)-coef_f1_q_sym12*q0*(q1*q1*q1)-coef_f1_q_sym19*q0*(q2*q2*q2)-coef_f1_q_sym23*q0*(q3*q3*q3)+coef_f0_q_sym2*(q0*q0)*(q1*q1)-coef_f1_q_sym5*(q0*q0)*(q1*q1)+coef_f0_q_sym15*(q1*q1)*(q2*q2)-coef_f1_q_sym8*(q0*q0)*(q2*q2)+coef_f0_q_sym17*(q1*q1)*(q3*q3)-coef_f1_q_sym10*(q0*q0)*(q3*q3)+coef_f0_q_sym11*q0*q1+coef_f0_q_sym22*q1*q2+coef_f0_q_sym24*q1*q3-coef_f1_q_sym18*q0*q1-coef_f1_q_sym22*q0*q2-coef_f1_q_sym24*q0*q3+coef_f0_q_sym3*(q0*q0)*q1*q2+coef_f0_q_sym4*(q0*q0)*q1*q3+coef_f0_q_sym6*q0*(q1*q1)*q2+coef_f0_q_sym7*q0*(q1*q1)*q3+coef_f0_q_sym8*q0*q1*(q2*q2)+coef_f0_q_sym10*q0*q1*(q3*q3)-coef_f1_q_sym6*(q0*q0)*q1*q2-coef_f1_q_sym7*(q0*q0)*q1*q3+coef_f0_q_sym16*(q1*q1)*q2*q3-coef_f1_q_sym9*(q0*q0)*q2*q3+coef_f0_q_sym20*q1*(q2*q2)*q3+coef_f0_q_sym21*q1*q2*(q3*q3)-coef_f1_q_sym13*q0*(q1*q1)*q2-coef_f1_q_sym14*q0*(q1*q1)*q3-coef_f1_q_sym15*q0*q1*(q2*q2)-coef_f1_q_sym17*q0*q1*(q3*q3)-coef_f1_q_sym20*q0*(q2*q2)*q3-coef_f1_q_sym21*q0*q2*(q3*q3)+coef_f0_q_sym9*q0*q1*q2*q3-coef_f1_q_sym16*q0*q1*q2*q3;
    symobj[1][0] = coef_f0_q_sym19*(q2*q2*q2*q2)-coef_f2_q_sym1*(q0*q0*q0*q0)+coef_f0_q_sym22*(q2*q2)-coef_f2_q_sym11*(q0*q0)+coef_f0_q_sym1*(q0*q0*q0)*q2+coef_f0_q_sym8*q0*(q2*q2*q2)+coef_f0_q_sym12*(q1*q1*q1)*q2+coef_f0_q_sym15*q1*(q2*q2*q2)-coef_f2_q_sym2*(q0*q0*q0)*q1+coef_f0_q_sym20*(q2*q2*q2)*q3-coef_f2_q_sym3*(q0*q0*q0)*q2-coef_f2_q_sym4*(q0*q0*q0)*q3+coef_f0_q_sym23*q2*(q3*q3*q3)-coef_f2_q_sym12*q0*(q1*q1*q1)-coef_f2_q_sym19*q0*(q2*q2*q2)-coef_f2_q_sym23*q0*(q3*q3*q3)+coef_f0_q_sym3*(q0*q0)*(q2*q2)+coef_f0_q_sym13*(q1*q1)*(q2*q2)+coef_f0_q_sym21*(q2*q2)*(q3*q3)-coef_f2_q_sym5*(q0*q0)*(q1*q1)-coef_f2_q_sym8*(q0*q0)*(q2*q2)-coef_f2_q_sym10*(q0*q0)*(q3*q3)+coef_f0_q_sym11*q0*q2+coef_f0_q_sym18*q1*q2+coef_f0_q_sym24*q2*q3-coef_f2_q_sym18*q0*q1-coef_f2_q_sym22*q0*q2-coef_f2_q_sym24*q0*q3+coef_f0_q_sym2*(q0*q0)*q1*q2+coef_f0_q_sym5*q0*(q1*q1)*q2+coef_f0_q_sym4*(q0*q0)*q2*q3+coef_f0_q_sym6*q0*q1*(q2*q2)+coef_f0_q_sym9*q0*(q2*q2)*q3+coef_f0_q_sym10*q0*q2*(q3*q3)+coef_f0_q_sym14*(q1*q1)*q2*q3+coef_f0_q_sym16*q1*(q2*q2)*q3+coef_f0_q_sym17*q1*q2*(q3*q3)-coef_f2_q_sym6*(q0*q0)*q1*q2-coef_f2_q_sym7*(q0*q0)*q1*q3-coef_f2_q_sym9*(q0*q0)*q2*q3-coef_f2_q_sym13*q0*(q1*q1)*q2-coef_f2_q_sym14*q0*(q1*q1)*q3-coef_f2_q_sym15*q0*q1*(q2*q2)-coef_f2_q_sym17*q0*q1*(q3*q3)-coef_f2_q_sym20*q0*(q2*q2)*q3-coef_f2_q_sym21*q0*q2*(q3*q3)+coef_f0_q_sym7*q0*q1*q2*q3-coef_f2_q_sym16*q0*q1*q2*q3;
    symobj[2][0] = coef_f0_q_sym24*(q3*q3)+coef_f0_q_sym23*(q3*q3*q3*q3)-coef_f3_q_sym1*(q0*q0*q0*q0)-coef_f3_q_sym11*(q0*q0)+coef_f0_q_sym1*(q0*q0*q0)*q3+coef_f0_q_sym10*q0*(q3*q3*q3)+coef_f0_q_sym12*(q1*q1*q1)*q3+coef_f0_q_sym17*q1*(q3*q3*q3)+coef_f0_q_sym19*(q2*q2*q2)*q3+coef_f0_q_sym21*q2*(q3*q3*q3)-coef_f3_q_sym2*(q0*q0*q0)*q1-coef_f3_q_sym3*(q0*q0*q0)*q2-coef_f3_q_sym4*(q0*q0*q0)*q3-coef_f3_q_sym12*q0*(q1*q1*q1)-coef_f3_q_sym19*q0*(q2*q2*q2)-coef_f3_q_sym23*q0*(q3*q3*q3)+coef_f0_q_sym4*(q0*q0)*(q3*q3)+coef_f0_q_sym14*(q1*q1)*(q3*q3)+coef_f0_q_sym20*(q2*q2)*(q3*q3)-coef_f3_q_sym5*(q0*q0)*(q1*q1)-coef_f3_q_sym8*(q0*q0)*(q2*q2)-coef_f3_q_sym10*(q0*q0)*(q3*q3)+coef_f0_q_sym11*q0*q3+coef_f0_q_sym18*q1*q3+coef_f0_q_sym22*q2*q3-coef_f3_q_sym18*q0*q1-coef_f3_q_sym22*q0*q2-coef_f3_q_sym24*q0*q3+coef_f0_q_sym2*(q0*q0)*q1*q3+coef_f0_q_sym3*(q0*q0)*q2*q3+coef_f0_q_sym5*q0*(q1*q1)*q3+coef_f0_q_sym7*q0*q1*(q3*q3)+coef_f0_q_sym8*q0*(q2*q2)*q3+coef_f0_q_sym9*q0*q2*(q3*q3)+coef_f0_q_sym13*(q1*q1)*q2*q3+coef_f0_q_sym15*q1*(q2*q2)*q3+coef_f0_q_sym16*q1*q2*(q3*q3)-coef_f3_q_sym6*(q0*q0)*q1*q2-coef_f3_q_sym7*(q0*q0)*q1*q3-coef_f3_q_sym9*(q0*q0)*q2*q3-coef_f3_q_sym13*q0*(q1*q1)*q2-coef_f3_q_sym14*q0*(q1*q1)*q3-coef_f3_q_sym15*q0*q1*(q2*q2)-coef_f3_q_sym17*q0*q1*(q3*q3)-coef_f3_q_sym20*q0*(q2*q2)*q3-coef_f3_q_sym21*q0*q2*(q3*q3)+coef_f0_q_sym6*q0*q1*q2*q3-coef_f3_q_sym16*q0*q1*q2*q3;
    symobj[3][0] = q0*q0+q1*q1+q2*q2+q3*q3-1.0;

    for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 1; ++j)
            eq(i, j) = symobj[i][j];

    std::memset(symobj, 0, 4 * 4 * sizeof(double));
    symobj[0][0] = coef_f0_q_sym11*q1-coef_f1_q_sym11*q0*2.0-coef_f1_q_sym18*q1-coef_f1_q_sym22*q2-coef_f1_q_sym24*q3+coef_f0_q_sym5*(q1*q1*q1)-coef_f1_q_sym1*(q0*q0*q0)*4.0-coef_f1_q_sym12*(q1*q1*q1)-coef_f1_q_sym19*(q2*q2*q2)-coef_f1_q_sym23*(q3*q3*q3)+coef_f0_q_sym1*(q0*q0)*q1*3.0+coef_f0_q_sym2*q0*(q1*q1)*2.0+coef_f0_q_sym6*(q1*q1)*q2+coef_f0_q_sym7*(q1*q1)*q3+coef_f0_q_sym8*q1*(q2*q2)-coef_f1_q_sym2*(q0*q0)*q1*3.0+coef_f0_q_sym10*q1*(q3*q3)-coef_f1_q_sym3*(q0*q0)*q2*3.0-coef_f1_q_sym5*q0*(q1*q1)*2.0-coef_f1_q_sym4*(q0*q0)*q3*3.0-coef_f1_q_sym8*q0*(q2*q2)*2.0-coef_f1_q_sym10*q0*(q3*q3)*2.0-coef_f1_q_sym13*(q1*q1)*q2-coef_f1_q_sym14*(q1*q1)*q3-coef_f1_q_sym15*q1*(q2*q2)-coef_f1_q_sym17*q1*(q3*q3)-coef_f1_q_sym20*(q2*q2)*q3-coef_f1_q_sym21*q2*(q3*q3)+coef_f0_q_sym3*q0*q1*q2*2.0+coef_f0_q_sym4*q0*q1*q3*2.0+coef_f0_q_sym9*q1*q2*q3-coef_f1_q_sym6*q0*q1*q2*2.0-coef_f1_q_sym7*q0*q1*q3*2.0-coef_f1_q_sym9*q0*q2*q3*2.0-coef_f1_q_sym16*q1*q2*q3;
    symobj[0][1] = coef_f0_q_sym11*q0+coef_f0_q_sym18*q1*2.0+coef_f0_q_sym22*q2+coef_f0_q_sym24*q3-coef_f1_q_sym18*q0+coef_f0_q_sym1*(q0*q0*q0)-coef_f1_q_sym2*(q0*q0*q0)+coef_f0_q_sym12*(q1*q1*q1)*4.0+coef_f0_q_sym19*(q2*q2*q2)+coef_f0_q_sym23*(q3*q3*q3)+coef_f0_q_sym2*(q0*q0)*q1*2.0+coef_f0_q_sym3*(q0*q0)*q2+coef_f0_q_sym5*q0*(q1*q1)*3.0+coef_f0_q_sym4*(q0*q0)*q3+coef_f0_q_sym8*q0*(q2*q2)+coef_f0_q_sym10*q0*(q3*q3)+coef_f0_q_sym13*(q1*q1)*q2*3.0-coef_f1_q_sym5*(q0*q0)*q1*2.0+coef_f0_q_sym14*(q1*q1)*q3*3.0+coef_f0_q_sym15*q1*(q2*q2)*2.0-coef_f1_q_sym6*(q0*q0)*q2-coef_f1_q_sym7*(q0*q0)*q3+coef_f0_q_sym17*q1*(q3*q3)*2.0+coef_f0_q_sym20*(q2*q2)*q3+coef_f0_q_sym21*q2*(q3*q3)-coef_f1_q_sym12*q0*(q1*q1)*3.0-coef_f1_q_sym15*q0*(q2*q2)-coef_f1_q_sym17*q0*(q3*q3)+coef_f0_q_sym6*q0*q1*q2*2.0+coef_f0_q_sym7*q0*q1*q3*2.0+coef_f0_q_sym9*q0*q2*q3+coef_f0_q_sym16*q1*q2*q3*2.0-coef_f1_q_sym13*q0*q1*q2*2.0-coef_f1_q_sym14*q0*q1*q3*2.0-coef_f1_q_sym16*q0*q2*q3;
    symobj[0][2] = coef_f0_q_sym22*q1-coef_f1_q_sym22*q0-coef_f1_q_sym3*(q0*q0*q0)+coef_f0_q_sym13*(q1*q1*q1)+coef_f0_q_sym3*(q0*q0)*q1+coef_f0_q_sym6*q0*(q1*q1)-coef_f1_q_sym6*(q0*q0)*q1+coef_f0_q_sym15*(q1*q1)*q2*2.0+coef_f0_q_sym16*(q1*q1)*q3-coef_f1_q_sym8*(q0*q0)*q2*2.0+coef_f0_q_sym19*q1*(q2*q2)*3.0-coef_f1_q_sym9*(q0*q0)*q3+coef_f0_q_sym21*q1*(q3*q3)-coef_f1_q_sym13*q0*(q1*q1)-coef_f1_q_sym19*q0*(q2*q2)*3.0-coef_f1_q_sym21*q0*(q3*q3)+coef_f0_q_sym8*q0*q1*q2*2.0+coef_f0_q_sym9*q0*q1*q3+coef_f0_q_sym20*q1*q2*q3*2.0-coef_f1_q_sym15*q0*q1*q2*2.0-coef_f1_q_sym16*q0*q1*q3-coef_f1_q_sym20*q0*q2*q3*2.0;
    symobj[0][3] = coef_f0_q_sym24*q1-coef_f1_q_sym24*q0-coef_f1_q_sym4*(q0*q0*q0)+coef_f0_q_sym14*(q1*q1*q1)+coef_f0_q_sym4*(q0*q0)*q1+coef_f0_q_sym7*q0*(q1*q1)-coef_f1_q_sym7*(q0*q0)*q1+coef_f0_q_sym16*(q1*q1)*q2+coef_f0_q_sym17*(q1*q1)*q3*2.0-coef_f1_q_sym9*(q0*q0)*q2+coef_f0_q_sym20*q1*(q2*q2)+coef_f0_q_sym23*q1*(q3*q3)*3.0-coef_f1_q_sym10*(q0*q0)*q3*2.0-coef_f1_q_sym14*q0*(q1*q1)-coef_f1_q_sym20*q0*(q2*q2)-coef_f1_q_sym23*q0*(q3*q3)*3.0+coef_f0_q_sym9*q0*q1*q2+coef_f0_q_sym10*q0*q1*q3*2.0+coef_f0_q_sym21*q1*q2*q3*2.0-coef_f1_q_sym16*q0*q1*q2-coef_f1_q_sym17*q0*q1*q3*2.0-coef_f1_q_sym21*q0*q2*q3*2.0;
    symobj[1][0] = coef_f0_q_sym11*q2-coef_f2_q_sym11*q0*2.0-coef_f2_q_sym18*q1-coef_f2_q_sym22*q2-coef_f2_q_sym24*q3+coef_f0_q_sym8*(q2*q2*q2)-coef_f2_q_sym1*(q0*q0*q0)*4.0-coef_f2_q_sym12*(q1*q1*q1)-coef_f2_q_sym19*(q2*q2*q2)-coef_f2_q_sym23*(q3*q3*q3)+coef_f0_q_sym1*(q0*q0)*q2*3.0+coef_f0_q_sym3*q0*(q2*q2)*2.0+coef_f0_q_sym5*(q1*q1)*q2+coef_f0_q_sym6*q1*(q2*q2)+coef_f0_q_sym9*(q2*q2)*q3+coef_f0_q_sym10*q2*(q3*q3)-coef_f2_q_sym2*(q0*q0)*q1*3.0-coef_f2_q_sym3*(q0*q0)*q2*3.0-coef_f2_q_sym5*q0*(q1*q1)*2.0-coef_f2_q_sym4*(q0*q0)*q3*3.0-coef_f2_q_sym8*q0*(q2*q2)*2.0-coef_f2_q_sym10*q0*(q3*q3)*2.0-coef_f2_q_sym13*(q1*q1)*q2-coef_f2_q_sym14*(q1*q1)*q3-coef_f2_q_sym15*q1*(q2*q2)-coef_f2_q_sym17*q1*(q3*q3)-coef_f2_q_sym20*(q2*q2)*q3-coef_f2_q_sym21*q2*(q3*q3)+coef_f0_q_sym2*q0*q1*q2*2.0+coef_f0_q_sym4*q0*q2*q3*2.0+coef_f0_q_sym7*q1*q2*q3-coef_f2_q_sym6*q0*q1*q2*2.0-coef_f2_q_sym7*q0*q1*q3*2.0-coef_f2_q_sym9*q0*q2*q3*2.0-coef_f2_q_sym16*q1*q2*q3;
    symobj[1][1] = coef_f0_q_sym18*q2-coef_f2_q_sym18*q0+coef_f0_q_sym15*(q2*q2*q2)-coef_f2_q_sym2*(q0*q0*q0)+coef_f0_q_sym2*(q0*q0)*q2+coef_f0_q_sym6*q0*(q2*q2)+coef_f0_q_sym12*(q1*q1)*q2*3.0+coef_f0_q_sym13*q1*(q2*q2)*2.0+coef_f0_q_sym16*(q2*q2)*q3+coef_f0_q_sym17*q2*(q3*q3)-coef_f2_q_sym5*(q0*q0)*q1*2.0-coef_f2_q_sym6*(q0*q0)*q2-coef_f2_q_sym7*(q0*q0)*q3-coef_f2_q_sym12*q0*(q1*q1)*3.0-coef_f2_q_sym15*q0*(q2*q2)-coef_f2_q_sym17*q0*(q3*q3)+coef_f0_q_sym5*q0*q1*q2*2.0+coef_f0_q_sym7*q0*q2*q3+coef_f0_q_sym14*q1*q2*q3*2.0-coef_f2_q_sym13*q0*q1*q2*2.0-coef_f2_q_sym14*q0*q1*q3*2.0-coef_f2_q_sym16*q0*q2*q3;
    symobj[1][2] = coef_f0_q_sym11*q0+coef_f0_q_sym18*q1+coef_f0_q_sym22*q2*2.0+coef_f0_q_sym24*q3-coef_f2_q_sym22*q0+coef_f0_q_sym1*(q0*q0*q0)+coef_f0_q_sym12*(q1*q1*q1)+coef_f0_q_sym19*(q2*q2*q2)*4.0-coef_f2_q_sym3*(q0*q0*q0)+coef_f0_q_sym23*(q3*q3*q3)+coef_f0_q_sym2*(q0*q0)*q1+coef_f0_q_sym3*(q0*q0)*q2*2.0+coef_f0_q_sym5*q0*(q1*q1)+coef_f0_q_sym4*(q0*q0)*q3+coef_f0_q_sym8*q0*(q2*q2)*3.0+coef_f0_q_sym10*q0*(q3*q3)+coef_f0_q_sym13*(q1*q1)*q2*2.0+coef_f0_q_sym14*(q1*q1)*q3+coef_f0_q_sym15*q1*(q2*q2)*3.0+coef_f0_q_sym17*q1*(q3*q3)+coef_f0_q_sym20*(q2*q2)*q3*3.0+coef_f0_q_sym21*q2*(q3*q3)*2.0-coef_f2_q_sym6*(q0*q0)*q1-coef_f2_q_sym8*(q0*q0)*q2*2.0-coef_f2_q_sym9*(q0*q0)*q3-coef_f2_q_sym13*q0*(q1*q1)-coef_f2_q_sym19*q0*(q2*q2)*3.0-coef_f2_q_sym21*q0*(q3*q3)+coef_f0_q_sym6*q0*q1*q2*2.0+coef_f0_q_sym7*q0*q1*q3+coef_f0_q_sym9*q0*q2*q3*2.0+coef_f0_q_sym16*q1*q2*q3*2.0-coef_f2_q_sym15*q0*q1*q2*2.0-coef_f2_q_sym16*q0*q1*q3-coef_f2_q_sym20*q0*q2*q3*2.0;
    symobj[1][3] = coef_f0_q_sym24*q2-coef_f2_q_sym24*q0+coef_f0_q_sym20*(q2*q2*q2)-coef_f2_q_sym4*(q0*q0*q0)+coef_f0_q_sym4*(q0*q0)*q2+coef_f0_q_sym9*q0*(q2*q2)+coef_f0_q_sym14*(q1*q1)*q2+coef_f0_q_sym16*q1*(q2*q2)+coef_f0_q_sym21*(q2*q2)*q3*2.0+coef_f0_q_sym23*q2*(q3*q3)*3.0-coef_f2_q_sym7*(q0*q0)*q1-coef_f2_q_sym9*(q0*q0)*q2-coef_f2_q_sym10*(q0*q0)*q3*2.0-coef_f2_q_sym14*q0*(q1*q1)-coef_f2_q_sym20*q0*(q2*q2)-coef_f2_q_sym23*q0*(q3*q3)*3.0+coef_f0_q_sym7*q0*q1*q2+coef_f0_q_sym10*q0*q2*q3*2.0+coef_f0_q_sym17*q1*q2*q3*2.0-coef_f2_q_sym16*q0*q1*q2-coef_f2_q_sym17*q0*q1*q3*2.0-coef_f2_q_sym21*q0*q2*q3*2.0;
    symobj[2][0] = coef_f0_q_sym11*q3-coef_f3_q_sym11*q0*2.0-coef_f3_q_sym18*q1-coef_f3_q_sym22*q2-coef_f3_q_sym24*q3+coef_f0_q_sym10*(q3*q3*q3)-coef_f3_q_sym1*(q0*q0*q0)*4.0-coef_f3_q_sym12*(q1*q1*q1)-coef_f3_q_sym19*(q2*q2*q2)-coef_f3_q_sym23*(q3*q3*q3)+coef_f0_q_sym1*(q0*q0)*q3*3.0+coef_f0_q_sym4*q0*(q3*q3)*2.0+coef_f0_q_sym5*(q1*q1)*q3+coef_f0_q_sym7*q1*(q3*q3)+coef_f0_q_sym8*(q2*q2)*q3+coef_f0_q_sym9*q2*(q3*q3)-coef_f3_q_sym2*(q0*q0)*q1*3.0-coef_f3_q_sym3*(q0*q0)*q2*3.0-coef_f3_q_sym5*q0*(q1*q1)*2.0-coef_f3_q_sym4*(q0*q0)*q3*3.0-coef_f3_q_sym8*q0*(q2*q2)*2.0-coef_f3_q_sym10*q0*(q3*q3)*2.0-coef_f3_q_sym13*(q1*q1)*q2-coef_f3_q_sym14*(q1*q1)*q3-coef_f3_q_sym15*q1*(q2*q2)-coef_f3_q_sym17*q1*(q3*q3)-coef_f3_q_sym20*(q2*q2)*q3-coef_f3_q_sym21*q2*(q3*q3)+coef_f0_q_sym2*q0*q1*q3*2.0+coef_f0_q_sym3*q0*q2*q3*2.0+coef_f0_q_sym6*q1*q2*q3-coef_f3_q_sym6*q0*q1*q2*2.0-coef_f3_q_sym7*q0*q1*q3*2.0-coef_f3_q_sym9*q0*q2*q3*2.0-coef_f3_q_sym16*q1*q2*q3;
    symobj[2][1] = coef_f0_q_sym18*q3-coef_f3_q_sym18*q0+coef_f0_q_sym17*(q3*q3*q3)-coef_f3_q_sym2*(q0*q0*q0)+coef_f0_q_sym2*(q0*q0)*q3+coef_f0_q_sym7*q0*(q3*q3)+coef_f0_q_sym12*(q1*q1)*q3*3.0+coef_f0_q_sym14*q1*(q3*q3)*2.0+coef_f0_q_sym15*(q2*q2)*q3+coef_f0_q_sym16*q2*(q3*q3)-coef_f3_q_sym5*(q0*q0)*q1*2.0-coef_f3_q_sym6*(q0*q0)*q2-coef_f3_q_sym7*(q0*q0)*q3-coef_f3_q_sym12*q0*(q1*q1)*3.0-coef_f3_q_sym15*q0*(q2*q2)-coef_f3_q_sym17*q0*(q3*q3)+coef_f0_q_sym5*q0*q1*q3*2.0+coef_f0_q_sym6*q0*q2*q3+coef_f0_q_sym13*q1*q2*q3*2.0-coef_f3_q_sym13*q0*q1*q2*2.0-coef_f3_q_sym14*q0*q1*q3*2.0-coef_f3_q_sym16*q0*q2*q3;
    symobj[2][2] = coef_f0_q_sym22*q3-coef_f3_q_sym22*q0+coef_f0_q_sym21*(q3*q3*q3)-coef_f3_q_sym3*(q0*q0*q0)+coef_f0_q_sym3*(q0*q0)*q3+coef_f0_q_sym9*q0*(q3*q3)+coef_f0_q_sym13*(q1*q1)*q3+coef_f0_q_sym16*q1*(q3*q3)+coef_f0_q_sym19*(q2*q2)*q3*3.0+coef_f0_q_sym20*q2*(q3*q3)*2.0-coef_f3_q_sym6*(q0*q0)*q1-coef_f3_q_sym8*(q0*q0)*q2*2.0-coef_f3_q_sym9*(q0*q0)*q3-coef_f3_q_sym13*q0*(q1*q1)-coef_f3_q_sym19*q0*(q2*q2)*3.0-coef_f3_q_sym21*q0*(q3*q3)+coef_f0_q_sym6*q0*q1*q3+coef_f0_q_sym8*q0*q2*q3*2.0+coef_f0_q_sym15*q1*q2*q3*2.0-coef_f3_q_sym15*q0*q1*q2*2.0-coef_f3_q_sym16*q0*q1*q3-coef_f3_q_sym20*q0*q2*q3*2.0;
    symobj[2][3] = coef_f0_q_sym11*q0+coef_f0_q_sym18*q1+coef_f0_q_sym22*q2+coef_f0_q_sym24*q3*2.0-coef_f3_q_sym24*q0+coef_f0_q_sym1*(q0*q0*q0)+coef_f0_q_sym12*(q1*q1*q1)+coef_f0_q_sym19*(q2*q2*q2)+coef_f0_q_sym23*(q3*q3*q3)*4.0-coef_f3_q_sym4*(q0*q0*q0)+coef_f0_q_sym2*(q0*q0)*q1+coef_f0_q_sym3*(q0*q0)*q2+coef_f0_q_sym5*q0*(q1*q1)+coef_f0_q_sym4*(q0*q0)*q3*2.0+coef_f0_q_sym8*q0*(q2*q2)+coef_f0_q_sym10*q0*(q3*q3)*3.0+coef_f0_q_sym13*(q1*q1)*q2+coef_f0_q_sym14*(q1*q1)*q3*2.0+coef_f0_q_sym15*q1*(q2*q2)+coef_f0_q_sym17*q1*(q3*q3)*3.0+coef_f0_q_sym20*(q2*q2)*q3*2.0+coef_f0_q_sym21*q2*(q3*q3)*3.0-coef_f3_q_sym7*(q0*q0)*q1-coef_f3_q_sym9*(q0*q0)*q2-coef_f3_q_sym10*(q0*q0)*q3*2.0-coef_f3_q_sym14*q0*(q1*q1)-coef_f3_q_sym20*q0*(q2*q2)-coef_f3_q_sym23*q0*(q3*q3)*3.0+coef_f0_q_sym6*q0*q1*q2+coef_f0_q_sym7*q0*q1*q3*2.0+coef_f0_q_sym9*q0*q2*q3*2.0+coef_f0_q_sym16*q1*q2*q3*2.0-coef_f3_q_sym16*q0*q1*q2-coef_f3_q_sym17*q0*q1*q3*2.0-coef_f3_q_sym21*q0*q2*q3*2.0;
    symobj[3][0] = q0*2.0;
    symobj[3][1] = q1*2.0;
    symobj[3][2] = q2*2.0;
    symobj[3][3] = q3*2.0;

    for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 4; ++j)
            Jacob(i, j) = symobj[i][j];
}


void v_func_pTop(Eigen::Matrix<double, 28, 1>& v,
                 const Eigen::Vector4d& q)
{
    v(0, 0) = (q(0)*q(0)*q(0))*q(1);
    v(1, 0) = (q(0)*q(0)*q(0))*q(2);
    v(2, 0) = (q(0)*q(0)*q(0))*q(3);
    v(3, 0) = q(0)*(q(1)*q(1)*q(1));
    v(4, 0) = q(0)*(q(1)*q(1))*q(2);
    v(5, 0) = q(0)*(q(1)*q(1))*q(3);
    v(6, 0) = q(0)*q(1)*(q(2)*q(2));
    v(7, 0) = q(0)*q(1)*q(2)*q(3);
    v(8, 0) = q(0)*q(1)*(q(3)*q(3));
    v(9, 0) = q(0)*(q(2)*q(2)*q(2));
    v(10, 0) = q(0)*(q(2)*q(2))*q(3);
    v(11, 0) = q(0)*q(2)*(q(3)*q(3));
    v(12, 0) = q(0)*(q(3)*q(3)*q(3));
    v(13, 0) = q(1)*q(1)*q(1)*q(1);
    v(14, 0) = (q(1)*q(1)*q(1))*q(2);
    v(15, 0) = (q(1)*q(1)*q(1))*q(3);
    v(16, 0) = (q(1)*q(1))*(q(2)*q(2));
    v(17, 0) = (q(1)*q(1))*q(2)*q(3);
    v(18, 0) = (q(1)*q(1))*(q(3)*q(3));
    v(19, 0) = q(1)*(q(2)*q(2)*q(2));
    v(20, 0) = q(1)*(q(2)*q(2))*q(3);
    v(21, 0) = q(1)*q(2)*(q(3)*q(3));
    v(22, 0) = q(1)*(q(3)*q(3)*q(3));
    v(23, 0) = q(2)*q(2)*q(2)*q(2);
    v(24, 0) = (q(2)*q(2)*q(2))*q(3);
    v(25, 0) = (q(2)*q(2))*(q(3)*q(3));
    v(26, 0) = q(2)*(q(3)*q(3)*q(3));
    v(27, 0) = q(3)*q(3)*q(3)*q(3);
}



void y_func_pTop(Eigen::Matrix<double, 9, 1>& y,
                 const Eigen::Vector4d& q)
{
    y(0, 0) = q(0)*q(1);
    y(1, 0) = q(0)*q(2);
    y(2, 0) = q(0)*q(3);
    y(3, 0) = q(1)*q(1);
    y(4, 0) = q(1)*q(2);
    y(5, 0) = q(1)*q(3);
    y(6, 0) = q(2)*q(2);
    y(7, 0) = q(2)*q(3);
    y(8, 0) = q(3)*q(3);
}


void jv_func_pTop(Eigen::MatrixXd& jv,
                  const Eigen::Vector4d& q)
{
    jv(0, 0) = (q(0)*q(0))*q(1)*3.0;
    jv(0, 1) = q(0)*q(0)*q(0);
    jv(1, 0) = (q(0)*q(0))*q(2)*3.0;
    jv(1, 2) = q(0)*q(0)*q(0);
    jv(2, 0) = (q(0)*q(0))*q(3)*3.0;
    jv(2, 3) = q(0)*q(0)*q(0);
    jv(3, 0) = q(1)*q(1)*q(1);
    jv(3, 1) = q(0)*(q(1)*q(1))*3.0;
    jv(4, 0) = (q(1)*q(1))*q(2);
    jv(4, 1) = q(0)*q(1)*q(2)*2.0;
    jv(4, 2) = q(0)*(q(1)*q(1));
    jv(5, 0) = (q(1)*q(1))*q(3);
    jv(5, 1) = q(0)*q(1)*q(3)*2.0;
    jv(5, 3) = q(0)*(q(1)*q(1));
    jv(6, 0) = q(1)*(q(2)*q(2));
    jv(6, 1) = q(0)*(q(2)*q(2));
    jv(6, 2) = q(0)*q(1)*q(2)*2.0;
    jv(7, 0) = q(1)*q(2)*q(3);
    jv(7, 1) = q(0)*q(2)*q(3);
    jv(7, 2) = q(0)*q(1)*q(3);
    jv(7, 3) = q(0)*q(1)*q(2);
    jv(8, 0) = q(1)*(q(3)*q(3));
    jv(8, 1) = q(0)*(q(3)*q(3));
    jv(8, 3) = q(0)*q(1)*q(3)*2.0;
    jv(9, 0) = q(2)*q(2)*q(2);
    jv(9, 2) = q(0)*(q(2)*q(2))*3.0;
    jv(10, 0) = (q(2)*q(2))*q(3);
    jv(10, 2) = q(0)*q(2)*q(3)*2.0;
    jv(10, 3) = q(0)*(q(2)*q(2));
    jv(11, 0) = q(2)*(q(3)*q(3));
    jv(11, 2) = q(0)*(q(3)*q(3));
    jv(11, 3) = q(0)*q(2)*q(3)*2.0;
    jv(12, 0) = q(3)*q(3)*q(3);
    jv(12, 3) = q(0)*(q(3)*q(3))*3.0;
    jv(13, 1) = (q(1)*q(1)*q(1))*4.0;
    jv(14, 1) = (q(1)*q(1))*q(2)*3.0;
    jv(14, 2) = q(1)*q(1)*q(1);
    jv(15, 1) = (q(1)*q(1))*q(3)*3.0;
    jv(15, 3) = q(1)*q(1)*q(1);
    jv(16, 1) = q(1)*(q(2)*q(2))*2.0;
    jv(16, 2) = (q(1)*q(1))*q(2)*2.0;
    jv(17, 1) = q(1)*q(2)*q(3)*2.0;
    jv(17, 2) = (q(1)*q(1))*q(3);
    jv(17, 3) = (q(1)*q(1))*q(2);
    jv(18, 1) = q(1)*(q(3)*q(3))*2.0;
    jv(18, 3) = (q(1)*q(1))*q(3)*2.0;
    jv(19, 1) = q(2)*q(2)*q(2);
    jv(19, 2) = q(1)*(q(2)*q(2))*3.0;
    jv(20, 1) = (q(2)*q(2))*q(3);
    jv(20, 2) = q(1)*q(2)*q(3)*2.0;
    jv(20, 3) = q(1)*(q(2)*q(2));
    jv(21, 1) = q(2)*(q(3)*q(3));
    jv(21, 2) = q(1)*(q(3)*q(3));
    jv(21, 3) = q(1)*q(2)*q(3)*2.0;
    jv(22, 1) = q(3)*q(3)*q(3);
    jv(22, 3) = q(1)*(q(3)*q(3))*3.0;
    jv(23, 2) = (q(2)*q(2)*q(2))*4.0;
    jv(24, 2) = (q(2)*q(2))*q(3)*3.0;
    jv(24, 3) = q(2)*q(2)*q(2);
    jv(25, 2) = q(2)*(q(3)*q(3))*2.0;
    jv(25, 3) = (q(2)*q(2))*q(3)*2.0;
    jv(26, 2) = q(3)*q(3)*q(3);
    jv(26, 3) = q(2)*(q(3)*q(3))*3.0;
    jv(27, 3) = (q(3)*q(3)*q(3))*4.0;
}



void jy_func_pTop(Eigen::MatrixXd& jy,
                  const Eigen::Vector4d& q)
{
    jy(0, 0) = q(1);
    jy(0, 1) = q(0);
    jy(1, 0) = q(2);
    jy(1, 2) = q(0);
    jy(2, 0) = q(3);
    jy(2, 3) = q(0);
    jy(3, 1) = q(1)*2.0;
    jy(4, 1) = q(2);
    jy(4, 2) = q(1);
    jy(5, 1) = q(3);
    jy(5, 3) = q(1);
    jy(6, 2) = q(2)*2.0;
    jy(7, 2) = q(3);
    jy(7, 3) = q(2);
    jy(8, 3) = q(3)*2.0;
}
