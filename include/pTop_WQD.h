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

#ifndef LIBQPEP_PTOP_WQD_H
#define LIBQPEP_PTOP_WQD_H

#include <time.h>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Geometry> // For Quaternion

void pTop_WQD(Eigen::Matrix<double, 4, 64>& W,
              Eigen::Matrix<double, 4, 4>& Q,
              Eigen::Matrix<double, 3, 28>& D,
              Eigen::Matrix<double, 3, 9>& G,
              Eigen::Vector3d& c,
              Eigen::Matrix<double, 4, 24>& coef_f_q_sym,
              Eigen::Matrix<double, 1, 85>& coef_J_pure,
              Eigen::Matrix<double, 3, 11>& coefs_tq,
              Eigen::Matrix<double, 3, 3>& pinvG,
              const std::vector<Eigen::Vector3d>& rr,
              const std::vector<Eigen::Vector3d>& bb,
              const std::vector<Eigen::Vector3d>& nv);

#endif //LIBQPEP_PTOP_WQD_H
