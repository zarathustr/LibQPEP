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

#ifndef LIBQPEP_HAND_EYE_WQD_H
#define LIBQPEP_HAND_EYE_WQD_H

#include <time.h>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#ifndef NO_OMP
#include "omp.h"
#endif

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Geometry> // For Quaternion

void hand_eye_WQD(Eigen::Matrix<double, 4, 64>& W,
              Eigen::Matrix<double, 4, 4>& Q,
              Eigen::Matrix<double, 3, 28>& D,
              Eigen::Matrix<double, 3, 9>& G,
              Eigen::Vector3d& c,
              Eigen::Matrix<double, 4, 24>& coef_f_q_sym,
              Eigen::Matrix<double, 1, 85>& coef_J_pure,
              Eigen::Matrix<double, 3, 11>& coefs_tq,
              Eigen::Matrix<double, 3, 3>& pinvG,
              const std::vector<Eigen::Matrix4d>& As,
              const std::vector<Eigen::Matrix4d>& Bs);

#endif //LIBQPEP_HAND_EYE_WQD_H
