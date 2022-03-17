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

#ifndef LIBQPEP_PNP_WQD_H
#define LIBQPEP_PNP_WQD_H

#include <time.h>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Geometry> // For Quaternion


void pnp_WQD(Eigen::Matrix<double, 4, 64>& W,
             Eigen::Matrix<double, 4, 4>& Q,
             Eigen::Matrix<double, 3, 37>& D,
             Eigen::Matrix<double, 4, 24>& coef_f_q_sym,
             Eigen::Matrix<double, 1, 70>& coef_J_pure,
             Eigen::Matrix<double, 3, 10>& coefs_tq,
             Eigen::Matrix<double, 3, 3>& pinvG,
             const std::vector<Eigen::Vector2d>& image_pt,
             const std::vector<Eigen::Vector3d>& world_pt,
             const Eigen::Matrix3d& K,
             const double& scale);

#endif //LIBQPEP_PNP_WQD_H
