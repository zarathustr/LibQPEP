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

#ifndef LIBQPEP_QPEP_LM_SINGLE_H
#define LIBQPEP_QPEP_LM_SINGLE_H

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
#include "QPEP_grobner.h"

struct QPEP_runtime QPEP_lm_single(Eigen::Matrix3d& R,
                    Eigen::Vector3d& t,
                    Eigen::Matrix4d& X,
                    const Eigen::Vector4d& q0,
                    const int& max_iter,
                    const double& mu,
                    const eq_Jacob_func_handle& eq_Jacob_func,
                    const t_func_handle& t_func,
                    const Eigen::MatrixXd& coef_f_q_sym,
                    const Eigen::MatrixXd& coefs_tq,
                    const Eigen::MatrixXd& pinvG,
                    const struct QPEP_runtime& stat_);

struct QPEP_runtime QPEP_lm_fsolve(Eigen::Matrix3d& R,
                                   Eigen::Vector3d& t,
                                   Eigen::Matrix4d& X,
                                   const Eigen::Vector4d& q0,
                                   const int& max_iter,
                                   const double& mu,
                                   const eq_Jacob_func_handle& eq_Jacob_func,
                                   const t_func_handle& t_func,
                                   const Eigen::MatrixXd& coef_f_q_sym,
                                   const Eigen::MatrixXd& coefs_tq,
                                   const Eigen::MatrixXd& pinvG,
                                   const Eigen::MatrixXd& W,
                                   const Eigen::MatrixXd& Q,
                                   const struct QPEP_runtime& stat_);

#endif //LIBQPEP_QPEP_LM_SINGLE_H
