//
// Created by Jin Wu on 1/11/2020.
//

#ifndef LIBQPEP_QPEP_LM_SINGLE_H
#define LIBQPEP_QPEP_LM_SINGLE_H

#include <time.h>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>

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
