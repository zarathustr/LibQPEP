//
// Created by Jin Wu on 1/11/2020.
//

#ifndef LIBQPEP_MISC_PTOP_FUNCS_H
#define LIBQPEP_MISC_PTOP_FUNCS_H

#include <time.h>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Geometry> // For Quaternion


Eigen::Vector3d t_pTop_func(const Eigen::MatrixXd& pinvG, const Eigen::MatrixXd& coefs_tq, const Eigen::Quaterniond& q);
std::vector<double> mon_J_pure_pTop_func(const Eigen::Quaterniond& q, const Eigen::Vector3d& t);
void eq_Jacob_pTop_func(Eigen::Matrix<double, 4, 1>& eq,
                        Eigen::Matrix<double, 4, 4>& Jacob,
                        const Eigen::MatrixXd& coef_f_q_sym,
                        const Eigen::Vector4d& q);

#endif //LIBQPEP_MISC_PTOP_FUNCS_H
