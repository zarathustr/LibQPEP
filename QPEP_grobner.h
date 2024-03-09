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
#ifndef LIBQPEP_QPEP_GROBNER_H
#define LIBQPEP_QPEP_GROBNER_H

#include <time.h>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Geometry> // For Quaternion
#include "QPEP.h"

typedef struct QPEP_runtime (* solver_func_handle)(Eigen::MatrixXcd& sol, const Eigen::VectorXd& data, const struct QPEP_options& opt);
typedef std::vector<double> (* mon_J_pure_func_handle)(const Eigen::Quaterniond& q, const Eigen::Vector3d& t);
typedef Eigen::Vector3d (* t_func_handle)(const Eigen::MatrixXd& pinvG, const Eigen::MatrixXd& coefs_tq, const Eigen::Quaterniond& q);
typedef void (* eq_Jacob_func_handle)(Eigen::Matrix<double, 4, 1>& eq, Eigen::Matrix<double, 4, 4>& Jacob, const Eigen::MatrixXd& coef_f_q_sym, const Eigen::Vector4d& q);

struct QPEP_runtime QPEP_WQ_grobner(Eigen::Matrix3d& R,
                                    Eigen::Vector3d& t,
                                    Eigen::Matrix4d& X,
                                    double* min,
                                    const Eigen::MatrixXd& W,
                                    const Eigen::MatrixXd& Q,
                                    const solver_func_handle& solver_func,
                                    const mon_J_pure_func_handle& mon_J_pure_func,
                                    const t_func_handle& t_func,
                                    const Eigen::MatrixXd& coef_J_pure,
                                    const Eigen::MatrixXd& coefs_tq,
                                    const Eigen::MatrixXd& pinvG,
                                    const int* perm,
                                    const struct QPEP_options& opt);

#endif //LIBQPEP_QPEP_GROBNER_H
