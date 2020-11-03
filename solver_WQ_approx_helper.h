//
// Created by Jin Wu on 1/11/2020.
//

#ifndef LIBQPEP_SOLVER_WQ_APPROX_HELPER_H
#define LIBQPEP_SOLVER_WQ_APPROX_HELPER_H

#include <time.h>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry> // For Quaternion


void data_func_WQ_approx_sparse(Eigen::SparseMatrix<double>& tmp,
                                          Eigen::MatrixXd& tmp2,
                                          const Eigen::VectorXd& data);


#endif //LIBQPEP_SOLVER_WQ_APPROX_HELPER_H
