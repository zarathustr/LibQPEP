//
// Created by Jin Wu on 3/11/2020.
//

#ifndef LIBQPEP_SOLVER_HELPER_H
#define LIBQPEP_SOLVER_HELPER_H


#include <time.h>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Geometry> // For Quaternion
#include <Eigen/Sparse>
#include "QPEP.h"

void data_func_solver_WQ_func(Eigen::SparseMatrix<double>& C1,
                              Eigen::MatrixXd& C2,
                              const Eigen::VectorXd& data);

#endif //LIBQPEP_SOLVER_HELPER_H
