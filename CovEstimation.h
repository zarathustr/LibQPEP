//
// Created by Jin Wu on 6/11/2020.
//

#ifndef LIBQPEP_COVESTIMATION_H
#define LIBQPEP_COVESTIMATION_H


#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Geometry> // For Quaternion
#include <Eigen/Sparse>



void csdp_cov(Eigen::Matrix4d& cov,
              const Eigen::MatrixXd& F,
              const Eigen::Matrix3d& cov_left,
              const Eigen::Vector4d& q);

#endif //LIBQPEP_COVESTIMATION_H
