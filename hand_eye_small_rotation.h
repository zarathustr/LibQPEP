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

#ifndef LIBQPEP_HAND_EYE_SMALL_ROTATION_H
#define LIBQPEP_HAND_EYE_SMALL_ROTATION_H

#include <time.h>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Geometry> // For Quaternion

void hand_eye_small_rotation(Eigen::Matrix4d& X,
                             const std::vector<Eigen::Matrix4d>& A,
                             const std::vector<Eigen::Matrix4d>& B);

void hand_eye_small_rotation_refine(Eigen::Matrix4d& X,
                                    const std::vector<Eigen::Matrix4d>& A,
                                    const std::vector<Eigen::Matrix4d>& B,
                                    const Eigen::Matrix4d& X0,
                                    const int& num);

#endif