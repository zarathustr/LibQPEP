#ifndef LIBQPEP_GENERATEPROJECTEDPOINTS_H
#define LIBQPEP_GENERATEPROJECTEDPOINTS_H

#include <time.h>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Geometry> // For Quaternion

void generateProjectedPoints(std::vector<Eigen::Vector2d>& image_pt,
                             std::vector<double>& s,
                             const std::vector<Eigen::Vector3d>& world_pt,
                             const Eigen::Matrix3d& K,
                             const Eigen::Matrix3d& R,
                             const Eigen::Vector3d& t);

#endif //LIBQPEP_GENERATEPROJECTEDPOINTS_H
