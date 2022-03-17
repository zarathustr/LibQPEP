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
//
//
// generateProjectedPoints.cpp: Generate image points for PnP problems


#include "generateProjectedPoints.h"

void generateProjectedPoints(std::vector<Eigen::Vector2d>& image_pt,
                             std::vector<double>& s,
                             const std::vector<Eigen::Vector3d>& world_pt,
                             const Eigen::Matrix3d& K,
                             const Eigen::Matrix3d& R,
                             const Eigen::Vector3d& t)
{
    int numPoints = world_pt.size();

    for(int i = 0; i < numPoints; ++i)
    {
        s.push_back(0);
    }

    for(int i = 0; i < numPoints; ++i)
    {
        Eigen::Vector4d world_point(world_pt[i](0), world_pt[i](1), world_pt[i](2), 1);
        Eigen::Matrix<double, 3, 4> cameraMatrix;
        cameraMatrix << R.transpose(), t;
        cameraMatrix = K * cameraMatrix;
        Eigen::Vector3d projectedPoint = cameraMatrix * world_point;
        s[i] = projectedPoint(2);
        Eigen::Vector2d projectedPoint_(projectedPoint(0), projectedPoint(1));
        image_pt.push_back(projectedPoint_ / s[i]);
    }
}
